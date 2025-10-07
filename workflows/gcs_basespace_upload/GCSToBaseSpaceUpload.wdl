version 1.0

workflow GCSToBaseSpaceUpload {
    input {
        String gs_bucket_path
        String basespace_project_name
        String? default_project_id
        Float max_batch_size_gb = 50.0
        Int max_concurrent_uploads_per_vm = 10
        Int chunk_size_mb = 25
        Int check_interval_seconds = 60
        String? basespace_access_token
        String docker_image = "google/cloud-sdk:latest"
        Boolean upload_r1_only = false
    }

    # First create basespace config file if token is provided
    if (defined(basespace_access_token)) {
        call CreateConfigFile {
            input:
                basespace_access_token = select_first([basespace_access_token, ""]),
                docker_image = docker_image
        }
    }

    # Check if BaseSpace project exists, if not create it
    call CheckOrCreateBaseSpaceProject {
        input:
            project_name = basespace_project_name,
            basespace_access_token = basespace_access_token,
            docker_image = docker_image,
            config_file = select_first([CreateConfigFile.config_file, ""])
    }

    # Determine which project ID to use - use default if returned value is empty
    String effective_project_id = if (CheckOrCreateBaseSpaceProject.project_id == "") then select_first([default_project_id, ""]) else CheckOrCreateBaseSpaceProject.project_id

    # List all files in the GCS bucket with their sizes and paths
    call ListGCSFilesWithSizes {
        input:
            gs_bucket_path = gs_bucket_path,
            docker_image = docker_image
    }

    # Group files into batches based on size
    call CreateFileBatches {
        input:
            file_manifest = ListGCSFilesWithSizes.file_manifest,
            max_batch_size_gb = max_batch_size_gb,
            docker_image = docker_image,
            upload_r1_only = upload_r1_only
    }

    # Process each batch in parallel
    scatter (batch_file in CreateFileBatches.batch_files) {
        call ProcessBatch {
            input:
                batch_file = batch_file,
                gs_bucket_path = gs_bucket_path,
                project_id = effective_project_id,
                max_concurrent_uploads = max_concurrent_uploads_per_vm,
                chunk_size_mb = chunk_size_mb,
                check_interval_seconds = check_interval_seconds,
                basespace_access_token = basespace_access_token,
                docker_image = docker_image,
                config_file = select_first([CreateConfigFile.config_file, ""])
        }
    }

    # Combine all results to verify completion
    call VerifyAllUploads {
        input:
            batch_results = ProcessBatch.batch_result,
            project_id = effective_project_id,
            basespace_access_token = basespace_access_token,
            docker_image = docker_image,
            config_file = select_first([CreateConfigFile.config_file, ""])
    }

    output {
        String basespace_project_id = effective_project_id
        String basespace_project_url = CheckOrCreateBaseSpaceProject.project_url
        Array[File] batch_upload_logs = ProcessBatch.upload_log
        String upload_summary = VerifyAllUploads.summary
        File detailed_results = VerifyAllUploads.detailed_results
        File orphaned_r2_files = CreateFileBatches.orphaned_r2_files
        File orphaned_r2_summary = CreateFileBatches.orphaned_r2_summary
    }
}

task CreateConfigFile {
    input {
        String basespace_access_token
        String docker_image
    }

    command <<<
        # Create BaseSpace config directory
        mkdir -p ~/.basespace

        # Create config file
        echo "apiServer   = https://api.basespace.illumina.com" > default.cfg
        echo "accessToken = ~{basespace_access_token}" >> default.cfg
    >>>

    output {
        File config_file = "default.cfg"
    }

    runtime {
        docker: docker_image
    }
}

task CheckOrCreateBaseSpaceProject {
    input {
        String project_name
        String? basespace_access_token
        String docker_image
        String config_file = ""
    }

    command <<<
        # Install necessary tools
        apt-get update && apt-get install -y wget curl jq unzip

        # Install BaseSpace CLI
        wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs"
        chmod +x bs
        mv bs /usr/local/bin/

        # Set up authentication if token is provided
        if [ ! -z "~{basespace_access_token}" ]; then
            # Create BaseSpace config directory
            mkdir -p ~/.basespace
            
            # Copy config file if provided
            if [ ! -z "~{config_file}" ]; then
                gsutil cp ~{config_file} ~/.basespace/default.cfg
            fi
            
            # Authenticated with BaseSpace
            echo "Authenticatied with BaseSpace..."
        else
            echo "Using existing BaseSpace authentication..."
        fi

        # Check if project with the same name already exists
        echo "Checking if BaseSpace project already exists: ~{project_name}"
        EXISTING_PROJECTS=$(bs list projects --format=json 2>/dev/null || echo '[]')

        # Validate JSON before processing
        if ! echo "$EXISTING_PROJECTS" | jq empty 2>/dev/null; then
            echo "Warning: Invalid JSON response when listing projects. Proceeding with creating a new project."
            EXISTING_PROJECTS='[]'
        fi

        PROJECT_EXISTS=$(echo "$EXISTING_PROJECTS" | jq -r --arg name "~{project_name}" '.[] | select(.Name == $name)' 2>/dev/null || echo '')


        if [ -n "$PROJECT_EXISTS" ]; then
            PROJECT_ID=$(echo "$PROJECT_EXISTS" | jq -r '.Id')
            PROJECT_URL=$(echo "$PROJECT_EXISTS" | jq -r '.HrefBaseSpaceUI')
            echo "Project with name '~{project_name}' already exists with ID: $PROJECT_ID"
            echo "Using existing project instead of creating a new one"
        else
            # Create project
            echo "Creating BaseSpace project: ~{project_name}"
            RESPONSE=$(bs create project --name="~{project_name}" --format=json)
            PROJECT_ID=$(echo $RESPONSE | jq -r '.Id')
            PROJECT_URL=$(echo $RESPONSE | jq -r '.HrefBaseSpaceUI')
            echo "Created project with ID: $PROJECT_ID"
        fi

        # Handle potential null or empty project ID
        if [ -z "$PROJECT_ID" ]; then
            echo "Failed to get or create project ID"
            echo "" > project_id.txt
        else
            echo $PROJECT_ID > project_id.txt
        fi
        
        # Handle potential null or empty project URL
        if [ -z "$PROJECT_URL" ]; then
            echo "" > project_url.txt
        else
            echo $PROJECT_URL > project_url.txt
        fi
    >>>

    output {
        String project_id = read_string("project_id.txt")
        String project_url = read_string("project_url.txt")
    }

    runtime {
        docker: docker_image
        cpu: 2
        memory: "4 GB"
        disks: "local-disk 10 SSD"
    }
}

task ListGCSFilesWithSizes {
    input {
        String gs_bucket_path
        String docker_image
    }

    command <<<
        # Install required tools
        apt-get update && apt-get install -y jq bc

        # Ensure the bucket path ends with a slash for directory listing
        if [[ "~{gs_bucket_path}" != */ ]]; then
            GS_PATH="~{gs_bucket_path}/"
        else
            GS_PATH="~{gs_bucket_path}"
        fi

        # Get stats for all files in the bucket path
        echo "Listing files in: $GS_PATH"
        gsutil -m ls -l "${GS_PATH}**" | grep -v '/$' > all_files_raw.txt

        # Create a CSV file with file details
        echo "full_path,relative_path,size_bytes,size_gb" > files.csv

        # Process each line in the file listing
        while IFS= read -r line; do
            # Skip directory markers and summary lines
            if [[ "$line" == *"TOTAL:"* ]] || [[ "$line" == *"gs://"/$ ]]; then
                continue
            fi
            
            # The format of gsutil ls -l is typically:
            # <size>  <date/time>  <gs-url>
            
            # Extract the GCS path - it's the last part of the line after the timestamp
            if [[ $line =~ (gs://.+)$ ]]; then
                full_path="${BASH_REMATCH[1]}"
                # Extract size - it's the first number in the line
                if [[ $line =~ ^[[:space:]]*([0-9]+) ]]; then
                    size_bytes="${BASH_REMATCH[1]}"
                    
                    # Skip if not a valid path
                    if [[ -z "$full_path" ]]; then
                        continue
                    fi
                    
                    # Clean up any leading/trailing whitespace
                    full_path=$(echo "$full_path" | xargs)
                    
                    # Calculate relative path by removing bucket path prefix
                    relative_path=${full_path#"$GS_PATH"}
                    
                    # Calculate size in GB with 9 decimal places for precision
                    size_gb=$(echo "scale=9; $size_bytes / 1073741824" | bc)
                    
                    # Add to CSV
                    echo "$full_path,$relative_path,$size_bytes,$size_gb" >> files.csv
                fi
            fi
        done < all_files_raw.txt

        # Convert CSV to JSON manifest
        echo "[" > file_manifest.json
        first_line=true

        # Skip header line and process each file entry
        tail -n +2 files.csv | while IFS=, read -r full_path relative_path size_bytes size_gb; do
            # Add comma separator for all but the first entry
            if [ "$first_line" = true ]; then
                first_line=false
            else
                echo "," >> file_manifest.json
            fi

            # Create JSON entry for this file
            cat <<EOF >> file_manifest.json
  {
    "full_path": "$full_path",
    "relative_path": "$relative_path",
    "size_bytes": $size_bytes,
    "size_gb": $size_gb
  }
EOF
        done

        echo "]" >> file_manifest.json

        # Print summary
        total_files=$(grep -c "full_path" file_manifest.json)
        total_bytes=$(awk 'BEGIN {total=0} /"size_bytes":/ {gsub(/[^0-9]/, "", $2); total += $2} END {print total}' file_manifest.json)
        total_gb=$(echo "scale=2; $total_bytes / 1073741824" | bc)
        echo "Found $total_files files totaling $total_gb GB"
    >>>

    output {
        File file_manifest = "file_manifest.json"
    }

    runtime {
        docker: docker_image
        cpu: 2
        memory: "4 GB"
        disks: "local-disk 10 SSD"
    }
}


task CreateFileBatches {
    input {
        File file_manifest
        Float max_batch_size_gb
        String docker_image
        Boolean upload_r1_only = false
    }

    command <<<
        # Install required tools
        apt-get update && apt-get install -y jq bc

        # Convert max_batch_size_gb to bytes for comparison (ensure integer result)
        max_batch_size_bytes=$(echo "~{max_batch_size_gb} * 1073741824 / 1" | bc | cut -d'.' -f1)
        echo "Maximum batch size: ~{max_batch_size_gb} GB ($max_batch_size_bytes bytes)"

        # Create a working directory for batch files
        mkdir -p batches

        # Group paired-end files and identify orphaned R2 files
        echo "Analyzing paired-end fastq files..."
        
        # Extract base names for R1 files
        jq -r '.[] | select(.relative_path | endswith("_R1_001.fastq.gz")) | 
            (.relative_path | gsub("_R1_001\\.fastq\\.gz$"; ""))' ~{file_manifest} | 
            sort | uniq > r1_base_names.txt
            
        # Extract base names for R2 files
        jq -r '.[] | select(.relative_path | endswith("_R2_001.fastq.gz")) | 
            (.relative_path | gsub("_R2_001\\.fastq\\.gz$"; ""))' ~{file_manifest} | 
            sort | uniq > r2_base_names.txt
            
        echo "Found $(wc -l < r1_base_names.txt) R1 files"
        echo "Found $(wc -l < r2_base_names.txt) R2 files"

        # Find orphaned R2 files (R2 files without corresponding R1)
        echo "Identifying orphaned R2 files..."
        comm -23 r2_base_names.txt r1_base_names.txt > orphaned_r2_base_names.txt
        orphaned_r2_count=$(wc -l < orphaned_r2_base_names.txt)
        
        if [ $orphaned_r2_count -gt 0 ]; then
            echo "WARNING: Found $orphaned_r2_count orphaned R2 files (R2 files without corresponding R1 files)"
            
            # Create detailed list of orphaned R2 files
            echo "[]" > orphaned_r2_files.json
            while read -r base_name; do
                # Find the R2 file for this base
                r2_file=$(jq -r '.[] | select(.relative_path == "'"$base_name"'_R2_001.fastq.gz")' ~{file_manifest})
                
                if [[ -n "$r2_file" && "$r2_file" != "null" ]]; then
                    # Add to orphaned R2 files list
                    orphaned_r2_files_temp=$(cat orphaned_r2_files.json)
                    if [ "$orphaned_r2_files_temp" = "[]" ]; then
                        echo "[$r2_file]" > orphaned_r2_files.json
                    else
                        # Remove closing bracket, add comma and new entry, then close bracket
                        sed -i '$ s/]$/,/' orphaned_r2_files.json
                        echo "$r2_file]" >> orphaned_r2_files.json
                    fi
                    echo "Found orphaned R2 file: ${base_name}_R2_001.fastq.gz"
                fi
            done < orphaned_r2_base_names.txt
            
            # Create human-readable summary of orphaned R2 files
            echo "Orphaned R2 Files (files without corresponding R1 files):" > orphaned_r2_summary.txt
            echo "=========================================================" >> orphaned_r2_summary.txt
            echo "Total orphaned R2 files: $orphaned_r2_count" >> orphaned_r2_summary.txt
            echo "" >> orphaned_r2_summary.txt
            echo "File list:" >> orphaned_r2_summary.txt
            jq -r '.[] | "- " + .relative_path + " (" + (.size_gb | tostring) + " GB)"' orphaned_r2_files.json >> orphaned_r2_summary.txt
            echo "" >> orphaned_r2_summary.txt
            echo "These files have been excluded from batches and require manual review." >> orphaned_r2_summary.txt
            
        else
            echo "No orphaned R2 files found"
            echo "[]" > orphaned_r2_files.json
            echo "No orphaned R2 files found." > orphaned_r2_summary.txt
        fi

        # Filter out R2 files if upload_r1_only is true
        if [ "~{upload_r1_only}" = "true" ]; then
            echo "R1-only mode enabled: excluding all R2 files from upload"
            
            # Get all R2 file paths
            jq -r '.[] | select(.relative_path | endswith("_R2_001.fastq.gz")) | .full_path' ~{file_manifest} > r2_files_to_exclude.txt
            
            # Add R2 files to the excluded files list
            cat r2_files_to_exclude.txt >> excluded_files_list.txt
            
            # Update paired_files_map to only include R1 files
            echo "{" > paired_files_map_r1_only.json
            first_entry=true
            
            while read -r base_name; do
                # Find only R1 file for this base
                r1_file=$(jq -r '.[] | select(.relative_path == "'"$base_name"'_R1_001.fastq.gz")' ~{file_manifest})
                
                if [[ -n "$r1_file" && "$r1_file" != "null" ]]; then
                    # Add to map (single file instead of pair)
                    if [ "$first_entry" = true ]; then
                        first_entry=false
                    else
                        echo "," >> paired_files_map_r1_only.json
                    fi
                    
                    echo "  \"$base_name\": [" >> paired_files_map_r1_only.json
                    echo "    $r1_file" >> paired_files_map_r1_only.json
                    echo "  ]" >> paired_files_map_r1_only.json
                    
                    echo "Added R1-only for: $base_name"
                fi
            done < paired_base_names.txt
            
            echo "}" >> paired_files_map_r1_only.json
            
            # Replace the paired_files_map with R1-only version
            mv paired_files_map_r1_only.json paired_files_map.json
            
            echo "Excluded $(wc -l < r2_files_to_exclude.txt) R2 files from upload"
        fi

        # Find properly paired files (both R1 and R2 exist)
        echo "Identifying properly paired files..."
        comm -12 r1_base_names.txt r2_base_names.txt > paired_base_names.txt
        paired_count=$(wc -l < paired_base_names.txt)
        echo "Found $paired_count complete paired file sets"

        # Create paired files map (only for complete pairs)
        echo "{" > paired_files_map.json
        first_entry=true
        
        while read -r base_name; do
            # Find R1 and R2 files for this base
            r1_file=$(jq -r '.[] | select(.relative_path == "'"$base_name"'_R1_001.fastq.gz")' ~{file_manifest})
            r2_file=$(jq -r '.[] | select(.relative_path == "'"$base_name"'_R2_001.fastq.gz")' ~{file_manifest})
            
            # Skip if either file is missing (shouldn't happen with comm -12, but safety check)
            if [[ -z "$r1_file" || -z "$r2_file" || "$r1_file" == "null" || "$r2_file" == "null" ]]; then
                echo "Warning: Skipping incomplete pair for base: $base_name"
                continue
            fi
            
            # Add to map
            if [ "$first_entry" = true ]; then
                first_entry=false
            else
                echo "," >> paired_files_map.json
            fi
            
            echo "  \"$base_name\": [" >> paired_files_map.json
            echo "    $r1_file," >> paired_files_map.json
            echo "    $r2_file" >> paired_files_map.json
            echo "  ]" >> paired_files_map.json
            
            echo "Added complete pair for: $base_name"
        done < paired_base_names.txt
        
        echo "}" >> paired_files_map.json

        # Create list of all paired file paths and orphaned R2 file paths
        echo "Creating exclusion list..."
        jq -r 'if . == {} then [] else to_entries[] | .value[] | .full_path end' paired_files_map.json > paired_files_list.txt
        jq -r '.[] | .full_path' orphaned_r2_files.json > orphaned_r2_files_list.txt
        
        # Combine both lists for exclusion
        cat paired_files_list.txt orphaned_r2_files_list.txt > excluded_files_list.txt

        # Process non-paired files (excluding orphaned R2 files)
        echo "Processing non-paired files (excluding orphaned R2 files)..."
        
        # Initialize non_paired_files.json
        echo "[]" > non_paired_files.json
        
        # Get all files that are not in the excluded files list
        jq -r '.[] | @json' ~{file_manifest} | while read -r file_json; do
            file_path=$(echo "$file_json" | jq -r '.full_path')
            
            # Check if file is in excluded files list
            if ! grep -q "^$file_path$" excluded_files_list.txt 2>/dev/null; then
                # Add to non-paired files
                non_paired_files_temp=$(cat non_paired_files.json)
                if [ "$non_paired_files_temp" = "[]" ]; then
                    echo "[$file_json]" > non_paired_files.json
                else
                    # Remove closing bracket, add comma and new entry, then close bracket
                    sed -i '$ s/]$/,/' non_paired_files.json
                    echo "$file_json]" >> non_paired_files.json
                fi
                echo "Added non-paired file: $file_path"
            fi
        done
        
        # Sort non-paired files by size (largest first) for better packing
        if [ -s non_paired_files.json ] && [ "$(jq '. | length' non_paired_files.json)" -gt 0 ]; then
            echo "Sorting $(jq '. | length' non_paired_files.json) non-paired files by size..."
            jq 'sort_by(.size_bytes) | reverse' non_paired_files.json > sorted_non_paired.json
        else
            echo "No non-paired files found, creating empty sorted file"
            echo "[]" > sorted_non_paired.json
        fi

        # Initialize batch variables
        batch_num=0
        current_batch_file="batches/batch_$(printf "%04d" $batch_num).json"
        current_batch_size=0
        batch_count=0

        # Start the first batch
        echo "[" > "$current_batch_file"
        first_in_batch=true

        # First, handle paired-end files which must stay together
        if [ -s paired_files_map.json ] && [ "$(jq 'keys | length' paired_files_map.json)" -gt 0 ]; then
            echo "Processing paired-end files..."
            
            while read -r pair_entry; do
                base_name=$(echo "$pair_entry" | jq -r '.key')
                
                # Get paired files info
                r1_file=$(echo "$pair_entry" | jq -r '.value[0]')
                r2_file=$(echo "$pair_entry" | jq -r '.value[1]')
                
                # Check if this is R1-only mode (array will have only 1 element)
                num_files=$(echo "$pair_entry" | jq -r '.value | length')
                if [ "$num_files" -eq 1 ]; then
                    # R1-only mode
                    pair_total_size=$r1_size
                else
                    # Normal paired mode
                    r2_file=$(echo "$pair_entry" | jq -r '.value[1]')
                    r2_size=$(echo "$r2_file" | jq -r '.size_bytes')
                    pair_total_size=$((r1_size + r2_size))
                fi
                
                # If current batch isn't empty and adding this pair would exceed max size,
                # close current batch and start a new one
                pair_check=$((current_batch_size + pair_total_size))
                if [ "$first_in_batch" = false ] && [ "$pair_check" -gt "$max_batch_size_bytes" ]; then
                    echo "]" >> "$current_batch_file"
                    batch_count=$((batch_count + 1))
                    
                    # Start a new batch
                    batch_num=$((batch_num + 1))
                    current_batch_file="batches/batch_$(printf "%04d" $batch_num).json"
                    echo "[" > "$current_batch_file"
                    first_in_batch=true
                    current_batch_size=0
                fi
                
                # Add files to the current batch
                if [ "$first_in_batch" = false ]; then
                    echo "," >> "$current_batch_file"
                else
                    first_in_batch=false
                fi

                # Add R1 file
                echo "$r1_file" >> "$current_batch_file"

                # Add R2 file only if it exists (not in R1-only mode)
                if [ "$num_files" -gt 1 ]; then
                    echo "," >> "$current_batch_file"
                    echo "$r2_file" >> "$current_batch_file"
                fi

                current_batch_size=$((current_batch_size + pair_total_size))

                if [ "$num_files" -eq 1 ]; then
                    echo "Added R1-only file $base_name (${pair_total_size} bytes) to batch $batch_num"
                else
                    echo "Added paired files $base_name (${pair_total_size} bytes) to batch $batch_num"
                fi
            done < <(jq -r 'to_entries[] | @json' paired_files_map.json)
        else
            echo "No paired files found to process"
        fi

        # Now process individual non-paired files (orphaned R2 files already excluded)
        echo "Processing non-paired files..."
        if [ -s sorted_non_paired.json ] && [ "$(jq '. | length' sorted_non_paired.json)" -gt 0 ]; then
            # Use process substitution to avoid subshell issues
            while read -r file_json; do
                # Skip empty lines
                if [ -z "$file_json" ]; then
                    continue
                fi
                
                # Extract file size and path
                file_size=$(echo "$file_json" | jq -r '.size_bytes')
                file_path=$(echo "$file_json" | jq -r '.full_path')
                
                echo "Processing file: $file_path (size: $file_size bytes)"
                
                # If this single file is larger than max batch size, it gets its own batch
                if [ "$file_size" -gt "$max_batch_size_bytes" ]; then
                    echo "File exceeds max batch size, creating dedicated batch"
                    
                    # If we've already started a batch with content, close it
                    if [ "$first_in_batch" = false ]; then
                        echo "]" >> "$current_batch_file"
                        batch_count=$((batch_count + 1))
                        echo "Closed batch $batch_num with size ${current_batch_size} bytes"
                        
                        # Start a new batch for this single large file
                        batch_num=$((batch_num + 1))
                        current_batch_file="batches/batch_$(printf "%04d" $batch_num).json"
                        echo "[" > "$current_batch_file"
                        first_in_batch=true
                        current_batch_size=0
                    fi
                    
                    # Add this file to its own batch
                    echo "$file_json" >> "$current_batch_file"
                    current_batch_size=$((current_batch_size + file_size))
                    first_in_batch=false
                    
                    # Close this batch and start a new one
                    echo "]" >> "$current_batch_file"
                    batch_count=$((batch_count + 1))
                    echo "Created dedicated batch $batch_num for large file ${file_path}"
                    
                    batch_num=$((batch_num + 1))
                    current_batch_file="batches/batch_$(printf "%04d" $batch_num).json"
                    echo "[" > "$current_batch_file"
                    first_in_batch=true
                    current_batch_size=0
                    
                    continue
                fi
                
                # If adding this file would exceed the batch size, start a new batch
                total_size_check=$((current_batch_size + file_size))
                if [ "$total_size_check" -gt "$max_batch_size_bytes" ] && [ "$first_in_batch" = false ]; then
                    echo "Adding file would exceed batch size, creating new batch"
                    echo "]" >> "$current_batch_file"
                    batch_count=$((batch_count + 1))
                    echo "Closed batch $batch_num with size ${current_batch_size} bytes"
                    
                    # Start a new batch
                    batch_num=$((batch_num + 1))
                    current_batch_file="batches/batch_$(printf "%04d" $batch_num).json"
                    echo "[" > "$current_batch_file"
                    first_in_batch=true
                    current_batch_size=0
                fi
                
                # Add comma separator if not the first file in the batch
                if [ "$first_in_batch" = false ]; then
                    echo "," >> "$current_batch_file"
                else
                    first_in_batch=false
                fi
                
                # Add file to current batch
                echo "$file_json" >> "$current_batch_file"
                current_batch_size=$((current_batch_size + file_size))
                
                echo "Added file $file_path (${file_size} bytes) to batch $batch_num"
            done < <(jq -c '.[]' sorted_non_paired.json)
        else
            echo "No non-paired files found to process"
        fi

        # Close the last batch if not empty
        if [ "$first_in_batch" = false ]; then
            echo "]" >> "$current_batch_file"
            batch_count=$((batch_count + 1))
            echo "Closed final batch $batch_num with size ${current_batch_size} bytes"
        else
            # Remove the empty file
            rm -f "$current_batch_file"
            echo "No files added to final batch, removing empty batch file"
        fi

        # List all batch files
        echo "Checking for created batch files..."
        batch_files=$(find batches -name "batch_*.json" 2>/dev/null | sort)
        if [ -n "$batch_files" ]; then
            echo "$batch_files" > batch_files.txt
            
            # Generate summary for each batch
            echo "Batch summary:"
            for batch_file in $(cat batch_files.txt); do
                if [ -s "$batch_file" ]; then
                    file_count=$(jq 'length' "$batch_file" 2>/dev/null || echo "0")
                    batch_size_bytes=$(jq 'map(.size_bytes) | add // 0' "$batch_file" 2>/dev/null || echo "0")
                    if [ "$batch_size_bytes" != "0" ] && [ -n "$batch_size_bytes" ]; then
                        batch_size_gb=$(echo "scale=2; $batch_size_bytes / 1073741824" | bc)
                    else
                        batch_size_gb="0.00"
                    fi
                    
                    # Check for paired files in this batch
                    paired_count=$(jq '[.[] | select(.relative_path | endswith("_R1_001.fastq.gz") or endswith("_R2_001.fastq.gz"))] | length' "$batch_file" 2>/dev/null || echo "0")
                    if [ "$paired_count" -gt 0 ] 2>/dev/null; then
                        pair_sets=$((paired_count / 2))
                        echo "Batch $(basename $batch_file): $file_count files ($pair_sets paired sets), $batch_size_gb GB"
                    else
                        echo "Batch $(basename $batch_file): $file_count files, $batch_size_gb GB"
                    fi
                else
                    echo "Batch $(basename $batch_file): empty or corrupted"
                fi
            done
        else
            echo "No batch files were created!"
            # Create an empty output to avoid workflow failure
            echo "[]" > batches/batch_0000.json
        fi

        echo "Created $batch_count batches"
        
        # Final summary
        echo ""
        echo "=== BATCH CREATION SUMMARY ==="
        echo "Total batches created: $batch_count"
        echo "Complete paired file sets processed: $paired_count"
        echo "Orphaned R2 files found: $orphaned_r2_count"
        if [ $orphaned_r2_count -gt 0 ]; then
            echo "WARNING: Orphaned R2 files have been excluded from batches and saved for manual review."
        fi
        echo "==============================="
    >>>

    output {
        Array[File] batch_files = glob("batches/batch_*.json")
        File orphaned_r2_files = "orphaned_r2_files.json"
        File orphaned_r2_summary = "orphaned_r2_summary.txt"
    }

    runtime {
        docker: docker_image
        cpu: 4
        memory: "8 GB"
        disks: "local-disk 20 SSD"
    }
}




task ProcessBatch {
    input {
        File batch_file
        String gs_bucket_path
        String project_id
        Int max_concurrent_uploads
        Int chunk_size_mb
        Int check_interval_seconds
        String? basespace_access_token
        String docker_image
        String config_file = ""
    }

    command <<<
        # Install necessary tools
        apt-get update && apt-get install -y wget curl jq unzip bc

        # Create temp directory for downloads
        mkdir -p downloaded logs

        # Install BaseSpace CLI
        wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs"
        chmod +x bs
        mv bs /usr/local/bin/

        # Set up authentication if token is provided
        if [ ! -z "~{basespace_access_token}" ]; then
            # Create BaseSpace config directory
            mkdir -p ~/.basespace
            
            # Copy config file if provided
            if [ ! -z "~{config_file}" ]; then
                gsutil cp ~{config_file} ~/.basespace/default.cfg
            fi
            
            # Authenticated with BaseSpace
            echo "Authenticated with BaseSpace..."
        else
            echo "Using existing BaseSpace authentication..."
        fi

        # Get batch name for logging
        batch_name=$(basename ~{batch_file})
        # Remove the .txt extension for dataset name
        dataset_name=${batch_name%.txt}
        echo "Processing batch: $batch_name"

        # Parse the batch file 
        file_count=$(jq length ~{batch_file})
        echo "Batch contains $file_count files"

        # Check if project_id is empty
        if [ -z "~{project_id}" ]; then
            echo "ERROR: No valid project ID available" >> upload_errors.log
            echo "FAILED" > result.txt
            echo "Failed to upload batch $batch_name - No valid project ID"
            
            # Create batch result JSON for failed upload
            timestamp=$(date '+%Y-%m-%d %H:%M:%S')
            cat <<EOF > batch_result.json
{
  "batch_name": "$batch_name",
  "total_files": $file_count,
  "result": "FAILED",
  "error": "No valid project ID available",
  "timestamp": "$timestamp"
}
EOF
            # Create empty log file
            echo "ERROR: No valid project ID available" > batch_upload.log
            exit 0
        fi

        # Extract all files in the batch
        jq -r '.[] | .full_path + "," + .relative_path' ~{batch_file} > batch_files.csv

        # Clean up any previous downloads
        rm -rf downloaded/*

        # Download all files in the batch and rename fastq/fastq.gz files
        while IFS=, read -r gs_path rel_path; do
            # Create the directory structure
            dir_path=$(dirname "$rel_path")
            mkdir -p "downloaded/$dir_path"
            
            echo "Downloading: $gs_path to downloaded/$rel_path"
            if ! gsutil cp "$gs_path" "downloaded/$rel_path"; then
                echo "ERROR: Failed to download $gs_path" >> download_errors.log
                continue
            fi
            
            # Check if file is fastq or fastq.gz and rename if needed
            if [[ "$rel_path" == *.fastq ]] || [[ "$rel_path" == *.fastq.gz ]]; then
                echo "Found fastq file: $rel_path - checking for renaming"
                
                # Get file name without directory
                file_name=$(basename "$rel_path")
                
                # First check: ensure file has at least 7 underscores total
                total_underscores=$(echo "$file_name" | tr -cd '_' | wc -c)
                
                if [ $total_underscores -ge 7 ]; then  # Only proceed if there are at least 7 underscores
                    echo "File has $total_underscores underscores, proceeding with renaming"
                    
                    # Strip everything before and including the 3rd underscore
                    # First, check if we have at least 3 underscores
                    if [ $total_underscores -ge 3 ]; then
                        # Extract part after third underscore
                        remainder=$(echo "$file_name" | cut -d'_' -f4-)
                        
                        # Define file extension pattern to support both .fastq and .fastq.gz
                        ext_pattern="\.fastq(\.gz)?"
                        
                        # Updated pattern for Illumina fastq files:
                        # - Include _S[0-9]+ as part of the suffix pattern to preserve
                        # - Support both .fastq and .fastq.gz extensions
                        # Pattern looks for: prefix_S[digits]_L[digits]_R[digits]_[digits].fastq[.gz]
                        if [[ "$remainder" =~ (.*)(_S[0-9]+_L[0-9]+_R[0-9]+_[0-9]+${ext_pattern}) ]]; then
                            prefix="${BASH_REMATCH[1]}"
                            suffix="${BASH_REMATCH[2]}"
                            
                            # Replace all underscores with hyphens in the prefix part
                            modified_prefix=${prefix//_/-}
                            
                            # Combine back with the suffix that keeps its underscores
                            new_name="${modified_prefix}${suffix}"
                            
                            # Debug output
                            echo "Original: $remainder"
                            echo "Prefix: $prefix (will replace _ with -)"
                            echo "Suffix: $suffix (will preserve _)"
                            echo "New name: $new_name"
                            
                            # Move file to new name
                            mv "downloaded/$rel_path" "downloaded/$dir_path/$new_name"
                        else
                            echo "Could not identify standard Illumina fastq pattern in $remainder, keeping original name after stripping first 3 fields"
                            mv "downloaded/$rel_path" "downloaded/$dir_path/$remainder"
                        fi
                    else
                        echo "File $file_name doesn't have at least 3 underscores, keeping original name"
                    fi
                else
                    echo "File $file_name has only $total_underscores underscores (needs at least 7), skipping renaming"
                fi
            fi
        done < batch_files.csv

        # Upload the entire batch to BaseSpace
        echo "Uploading batch to BaseSpace project ID: ~{project_id} with dataset name: $dataset_name"


        # Upload FASTQ files first (if any exist)
        FASTQ_FILES=$(find downloaded -type f -name "*.fastq*" | wc -l)
        if [ $FASTQ_FILES -gt 0 ]; then
            echo "Uploading $FASTQ_FILES FASTQ files to BaseSpace"
            # Create a list of FASTQ files to upload
            find downloaded -type f -name "*.fastq*" > fastq_list.txt
            
            # Use the bs CLI to upload all FASTQ files in one command
            if bs upload dataset \
                --project ~{project_id} \
                --name "${dataset_name}_fastq" \
                --recursive \
                downloaded/ > logs/fastq_upload.log 2>&1; then
                
                # Mark each file as successful
                while read -r file_path; do
                    rel_path=${file_path#downloaded/}
                    gs_path=$(grep ",$rel_path$" batch_files.csv | cut -d, -f1)
                    echo "$gs_path,SUCCESS" >> results.csv
                done < fastq_list.txt
                echo "SUCCESS" > result.txt
                echo "Successfully uploaded common files in batch $batch_name as dataset ${dataset_name}_common_files"
            else
                # Mark each file as failed
                while read -r file_path; do
                    rel_path=${file_path#downloaded/}
                    gs_path=$(grep ",$rel_path$" batch_files.csv | cut -d, -f1)
                    echo "$gs_path,UPLOAD_FAILED" >> results.csv
                done < fastq_list.txt
                echo "FAILED" > result.txt
                echo "Failed to upload batch ${batch_name}_fastq"
                cat logs/fastq_upload.log >> upload_errors.log
            fi
        else
            echo "No FASTQ files found to upload"
        fi
        
        # Upload all non-FASTQ files (if any exist)
        NON_FASTQ_FILES=$(find downloaded -type f ! -name "*.fastq*" | wc -l)
        if [ $NON_FASTQ_FILES -gt 0 ]; then
            echo "Uploading $NON_FASTQ_FILES non-FASTQ files to BaseSpace"
            # Create a list of non-FASTQ files to upload
            find downloaded -type f ! -name "*.fastq*" > common_list.txt
            
            # Use the bs CLI to upload all non-FASTQ files in one command
            if bs upload dataset \
                --project ~{project_id} \
                --name "${dataset_name}_common_files" \
                --recursive \
                --type common.files \
                --exclude "*.fastq*" \
                downloaded/ > logs/common_upload.log 2>&1; then
                
                # Mark each file as successful
                while read -r file_path; do
                    rel_path=${file_path#downloaded/}
                    gs_path=$(grep ",$rel_path$" batch_files.csv | cut -d, -f1)
                    echo "$gs_path,SUCCESS" >> results.csv
                done < common_list.txt
                echo "SUCCESS" > result.txt
                echo "Successfully uploaded common files in batch $batch_name as dataset ${dataset_name}_common_files"
            else
                # Mark each file as failed
                while read -r file_path; do
                    rel_path=${file_path#downloaded/}
                    gs_path=$(grep ",$rel_path$" batch_files.csv | cut -d, -f1)
                    echo "$gs_path,UPLOAD_FAILED" >> results.csv
                done < common_list.txt
                echo "FAILED" > result.txt
                echo "Failed to upload batch ${batch_name}_fastq"
                cat logs/common_upload.log >> upload_errors.log
            fi
        else
            echo "No non-FASTQ files found to upload"
        fi
        
        # Combine logs
        cat logs/*.log > batch_upload.log
        
        # Create result summary
        timestamp=$(date '+%Y-%m-%d %H:%M:%S')
        result=$(cat result.txt)
        
        # Create batch result JSON
        cat <<EOF > batch_result.json
{
  "batch_name": "$batch_name",
  "dataset_name": "$dataset_name",
  "total_files": $file_count,
  "result": "$result",
  "timestamp": "$timestamp"
}
EOF
    >>>

    output {
        File batch_result = "batch_result.json"
        File upload_log = "batch_upload.log"
        File results_detail = "results.csv"
    }

    runtime {
        docker: docker_image
        cpu: 8
        memory: "16 GB"
        disks: "local-disk 100 SSD"
        preemptible: 1
    }
}

task VerifyAllUploads {
    input {
        Array[File] batch_results
        String project_id
        String? basespace_access_token
        String docker_image
        String config_file = ""
    }

    command <<<
        # Install necessary tools
        apt-get update && apt-get install -y wget curl jq unzip bc

        # Install BaseSpace CLI
        wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs"
        chmod +x bs
        mv bs /usr/local/bin/

        # Set up authentication if token is provided
        if [ ! -z "~{basespace_access_token}" ]; then
            # Create BaseSpace config directory
            mkdir -p ~/.basespace
            
            # Copy config file if provided
            if [ ! -z "~{config_file}" ]; then
                gsutil cp ~{config_file} ~/.basespace/default.cfg
            fi
            
            # Authenticate with BaseSpace
            echo "Authenticating with BaseSpace..."
        else
            echo "Using existing BaseSpace authentication..."
        fi

        # Check if project_id is empty
        if [ -z "~{project_id}" ]; then
            echo "ERROR: No valid project ID available for verification" >> verification_errors.log
            
            # Create empty datasets file
            echo "[]" > datasets.json
            
            # Create result summary with error
            cat <<EOF > summary.txt
===== Upload Summary =====
ERROR: No valid project ID available for verification
Total batches processed: 0
Successfully uploaded batches: 0
Failed batch uploads: 0
Total files processed: 0
Datasets in BaseSpace project: 0
========================
EOF

            # Create detailed results JSON
            cat <<EOF > detailed_results.json
{
  "batches": [],
  "summary": {
    "error": "No valid project ID available for verification",
    "total_batches": 0,
    "successful_batches": 0,
    "failed_batches": 0,
    "total_files": 0,
    "datasets_in_basespace": 0
  }
}
EOF
            cat summary.txt
            exit 0
        fi

        # Check project status
        echo "Verifying uploads in BaseSpace project ID: ~{project_id}"
        bs list datasets --project-id=~{project_id} --format=json > datasets.json

        # Count datasets in BaseSpace
        if [ -s datasets.json ]; then
            dataset_count=$(jq 'length' datasets.json)
        else
            dataset_count=0
        fi

        # Combine and analyze batch results
        total_batches=0
        successful_batches=0
        failed_batches=0
        total_files=0

        # Create detailed results JSON
        echo "{" > detailed_results.json
        echo "  \"batches\": [" >> detailed_results.json

        # Process each batch result file
        first=true
        for batch_file in ~{sep=" " batch_results}; do
            total_batches=$((total_batches + 1))
            
            # Extract result from the batch
            result=$(jq -r '.result' "$batch_file")
            batch_files=$(jq -r '.total_files' "$batch_file")
            total_files=$((total_files + batch_files))
            
            if [[ "$result" == "SUCCESS" ]]; then
                successful_batches=$((successful_batches + 1))
            else
                failed_batches=$((failed_batches + 1))
            fi
            
            # Add batch details to JSON
            if [ "$first" = true ]; then
                first=false
            else
                echo "," >> detailed_results.json
            fi
            
            cat "$batch_file" >> detailed_results.json
        done

        # Close the batches array
        echo "  ]," >> detailed_results.json
        
        # Add summary section
        cat <<EOF >> detailed_results.json
  "summary": {
    "total_batches": $total_batches,
    "successful_batches": $successful_batches,
    "failed_batches": $failed_batches,
    "total_files": $total_files,
    "datasets_in_basespace": $dataset_count
  }
}
EOF

        # Create human-readable summary
        cat <<EOF > summary.txt
===== Upload Summary =====
Total batches processed: $total_batches
Successfully uploaded batches: $successful_batches
Failed batch uploads: $failed_batches
Total files processed: $total_files
EOF

        # Calculate success rate
        if [ $total_batches -gt 0 ]; then
            success_rate=$(echo "scale=2; ($successful_batches * 100) / $total_batches" | bc)
            echo "Success rate: ${success_rate}%" >> summary.txt
        fi
        echo "Datasets in BaseSpace project: $dataset_count" >> summary.txt
        echo "========================" >> summary.txt

        cat summary.txt
    >>>

    output {
        String summary = read_string("summary.txt")
        File detailed_results = "detailed_results.json"
        File basespace_datasets = "datasets.json"
    }

    runtime {
        docker: docker_image
        cpu: 2
        memory: "4 GB"
        disks: "local-disk 10 SSD"
    }
}