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
            docker_image = docker_image
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
            echo "Creating BaseSpace project: $project_name"
            RESPONSE=$(bs create project --name="$project_name" --format=json)
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
    }

    command <<<
        # Install required tools
        apt-get update && apt-get install -y jq bc

        # Convert max_batch_size_gb to bytes for comparison
        max_batch_size_bytes=$(echo "~{max_batch_size_gb} * 1073741824" | bc -l)
        echo "Max batch size: ~{max_batch_size_gb} GB ($max_batch_size_bytes bytes)"

        # Create a working directory for batch files
        mkdir -p batches

        # Debug: Print the file manifest structure
        echo "File manifest structure:"
        jq '.' ~{file_manifest} | head -n 20
        echo "..."

        # Count total files
        total_files=$(jq 'length' ~{file_manifest})
        echo "Total files in manifest: $total_files"

        # Sort files by size (largest first) to optimize packing
        jq -r 'sort_by(-.size_bytes) | .[] | [.full_path, .relative_path, .size_bytes, .size_gb] | @csv' ~{file_manifest} > sorted_files.csv
        echo "Sorted files by size (largest first)"

        # Initialize batching variables
        batch_num=0
        current_batch_file="batches/batch_$(printf "%04d" $batch_num).txt"
        current_batch_size=0
        file_count_in_current_batch=0
        total_batches=0

        # Create the first batch file
        echo "[]" > "$current_batch_file"

        # Process each file from the sorted list
        while IFS=, read -r full_path relative_path size_bytes size_gb; do
            # Remove quotes if present
            full_path=$(echo "$full_path" | sed 's/"//g')
            relative_path=$(echo "$relative_path" | sed 's/"//g')
            size_bytes=$(echo "$size_bytes" | sed 's/"//g')
            
            echo "Processing file: $relative_path ($size_bytes bytes)"
            
            # If this single file is larger than max batch size, it gets its own batch
            # Use bc for comparison, not bash arithmetic
            is_larger=$(echo "$size_bytes > $max_batch_size_bytes" | bc -l)
            if [ "$is_larger" -eq 1 ]; then
                echo "Large file detected: $relative_path exceeds batch size limit"
                
                # If we've already added files to the current batch, close it
                if [ "$file_count_in_current_batch" -gt 0 ]; then
                    echo "Closing current batch with $file_count_in_current_batch files"
                    total_batches=$((total_batches + 1))
                    batch_num=$((batch_num + 1))
                    current_batch_file="batches/batch_$(printf "%04d" $batch_num).txt"
                    echo "[]" > "$current_batch_file"
                    current_batch_size=0
                    file_count_in_current_batch=0
                fi

                # Create a JSON entry for this large file batch
                cat <<EOF > "$current_batch_file"
[
  {
    "full_path": "$full_path",
    "relative_path": "$relative_path",
    "size_bytes": $size_bytes
  }
]
EOF
                
                # Move to next batch
                total_batches=$((total_batches + 1))
                batch_num=$((batch_num + 1))
                current_batch_file="batches/batch_$(printf "%04d" $batch_num).txt"
                echo "[]" > "$current_batch_file"
                current_batch_size=0
                file_count_in_current_batch=0
                continue
            fi

            # Calculate new batch size
            new_batch_size=$(echo "$current_batch_size + $size_bytes" | bc -l)
            
            # Check if adding this file would exceed the batch size
            # Use bc for comparison, not bash arithmetic
            would_exceed=$(echo "$new_batch_size > $max_batch_size_bytes" | bc -l)
            if [ "$would_exceed" -eq 1 ] && [ "$file_count_in_current_batch" -gt 0 ]; then
                echo "Batch size limit reached. Creating new batch."
                total_batches=$((total_batches + 1))
                batch_num=$((batch_num + 1))
                current_batch_file="batches/batch_$(printf "%04d" $batch_num).txt"
                echo "[]" > "$current_batch_file"
                current_batch_size=0
                file_count_in_current_batch=0
            fi

            # Add file to current batch
            if [ "$file_count_in_current_batch" -eq 0 ]; then
                # First file in the batch
                cat <<EOF > "$current_batch_file"
[
  {
    "full_path": "$full_path",
    "relative_path": "$relative_path",
    "size_bytes": $size_bytes
  }
]
EOF
            else
                # Add file to existing batch
                # Remove the last line (closing bracket)
                sed -i '$d' "$current_batch_file"
                
                # Add comma and new file entry
                cat <<EOF >> "$current_batch_file"
,
  {
    "full_path": "$full_path",
    "relative_path": "$relative_path",
    "size_bytes": $size_bytes
  }
]
EOF
            fi
            
            # Update current batch size using bc for floating point arithmetic
            current_batch_size=$(echo "$current_batch_size + $size_bytes" | bc -l)
            file_count_in_current_batch=$((file_count_in_current_batch + 1))
            echo "Added file to batch $batch_num. Current batch size: $current_batch_size bytes with $file_count_in_current_batch files"
        done < sorted_files.csv

        # Update the final batch count if the last batch has files
        if [ "$file_count_in_current_batch" -gt 0 ]; then
            total_batches=$((total_batches + 1))
        else
            # Remove empty batch file
            rm -f "$current_batch_file"
        fi

        # List all batch files and count the number of files in each
        echo "Final batch summary:"
        for batch_file in batches/batch_*.txt; do
            if [ -f "$batch_file" ]; then
                file_count=$(jq '. | length' "$batch_file")
                batch_size=$(jq 'reduce .[] as $item (0; . + ($item.size_bytes // 0))' "$batch_file" 2>/dev/null || echo "Unknown")
                batch_size_gb=$(echo "scale=4; $batch_size / 1073741824" | bc -l 2>/dev/null || echo "Unknown")
                echo "Batch $(basename $batch_file): $file_count files, approx. $batch_size bytes ($batch_size_gb GB)"
            fi
        done

        echo "Created $total_batches batches in total"
        find batches -name "batch_*.txt" | sort > batch_files.txt
        echo "Found $(wc -l < batch_files.txt) batch files"
    >>>

    output {
        Array[File] batch_files = glob("batches/batch_*.txt")
    }

    runtime {
        docker: docker_image
        cpu: 2
        memory: "4 GB"
        disks: "local-disk 10 SSD"
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

        # Download all files in the batch
        while IFS=, read -r gs_path rel_path; do
            # Create the directory structure
            dir_path=$(dirname "$rel_path")
            mkdir -p "downloaded/$dir_path"
            
            echo "Downloading: $gs_path to downloaded/$rel_path"
            if ! gsutil cp "$gs_path" "downloaded/$rel_path"; then
                echo "ERROR: Failed to download $gs_path" >> download_errors.log
            fi
        done < batch_files.csv

        # Upload the entire batch to BaseSpace
        echo "Uploading batch to BaseSpace project ID: ~{project_id}"
        
        if bs upload dataset \
            --project ~{project_id} \
            --recursive \
            --type common.files \
            downloaded/ > "batch_upload.log" 2>&1; then
            
            echo "SUCCESS" > result.txt
            echo "Successfully uploaded batch $batch_name"
        else
            echo "FAILED" > result.txt
            echo "Failed to upload batch $batch_name"
            cat batch_upload.log >> upload_errors.log
        fi
        
        # Create result summary
        timestamp=$(date '+%Y-%m-%d %H:%M:%S')
        result=$(cat result.txt)
        
        # Create batch result JSON
        cat <<EOF > batch_result.json
{
  "batch_name": "$batch_name",
  "total_files": $file_count,
  "result": "$result",
  "timestamp": "$timestamp"
}
EOF
    >>>

    output {
        File batch_result = "batch_result.json"
        File upload_log = "batch_upload.log"
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