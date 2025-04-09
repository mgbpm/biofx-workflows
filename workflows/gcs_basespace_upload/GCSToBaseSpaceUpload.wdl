version 1.0

workflow GCSToBaseSpaceUpload {
    input {
        String gs_bucket_path
        String basespace_project_name
        Float max_batch_size_gb = 50.0
        Int max_concurrent_uploads_per_vm = 10
        Int chunk_size_mb = 25
        Int check_interval_seconds = 60
        String basespace_access_token
        String docker_image = "google/cloud-sdk:latest"
    }

   #First create basespace config file
    call create_config_file {
        input:
            basespace_access_token = basespace_access_token
    }

    # First create the BaseSpace project
    call CreateBaseSpaceProject {
        input:
            project_name = basespace_project_name,
            basespace_access_token = basespace_access_token,
            docker_image = docker_image,
            config_file = create_config_file.config_file
    }

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
                project_id = CreateBaseSpaceProject.project_id,
                max_concurrent_uploads = max_concurrent_uploads_per_vm,
                chunk_size_mb = chunk_size_mb,
                check_interval_seconds = check_interval_seconds,
                basespace_access_token = basespace_access_token,
                docker_image = docker_image,
                config_file = create_config_file.config_file
        }
    }

    # Combine all results to verify completion
    call VerifyAllUploads {
        input:
            batch_results = ProcessBatch.batch_result,
            project_id = CreateBaseSpaceProject.project_id,
            basespace_access_token = basespace_access_token,
            docker_image = docker_image,
            config_file = create_config_file.config_file
    }

    output {
        String basespace_project_id = CreateBaseSpaceProject.project_id
        String basespace_project_url = CreateBaseSpaceProject.project_url
        Array[File] batch_upload_logs = ProcessBatch.upload_log
        String upload_summary = VerifyAllUploads.summary
        File detailed_results = VerifyAllUploads.detailed_results
    }
}

task CreateConfigFile {
    input {
        String basespace_access_token
    }

    command <<<
        echo "apiServer   = https://api.basespace.illumina.com" > default.cfg
        echo "accessToken = ~{basespace_access_token}" >> default.cfg
    >>>

    output {
        File config_file = "default.cfg"
    }

    runtime {
        docker: "ubuntu:latest"
    }
}

task CreateBaseSpaceProject {
    input {
        String project_name
        String? basespace_access_token
        String docker_image
        File config_file
    }

    command <<<
        # Install BaseSpace CLI
        apt-get update && apt-get install -y wget curl jq unzip

        # Install BaseSpace CLI
        wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs"
        #get "https://launch.basespace.illumina.com/CLI/latest/amd64-osx/bs"
        mv bs /usr/local/bin/
        chmod +x /usr/local/bin/bs

        # Set up authentication if token is provided
        if [ ! -z "~{basespace_access_token}" ]; then
            # Copy config to expected working dir
            cp ~{config_file} ~/.basespace/default.cfg
        else
            echo "No access token provided, assuming BaseSpace CLI is already configured"
        fi

        # Create project
        echo "Creating BaseSpace project: ~{project_name}"
        RESPONSE=$(bs create project --name="~{project_name}" --format=json)
        echo $RESPONSE > project_response.json

        # Extract project ID and URL
        PROJECT_ID=$(echo $RESPONSE | jq -r '.Response.Id')
        PROJECT_URL=$(echo $RESPONSE | jq -r '.Response.HrefBaseSpaceUI')

        echo $PROJECT_ID > project_id.txt
        echo $PROJECT_URL > project_url.txt
        echo "Created project with ID: $PROJECT_ID"
    >>>

    output {
        String project_id = read_string("project_id.txt")
        String project_url = read_string("project_url.txt")
        File project_details = "project_response.json"
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
        apt-get update && apt-get install -y jq

        # Ensure the bucket path ends with a slash for directory listing
        if [[ "~{gs_bucket_path}" != */ ]]; then
            GS_PATH="~{gs_bucket_path}/"
        else
            GS_PATH="~{gs_bucket_path}"
        fi

        # Get stats for all files in the bucket path
        echo "Listing files in: $GS_PATH"
        gsutil -m ls -l "${GS_PATH}**" | grep -v '/$' > all_files_raw.txt

        # Process the files with sizes and create JSON manifest
        echo "[" > file_manifest.json
        first_line=true

        while IFS= read -r line; do
            # Skip directory markers and summary lines
            if [[ "$line" == *"TOTAL:"* ]] || [[ "$line" == *"gs://"/$ ]]; then
                continue
            fi

            # Extract size and path
            size_bytes=$(echo "$line" | awk '{print $1}')
            full_path=$(echo "$line" | awk '{print $NF}')
            
            # Skip if not a valid size or path
            if [[ -z "$size_bytes" || -z "$full_path" || ! "$size_bytes" =~ ^[0-9]+$ ]]; then
                continue
            fi

            # Calculate relative path by removing bucket path prefix
            relative_path=${full_path#$GS_PATH}
            
            # Calculate size in GB
            size_gb=$(echo "scale=9; $size_bytes / 1073741824" | bc)

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
        done < all_files_raw.txt

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
        max_batch_size_bytes=$(echo "~{max_batch_size_gb} * 1073741824" | bc)

        # Create a working directory for batch files
        mkdir -p batches

        # First sort files by size (largest first) for better packing
        jq -r 'sort_by(.size_bytes) | reverse | .[]' ~{file_manifest} > sorted_files.json

        batch_num=0
        batch_files=()
        current_batch_file="batches/batch_$(printf "%04d" $batch_num).json"
        current_batch_size=0
        batch_count=0

        # Start the first batch
        echo "[" > "$current_batch_file"
        first_in_batch=true

        # Process each file
        while IFS= read -r file_json; do
            # Extract file size
            file_size=$(echo "$file_json" | jq -r '.size_bytes')

            # If this single file is larger than max batch size, it gets its own batch
            if (( file_size > max_batch_size_bytes )); then
                # If we've already started a batch, close it
                if [ "$first_in_batch" = false ]; then
                    echo "]" >> "$current_batch_file"
                    batch_count=$((batch_count + 1))
                else
                    # Remove the empty batch file we created
                    rm "$current_batch_file"
                fi

                # Create a new batch for this single file
                batch_num=$((batch_num + 1))
                current_batch_file="batches/batch_$(printf "%04d" $batch_num).json"
                echo "[" > "$current_batch_file"
                echo "$file_json" >> "$current_batch_file"
                echo "]" >> "$current_batch_file"
                batch_count=$((batch_count + 1))

                # Start a new batch
                batch_num=$((batch_num + 1))
                current_batch_file="batches/batch_$(printf "%04d" $batch_num).json"
                echo "[" > "$current_batch_file"
                first_in_batch=true
                current_batch_size=0
                continue
            fi

            # If adding this file would exceed the batch size, start a new batch
            if (( current_batch_size + file_size > max_batch_size_bytes )) && [ "$first_in_batch" = false ]; then
                echo "]" >> "$current_batch_file"
                batch_count=$((batch_count + 1))

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
        done < <(jq -c '.[]' sorted_files.json)

        # Close the last batch if not empty
        if [ "$first_in_batch" = false ]; then
            echo "]" >> "$current_batch_file"
            batch_count=$((batch_count + 1))
        else
            # Remove the empty file
            rm "$current_batch_file"
        fi

        # List all batch files
        find batches -name "batch_*.json" | sort > batch_files.txt

        # Generate summary for each batch
        for batch_file in $(cat batch_files.txt); do
            file_count=$(jq 'length' "$batch_file")
            batch_size_bytes=$(jq 'map(.size_bytes) | add' "$batch_file")
            batch_size_gb=$(echo "scale=2; $batch_size_bytes / 1073741824" | bc)
            echo "Batch $(basename $batch_file): $file_count files, $batch_size_gb GB"
        done

        echo "Created $batch_count batches"
    >>>

    output {
        Array[File] batch_files = glob("batches/batch_*.json")
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
        File config_file
    }

    command <<<
        # Install necessary tools
        apt-get update && apt-get install -y wget curl jq unzip

        # Install BaseSpace CLI
        wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs"
        #wget "https://launch.basespace.illumina.com/CLI/latest/amd64-osx/bs"
        mv bs /usr/local/bin/
        chmod +x /usr/local/bin/bs

        # Set up authentication if token is provided
        if [ ! -z "~{basespace_access_token}" ]; then
            # Copy config to expected working dir
            cp ~{config_file} ~/.basespace/default.cfg
        else
            echo "No access token provided, assuming BaseSpace CLI is already configured"
        fi

        # Create directories
        mkdir -p downloaded logs

        # Extract file paths from batch file
        jq -r '.[] | .full_path + "," + .relative_path' ~{batch_file} > batch_files.csv

        # Initialize results file
        echo "file_path,result" > results.csv

        # Helper function to upload a file
        upload_file() {
            local gs_path=$1
            local rel_path=$2

            # Extract filename and create directory structure
            local filename=$(basename "$gs_path")
            local dir_path=$(dirname "$rel_path")

            # Create local directory for the file
            mkdir -p "downloaded/$dir_path"
            local local_path="downloaded/$rel_path"

            echo "Downloading: $gs_path to $local_path"
            if ! gsutil cp "$gs_path" "$local_path"; then
                echo "$gs_path,DOWNLOAD_FAILED" >> results.csv
                return 1
            fi

            echo "Uploading: $local_path to BaseSpace project"
            # Use the bs CLI to upload the file, preserving directory structure
            if bs upload dataset \
                --project-id ~{project_id} \
                --chunk-size ~{chunk_size_mb} \
                --remote-path "$dir_path" \
                --no-prompt \
                --verbose \
                "$local_path" > "logs/${filename}.log" 2>&1; then

                echo "$gs_path,SUCCESS" >> results.csv
                # Clean up downloaded file to save disk space
                rm "$local_path"
            else
                echo "$gs_path,UPLOAD_FAILED" >> results.csv
                # Keep the log for debugging
                cat "logs/${filename}.log" >> upload_errors.log
            fi
        }

        export -f upload_file

        # Process all files in parallel with xargs 
        # using the number of CPUs for parallelization
        cat batch_files.csv | \
        xargs -P ~{max_concurrent_uploads} -I{} bash -c 'upload_file $(echo {} | cut -d, -f1) $(echo {} | cut -d, -f2)'

        # Create summary of batch results
        echo "{" > batch_result.json
        echo "  \"total_files\": $(cat batch_files.csv | wc -l)," >> batch_result.json
        echo "  \"processed_files\": $(cat results.csv | grep -v "file_path" | wc -l)," >> batch_result.json
        echo "  \"successful_uploads\": $(cat results.csv | grep SUCCESS | wc -l)," >> batch_result.json
        echo "  \"failed_uploads\": $(cat results.csv | grep -v SUCCESS | grep -v "file_path" | wc -l)" >> batch_result.json
        echo "}" >> batch_result.json

        # Combine logs
        cat logs/*.log > batch_upload.log
    >>>

    output {
        File batch_result = "batch_result.json"
        File results_detail = "results.csv"
        File upload_log = "batch_upload.log"
    }

    runtime {
        docker: docker_image
        cpu: 8
        memory: "16 GB"
        disks: "local-disk 500 SSD"
        preemptible: 1
    }
}

task VerifyAllUploads {
    input {
        Array[File] batch_results
        String project_id
        String? basespace_access_token
        String docker_image
        File config_file
    }

    command <<<
        # Install necessary tools
        apt-get update && apt-get install -y wget curl jq unzip

        # Install BaseSpace CLI
        wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs"
        #wget "https://launch.basespace.illumina.com/CLI/latest/amd64-osx/bs"
        mv bs /usr/local/bin/
        chmod +x /usr/local/bin/bs

        # Set up authentication if token is provided
        if [ ! -z "~{basespace_access_token}" ]; then
            # Copy config to expected working dir
            cp ~{config_file} ~/.basespace/default.cfg
        else
            echo "No access token provided, assuming BaseSpace CLI is already configured"
        fi

        # Check project status
        bs list datasets --project-id=~{project_id} --format=json > datasets.json

        # Combine and analyze batch results
        total_files=0
        processed_files=0
        successful_uploads=0
        failed_uploads=0

        # Process each batch result file
        for batch_file in ~{sep=" " batch_results}; do
            # Extract metrics from the batch result JSON
            batch_total=$(jq -r '.total_files' "$batch_file")
            batch_processed=$(jq -r '.processed_files' "$batch_file")
            batch_success=$(jq -r '.successful_uploads' "$batch_file")
            batch_failed=$(jq -r '.failed_uploads' "$batch_file")

            # Add to totals
            total_files=$((total_files + batch_total))
            processed_files=$((processed_files + batch_processed))
            successful_uploads=$((successful_uploads + batch_success))
            failed_uploads=$((failed_uploads + batch_failed))
        done

        # Count datasets in BaseSpace
        if [ -s datasets.json ]; then
            dataset_count=$(jq 'length' datasets.json)
        else
            dataset_count=0
        fi

        # Create detailed results JSON
        cat <<EOF > detailed_results.json
{
  "summary": {
    "total_files": $total_files,
    "processed_files": $processed_files,
    "successful_uploads": $successful_uploads,
    "failed_uploads": $failed_uploads,
    "datasets_in_basespace": $dataset_count
  },
  "batch_results": [
EOF

        # Add each batch result
        first=true
        for batch_file in ~{sep=" " batch_results}; do
            if [ "$first" = true ]; then
                first=false
            else
                echo "," >> detailed_results.json
            fi
            cat "$batch_file" >> detailed_results.json
        done

        # Close the JSON
        echo "]}" >> detailed_results.json

        # Create human-readable summary
        cat <<EOF > summary.txt
Upload Summary:
----------------
Total files: $total_files
Files processed: $processed_files
Successfully uploaded: $successful_uploads
Failed uploads: $failed_uploads
Datasets in BaseSpace project: $dataset_count
EOF

        # Calculate success rate
        if [ $processed_files -gt 0 ]; then
            success_rate=$(echo "scale=2; ($successful_uploads * 100) / $processed_files" | bc)
            echo "Success rate: ${success_rate}%" >> summary.txt
        fi
        echo "----------------" >> summary.txt

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