#!/bin/bash

# Output file
output_file="job_details.txt"

# Initialize the output file with a header
echo "job-ID,submit_cmd,exec_host" > $output_file

# List all job IDs
job_ids=$(qstat -u joao | awk 'NR>2 {print $1}') # Adjust column if needed

# Loop through each job ID
for job_id in $job_ids; do
    # Extract the submit_cmd for each job
    submit_cmd=$(qstat -j $job_id | grep "submit_cmd" | awk -F ": " '{print $2}')
    
    # Extract the exec_host for each job, if available
    exec_host=$(qstat -j $job_id | grep "exec_host" | awk -F ": " '{print $2}')
    
    # Handle cases where exec_host is unavailable
    if [ -z "$exec_host" ]; then
        exec_host="N/A"
    fi
    
    # Append job-ID, submit_cmd, and exec_host to the output file
    echo "$job_id,$submit_cmd,$exec_host" >> $output_file
done

echo "Job details saved to $output_file"
