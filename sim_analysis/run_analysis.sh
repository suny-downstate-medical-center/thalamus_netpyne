#!/bin/bash

# Loop through all .py files in the current directory
for file in *.py; do
    if [ -f "$file" ]; then
        echo "Running $file..."
        python "$file"
    fi
done