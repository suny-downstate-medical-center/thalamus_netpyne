#!/bin/bash

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <path_to_folder> <subfolder_name> <string_to_match>"
    exit 1
fi

# Assign arguments to variables
FOLDER_PATH="$1"
SUBFOLDER_NAME="$2"
MATCH_STRING="$3"

# Create the subfolder if it doesn't exist
SUBFOLDER_PATH="$FOLDER_PATH/$SUBFOLDER_NAME"
mkdir -p "$SUBFOLDER_PATH"

# Copy matched files into the subfolder
for file in "$FOLDER_PATH"/*"$MATCH_STRING"*; do
    if [ -f "$file" ]; then  # Ensure it's a file, not a directory
        cp "$file" "$SUBFOLDER_PATH"
    fi
done

echo "Files containing '$MATCH_STRING' have been copied to '$SUBFOLDER_PATH'."
