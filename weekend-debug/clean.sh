#!/bin/bash

# Define the directory to clean (default is the current directory)
CLEAN_DIR=${1:-.}

# Define the file extensions to remove
EXTENSIONS=(".xmf" ".h5" ".res")

# Loop over the extensions and remove corresponding files
for ext in "${EXTENSIONS[@]}"; do
    find "$CLEAN_DIR" -type f -name "*$ext" -exec rm -f {} \;
done

echo "Cleanup complete."