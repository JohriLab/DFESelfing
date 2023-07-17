#!/bin/bash

output_file_fixed="/combine_output.fixed"
output_file_txt="/combine_output.txt"

# Function to combine lines without replacing "m1"
combine_lines_without_replace() {
    for file in "$1"/*.fixed; do
        grep "m1" "$file" >> "$directory""$output_file_fixed"
    done
    for file in "$1"/*.txt; do
        grep "m1" "$file" >> "$directory""$output_file_txt"
    done
}

# Function to combine lines with "m1" replaced by "m2"
combine_lines_with_replace() {
    for file in "$1"/*.fixed; do
        grep "m1" "$file" | sed 's/m1/m2/g' >> "$directory""$output_file_fixed"
    done
    for file in "$1"/*.txt; do
        grep "m1" "$file" | sed 's/m1/m2/g' >> "$directory""$output_file_txt"
    done
}

# Get target directory from user input
read -p "Enter the target directory: " directory

# Check if the directory exists
if [ ! -d "$directory" ]; then
    echo "Directory not found!"
    exit 1
fi

# Prompt for user option
echo "Choose an option:"
echo "1. Combine lines without replacing 'm1'"
echo "2. Combine lines with 'm1' replaced by 'm2'"
read -p "Enter your choice (1 or 2): " option

# Remove output files if they already exist
if [ -f "$directory""$output_file_fixed" ]; then
    rm "$directory""$output_file_fixed"
fi
if [ -f "$directory""$output_file_txt" ]; then
    rm "$directory""$output_file_txt"
fi
# Perform the selected option
case $option in
    1) combine_lines_without_replace "$directory" ;;
    2) combine_lines_with_replace "$directory" ;;
    *) echo "Invalid option!" ;;
esac

echo "Lines combined successfully. Output file: $directory$output_file_fixed"
echo "Here is what the first few lines of output look like (fixed mutations):"
head $directory$output_file_fixed

echo "Here is what the first few lines of output look like (segregating mutations):"
head $directory$output_file_txt