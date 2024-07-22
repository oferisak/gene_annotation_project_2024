#!/bin/bash

input_file=$1
output_folder=${2:-.}  # Default to current folder if no output folder is provided

# Extract the base name of the input file (strip path)
input_file_basename=$(basename "$input_file")

# Construct the output file path
output_file="${output_folder}/hg19_${input_file_basename%.gz}"  # Remove .gz extension if present

# Determine if the input file is gzipped
if [[ "$input_file" == *.gz ]]; then
    read_cmd="zcat"
else
    read_cmd="cat"
fi

# Process the file
$read_cmd "$input_file" | while IFS= read -r line || [[ -n "$line" ]]; do
    # Check if the line is a header or the column definition line
    if [[ "$line" == \##* ]] || [[ "$line" == \#CHROM* ]]; then
        # Write the header or column definition line as is to the output file
        echo "$line" >> "$output_file"
    elif [[ "$line" =~ ^MT ]]; then
        # Replace 'MT' with 'chrM' and write to the output file
        echo "chrM${line:2}" >> "$output_file"
    elif [[ "$line" =~ ^[0-9XY] || "$line" =~ ^chr[0-9XY] ]]; then
        # Prepend 'chr' to numeric chromosome lines if not already present and write to the output file
        if [[ "$line" =~ ^chr ]]; then
            echo "$line" >> "$output_file"
        else
            echo "chr$line" >> "$output_file"
        fi
    # Exclude lines starting with any non-numeric value other than 'MT'
    fi
done

echo "Processing complete. Output saved to $output_file"
