#!/usr/bin/env python3
import sys

def process_data(input_file):
    # Read data from input file
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # Dictionary to store maximum values for each key in the first column
    max_values = {}

    # Process each line
    for line in lines:
        parts = line.split()
        key = parts[0]
        value = float(parts[2])

        # Update maximum value for the key
        if key in max_values:
            if value > max_values[key][1]:
                max_values[key] = (line.strip(), value)
        else:
            max_values[key] = (line.strip(), value)

    # Output filtered lines
    filtered_lines = [max_values[key][0] for key in max_values]

    # Print the filtered lines
    for line in filtered_lines:
        print(line)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python process_data.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    process_data(input_file)
