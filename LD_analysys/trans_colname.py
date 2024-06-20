import sys

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 3:
    print("Usage: python script.py input_file.txt output_file.txt")
    sys.exit(1)

# Get input and output file paths from command-line arguments
input_file_path = sys.argv[1]
output_file_path = sys.argv[2]

# Open the input file
with open(input_file_path, "r") as input_file:
    # Open a new file for writing the output
    with open(output_file_path, "w") as output_file:
        # Iterate through each line in the input file
        for line in input_file:
            # Split the line into columns based on whitespace
            columns = line.split()

            # Check if the first column contains "SV" and the second column contains "SNP"
            if "SV" in columns[0] and "SNP" in columns[1]:
                # Write the first column (with "SV"), the second column (with "SNP"), and the third column to the output file
                output_file.write(columns[0] + " " + columns[1] + " " + columns[2] + "\n")
            if "SNP" in columns[0] and "SV" in columns[1]:
                output_file.write(columns[1] + " " + columns[0] + " " + columns[2] + "\n")



print("Output has been written to", output_file_path)
