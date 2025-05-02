import os
import re

# Define the folder containing the files
data_folder = "MATLAB/data"

# Ensure the folder exists
if not os.path.exists(data_folder):
    print(f"Folder '{data_folder}' does not exist.")
    exit()

# Regular expression to match the file pattern
file_pattern = re.compile(r"corona_(\d+)_angle_([\d.]+)_alpha_([\d.]+)\.csv")

# Iterate through all files in the folder
for filename in os.listdir(data_folder):
    match = file_pattern.match(filename)
    if match:
        corona_value, angle_value, alpha_value = match.groups()
        # Check if corona_value is an integer (no decimal point)
        if corona_value.isdigit():
            # Create the new filename
            new_filename = f"corona_{corona_value}.00_angle_{angle_value}_alpha_{alpha_value}.csv"
            # Rename the file
            old_path = os.path.join(data_folder, filename)
            new_path = os.path.join(data_folder, new_filename)
            os.rename(old_path, new_path)
            print(f"Renamed: {filename} -> {new_filename}")