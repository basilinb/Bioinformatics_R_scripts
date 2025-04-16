import pandas as pd
import os
import shutil
import hashlib
import csv

import pandas as pd
import os
import shutil
import hashlib
import csv


# Function to compute MD5 hash of a file
def compute_md5(file_path):
    hasher = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):  # Read in chunks
            hasher.update(chunk)
    return hasher.hexdigest()

# Function to process files in a folder
def process_files(source_folder, destination_folder, extensions, copy = False):
    md5_data = []
    for root, _, files in os.walk(source_folder):
        for filename in files:
            if any(filename.endswith(ext) for ext in extensions):
                # skip if file name contains "empty" or "Undetermined"
                if "empty" in filename or "Undetermined" in filename:
                    continue
                src_path = os.path.join(root, filename)
                dst_path = os.path.join(destination_folder, filename)

                if not os.path.exists(dst_path):
                    if copy:
                        shutil.copy(src_path, dst_path)
                        print(f"Copied: {src_path} -> {dst_path}")
                    else:
                        os.symlink(src_path, dst_path)
                        print(f"Linked: {src_path} -> {dst_path}")
                else:
                    print(f"Skipped (already exists): {dst_path}")

                md5_hash = compute_md5(src_path)
                md5_data.append((filename, md5_hash))

    return md5_data

# Run for the v2 dataset (Novogene) --------------------------------------------
# Define source and destination folders
destination_folder = "/nfs/bioinformatics/workspace/geo_submissions/spribitzer_amishra_P503-1"

# Collect MD5 checksums
md5_data = []
# CSV file for MD5 checksums
md5_csv_path = os.path.join(destination_folder, "md5_checksums.csv")

# create destination folder
os.makedirs(destination_folder, exist_ok=True)

# Get fastqs
# Since the fastq files are in a different location that might not be accessible to Stephanie, we need to copy them
md5s = process_files(
    "/nfs/spribitzer/itn/bioinformatics/pipeline/Illumina/230302_VH00126_262_AAAVJGHHV/Unaligned",
    destination_folder,
    [".fastq", ".fastq.gz"],
    copy=True
)
md5_data.extend(md5s)

md5s = process_files(
    "/nfs/spribitzer/itn/bioinformatics/pipeline/Illumina/230118_VH00126_253_AAAV5VCHV/Unaligned",
    destination_folder,
    [".fastq", ".fastq.gz"],
    copy=True
)
md5_data.extend(md5s)

# Get matrix h5 files
md5s = process_files(
    "/nfs/bioinformatics/pipeline/Illumina/MultiFlowcellProjects/230302_VH00126_262_AAAVJGHHV/Project_P503-1Processed_bri_230320/peaks",
    destination_folder,
    [".bed", "narrowPeak"]
)
md5_data.extend(md5s)


# Write MD5 checksums to CSV
with open(md5_csv_path, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Filename", "MD5 Checksum"])  # Header
    writer.writerows(md5_data)

print(f"MD5 checksums saved to {md5_csv_path}")

