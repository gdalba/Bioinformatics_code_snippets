import os
import sys

print('-----------------------------------')
print('CTENOTATION V1 - 06/19/2023')
print('modify_fsatq_names.py')
print('Made by Gabriel Dall\'Alba')
print('For Browne data only')
print('provide input_folder')
print('-----------------------------------')

# Receive sysargs for input folder
input_folder = sys.argv[1]

# Get a list of all .fastq.gz files
fastq_files = [f for f in os.listdir(input_folder) if f.endswith('.fastq.gz')]

for fastq_file in fastq_files:
    if "_pair" in fastq_file:
        new_file_name = fastq_file.replace("R1_pair", "1").replace("R2_pair", "2")
        os.rename(os.path.join(input_folder, fastq_file), os.path.join(input_folder, new_file_name))
