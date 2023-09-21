import argparse
import os
import csv
import subprocess

##########################################################################################################################
# DEFINITION:

# bigWigAvgOutScript.py operates using the following parameters:

    # base_dir.............A base directory from which to access other directories.
    
    # bigWigLinksFile------A path to a .csv file containing links to .bigWig files to be processed.
    #                      Each row is of the form "target, fold_change_over_control_link, signal_p_link."
    #                      The links should be either saved .bigWig files or urls to a database.
    
    # TAD_bed..............A path to a bed file containing TAD boundary regions.
    
    # Control_bed----------A path to a matching control for the TAD boundary file.
    
    # bigWigAvg............A path to the bigWigAverageOverBed function to be passed in a bash command to the subprocesser.
    
    # header---------------True if bigWigLinksFile has a header, False otherwise.
    
    # run_control..........1 if control files are to be constructed, 0 otherwise.

# The function is to write out bed files indicating the average bigWig values over the TAD boundary regions
##########################################################################################################################
# IN PROGRESS:
# Reconfigure script to operate on files in parallel. []

# Define parameters for operation:
base_dir = "/home/groups/CEDAR/franksto/projects/TAD_project/" # This is just for efficiency.

# Create file with Histone targets and bigWig file links
bigWigLinksFile = base_dir + "data/bigWigLinks_Compartment_Score.csv"
header = True

# Obtain the bigWigAverageOverBed reference to pass to subprocess.
bigWigAvg = base_dir + "scripts/bigWigAverageOverBed"

# For switching between left and right adjacent windows. If windows are centered, set side='' (the empty string)
side = "/Right"

TAD_beds_folder = base_dir + "src/TAD_bd_Adjacent_Windows"+side+"/"

# Parse command-line arguments
parser = argparse.ArgumentParser(description='Process TAD boundary data and compute enrichments at TAD boundaries.')
parser.add_argument('--TAD_bed', type=str, default='',
                    help='TAD boundary bed file.')

args = parser.parse_args()

TAD_bed = args.TAD_bed
#"GM12878_TopDom_window_5_r10kbp_TAD_boundaries_IN_SITU"

Control_bed = base_dir + "data/Prior_TAD_Boundary_Beds/GM12878_TopDom_window_5_r10kbp_TAD_boundaries_CONTROL" # Leave this alone

run_control = 0 # 0 for False, 1 for True

# For each TAD_bed file in the TAD_beds_folder directory:
#   Obtain all enrichments 

# Run this in parallel across the different files? Does this matter? (Probably not, the function runs quickly)
# 8/29/2023: I think this is starting to matter.
with open(bigWigLinksFile, 'r') as file:
    csv_reader = csv.reader(file)
    if header:
        next(csv_reader)   
        
    for row in csv_reader:
        
        target_head = base_dir + "src/Compartments_bigWigAvgOuts"+side+"/"
        this_folder_path = target_head + TAD_bed

        if not os.path.exists(this_folder_path):
            os.makedirs(this_folder_path)
        
        target_out = this_folder_path + "/" + row[0]
        target_control = target_head + "bigWigAvgBedControls/control_" + row[0]
        
        signal_type = ["_fold_change_over_control", "_signal_p"]
        
        bed_type = [TAD_beds_folder+TAD_bed+'.bed', Control_bed]
        target_type = [target_out, target_control]
        
        for k in [1,2]:
            if not(not row[k]):
                for j in range(run_control + 1):
                    this_output = target_type[j] + signal_type[k-1] + ".bed"
                    this_command = f"{bigWigAvg} {row[k]} {bed_type[j]} {this_output}"
                    subprocess.run(this_command.split(' '))
