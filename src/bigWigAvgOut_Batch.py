#!/usr/bin/env python3

# Description
# This file will contain all the configurations to be made for each run.
#
import numpy as np
import os
import datetime
import torch
from torch import nn

print("Imported Config libraries\n")
###################################################################################################

# For switching between left and right adjacent windows. If windows are centered, set side='' (the empty string)
side = "/Right"

def generate_tree(params, current_index, current_params, all_param_sets):
    if current_index == len(params):
        all_param_sets.append(current_params.copy())
        return

    current_param_name = params[current_index]
    current_param_values = globals()[current_param_name + '_values']

    for value in current_param_values:
        current_params[current_param_name] = value
        generate_tree(params, current_index + 1, current_params, all_param_sets)


def create_script_file(script_filename, param_dict):
    # Prepare the script content with your modular script command and the parameter set
    script_content = f"""#!/bin/bash\n#SBATCH --partition gpu\n#SBATCH --output=OE_bigWiGAvgOutScript%j.out         ### File in which to store job output\n#SBATCH --error=OE_bigWiGAvgOutScript%j.err          ### File in which to store job error messages\n#SBATCH --cpus-per-task 2\n#SBATCH --mem 100G\n#SBATCH --time 48:00:00\n#SBATCH --job-name OE_bigWiGAvgOutScript\n#SBATCH --gres=gpu:p100:1\n\n# source activate pyTAD_analysis\n/usr/bin/time -v python /home/groups/CEDAR/franksto/projects/TAD_project/src/bigWigAvgOutScript.py"""

    for param_name, param_value in param_dict.items():
        script_content += f" --{param_name} {param_value}"

    # Write the script content to the script file
    with open(script_filename, 'w') as f:
        f.write(script_content)

# Base Configuration Class

# TESTING
if __name__ == '__main__':
    import subprocess

    folder_path = '/home/groups/CEDAR/franksto/projects/TAD_project/src/TAD_bd_Adjacent_Windows' + side
    file_names = os.listdir(folder_path)

    # Remove the '.bed' extension and create a list of values
    TAD_bed_values = [file_name.split('.')[0] for file_name in file_names if file_name.endswith('.bed')]

    all_param_sets = []
    generate_tree(['TAD_bed'], 0, {}, all_param_sets)

    for i, param_set in enumerate(all_param_sets):
        this_script_filename = f"script_{i}.sh"
        create_script_file(this_script_filename, param_set)

        subprocess.run(['sbatch', this_script_filename])

    