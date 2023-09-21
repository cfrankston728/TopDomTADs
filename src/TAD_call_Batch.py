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
    script_content = f"""#!/bin/bash\n#SBATCH --partition gpu\n#SBATCH --output=TopDomTADs_5_%j.out         ### File in which to store job output\n#SBATCH --error=TopDomTADs_5_%j.err          ### File in which to store job error messages\n#SBATCH --cpus-per-task 2\n#SBATCH --mem 100G\n#SBATCH --time 24:00:00\n#SBATCH --job-name TopDomTADs_5_\n#SBATCH --gres=gpu:p100:1\n\n# source activate pyTAD_analysis\n/usr/bin/time -v python /home/groups/CEDAR/franksto/projects/TAD_project/src/TAD_call.py"""

    for param_name, param_value in param_dict.items():
        script_content += f" --{param_name} {param_value}"

    # Write the script content to the script file
    with open(script_filename, 'w') as f:
        f.write(script_content)

# Base Configuration Class


# TESTING
if __name__ == '__main__':
    import subprocess

    window_values = [4,5, 9,10, 19,20, 39,40]
    hic_type_values = ['IN_TACT']#, 'IN_SITU']
    reso_values = [5000, 2000, 1000, 500, 200, 100, 50, 20, 10]#10000, 25000, 50000, 100000, 500000, 1000000]

    all_param_sets = []
    generate_tree(['window', 'hic_type', 'reso'], 0, {}, all_param_sets)

    for i, param_set in enumerate(all_param_sets):
        this_script_filename = f"script_{i}.sh"
        create_script_file(this_script_filename, param_set)

        subprocess.run(['sbatch', this_script_filename])