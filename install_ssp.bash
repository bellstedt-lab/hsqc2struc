#!/usr/bin/bash

# Get the path to the conda executable
conda_path=$(which conda)

# If conda is not found, exit
if [ -z "$conda_path" ]; then
    echo "Conda not found!"
    exit 1
fi

# Extract the base directory from the conda path
base_dir=$(dirname $(dirname $conda_path))

# Run initialization script
source "$base_dir/etc/profile.d/conda.sh"

conda create --name hsqc2struc python=3.10 --yes
conda activate hsqc2struc
conda install pip --yes

python3 -m pip install --upgrade -r requirements.txt ../ts_remote_api*.whl ../bruker_nmr_api*.whl --no-input
ln -s ssp.py ../ssp.py

python_path=$(which python)

echo "Installation (hopefully) finished"
echo ""
echo "Your Hsqc2Struc Python executable is located at: $python_path"
echo ""
echo " Please open TopSpin, enter 'set' in the command line, navigate to Python3 Environment, and paste in the above path. Further Information can be found in the Readme.md file."
