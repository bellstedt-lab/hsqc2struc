#!/usr/bin/bash

cd ..
conda create --name hsqc2struc python=3.10
conda activate hsqc2struc
conda install pip
git clone https://github.com/bellstedt-lab/hsqc2struc
python3 -m pip install hsqc2struc/requirements.txt ts_remote_api*.whl bruker_nmr_api*.whl
ln -s hsqc2struc/ssp.py ssp.py

python_path=$(which python)
echo "Your Hsqc2Struc Python executable is located at: $python_path"
echo "Installation finished. Please open TopSpin, enter 'set' in the commandline, navigate to Python3 Environment and paste the above path. Further Information can be found in the Readme.md file."
