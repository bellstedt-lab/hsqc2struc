# HSQC-to-STRUC
[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIT)

HSQC2STRUC is a Python 3 script to predict the secondary structure composition of proteins by using their unassigned 1H,15N-HSQC spectra. This program is intended to be used as a command line tool in TopSpin 4.2. Furthermore, the CatBoost model per se can be integrated into any custom script. However, due to simplicity, we recommend using our web service at https://hsqc2struc.bellstedt-lab.ch, which relies on exactly the same machine-learning model.

## Command line script
### Requirements
- python3
- conda

### Installation
1. Create new conda environment with `conda create --name hsqc2struc_cmdline`
2. Activate the new environment with `conda activate hsqc2struc_cmdline`
3. Install pip for the handling of the required packages: `conda install pip`
4. Install the required python packages: `python3 -m pip install numpy pandas catboost`
5. Execute `git clone https://github.com/bellstedt-lab/hsqc2struc`
6. Change into the new hsqc2struc directory and now you are ready to use our model (see next section)

### Usage
Run ```python3 hsqc2struc.py peak_list_ubiquitin.csv``` and the secondary structure content of Ubiquitin will be predicted using our CatBoost model. Feel free to adapt the example script and to run the command with your own HSQC peak list (*.csv). Be careful not (!) to include any peaks belonging to side chain NH2 groups. 

## Integration into Bruker TopSpin
### Requirements
- TopSpin 4.2
- conda
  
### Installation
Go to /opt/topspin4.2.0/python/examples/ (on MacOS & Linux) and execute ```git clone https://github.com/bellstedt-lab/hsqc2struc```. Change into the new hsqc2struc directory and run ```bash init.bash``` to install the requirements. 
now it should take some time to create a new conda environment and install all dependencies. If you observe an "xcrun error: invalid active developer path" error on MacOS, please run "xcode-select --install".
If successful, open TopSpin and type ```set``` in the command line and select the just created environment. 
To find the path of the conda environment, you need to enter ```conda activate hsqc2struc``` followed by ```which python```. Copy the given path and set it as a new Python environment.

### Usage
Load Protein HSQC spectrum, perform peak-picking, and ensure that the sidechain NH2 groups are not (!) picked. Type "xpy3 ssp" in the TopSpin command line. After a few seconds, a new pop-up window should appear with the predicted secondary structure composition in percent. If you also want to include the SHAP plots, type "xpy3 ssp shap".
 


## More information about the contributing research group
[![GitHub Logo](https://www.bellstedt-lab.ch/images/logo_blab_400px.png)](https://www.bellstedt-lab.ch)
