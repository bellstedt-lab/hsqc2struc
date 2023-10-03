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
- Please be sure to activate the conda environment before executing the script
- Run `python3 hsqc2struc.py test` to check if everything If the test is passed, continue.
- Run `python3 hsqc2struc.py BMRB_peak_list_ubiquitin.csv` or `python3 hsqc2struc.py peak_list_ubiquitin.csv` to predict the secondary structure content of Ubiquitin using our CatBoost model. Feel free to adapt the example script and run the command with your own HSQC peak list (*.csv). Be sure to have the same (header) structure in your own peak list that corresponds to either of the provided examples input files. Also, be careful not (!) to include any peaks belonging to side chain NH2 groups. Alternatively, you can collect HSQC-DEPT spectra and phase the spectrum so that the NH2 groups have negative intensities (which are then ignored by the predictor).

### Deinstallation
1. Deactivate the conda environment (if necessary) with `conda deactivate`
2. Delete the conda environment with `conda remove --name hsqc2struc_cmdline --all`
3. Remove the hsqc2struc directory with `rm -Rf hsqc2struc`

## Integration into Bruker TopSpin
### Requirements
- TopSpin >= 4.2
- conda
  
### Installation
1. Open a terminal and go to /opt/topspin4.X.Y/python/examples/ (on Linux). If you do not have this directory, please have a look at the Bruker Website / Documentation to reinstall the TopSpin Python Interface.
2. Create new conda environment with `conda create --name hsqc2struc python=3.10`
3. Activate the new environment with `conda activate hsqc2struc`
4. Install pip for the handling of the required packages: `conda install pip`
5. Execute `git clone https://github.com/bellstedt-lab/hsqc2struc` (but stay in the topspin example directory afterward)
6. Install the required Python packages:
   ```
   python3 -m pip install -r hsqc2struc/requirements.txt ts_remote_api*.whl bruker_nmr_api*.whl
   ```
8. Create a symbolic link so that TopSpin can find our predictor script: `ln -s hsqc2struc/ssp.py ssp.py`
3. Activate the environment with `conda activate hsqc2struc` and identify the location of the python executable with `which python` (copy the path into the clipboard or write it down)
5. Open TopSpin and type "set" in the command line, click on the "Change" Button next to "Select Python 3+ Environment" and paste or enter the location of the Python executable just identified.
6. Still in the preferences window, click on "Change" next to "Manage TopSpin Network Interface", enter your Admin password and Start the Network Interface. For your convenience, you can check the autostart option below.

### Usage
- Collect/Load Protein 1H,15N-HSQC spectrum
- Perform peak-picking (enter "pp" in TopSpin's command line), and ensure that the sidechain NH2 groups are not (!) picked. Alternatively, you can collect a 1H,15N-HSQC-DEPT spectrum and phase the spectrum so that the NH2 groups have negative intensities (which are then ignored by the predictor).
- Type "xpy3 ssp" in the TopSpin command line. After a few seconds, a new pop-up window should appear with the predicted secondary structure composition in percent.
- If you also want to include the SHAP plots, type "xpy3 ssp shap".

## Deinstallation
1. Delete the conda environment with `conda remove --name hsqc2struc --all`
2. Delete the hsqc2struc folder from the topspin subdirectory ( `rm -Rf /opt/topspin4.X.Y/python/examples/hsqc2struc` )
3. Remove the symbolic link created in the example directory with `rm /opt/topspin4.X.Y/python/examples/ssp.py`
 
## More information about the contributing research group
[![GitHub Logo](https://www.bellstedt-lab.ch/images/logo_blab_400px.png)](https://www.bellstedt-lab.ch)
