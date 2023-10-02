# HSQC-to-STRUC
[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIT)

HSQC2STRUC is a Python 3 script to predict the secondary structure composition of proteins by using their unassigned 1H,15N-HSQC spectra. This program is intended to be used as a command line tool in TopSpin 4.2. Furthermore, the CatBoost model can be integrated into any custom script. However, due to simplicity, we recommend using our web service at https://hsqc2struc.bellstedt-lab.ch, which relies on exactly the same machine learning model.

## Requirements
- TopSpin 4.2
- conda

## Installation
Go to /opt/topspin4.2.0/python/examples/ (on MacOS) and use 
```git clone https://github.com/joaldi2208/SSP.git .```
then
```bash init.bash```
now it should take some time to create a new conda environment and install all dependencies.
If successful open TopSpin and type ```set``` in the command line and select the just created environement. 
To find the path of the conda environment you need to enter ```conda activate Sec-struc-pred``` followed by ```which python````. Copy the given path and set it as new Python environment.

## Usage
Load Protein HSQC spectrum, perform peak-picking, and ensure that the sidechain NH2 groups are not (!) picked. Type "xpy3 ssp" in the TopSpin command line. After a few seconds, a new pop-up window should appear with the predicted secondary structure composition in percent. If you also want to include the SHAP plots, type "xpy3 ssp shap".
 


## More information about the supporting research groups
[![GitHub Logo](https://www.bellstedt-lab.ch/images/logo_blab_400px.png)](https://www.bellstedt-lab.ch)
