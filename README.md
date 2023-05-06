# SSP
[![License](https://img.shields.io/badge/License-MIT%202.0-blue.svg)](https://opensource.org/licenses/MIT)

![GitHub Logo](https://github.com/joaldi2208/SSP/blob/main/Logo.jpg?raw=true)
SSP stands for Secondary Structure Predictor and is a python 3 script to predict the secondary structure composition of proteins by using their N-HSQC spectra. This program is intended to be used as a command line tool in TopSpin 4.2.

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
To find the path of the conda environment you need to enter ```conda activate Sec-struc-pred``` followed by ```which python````. Copy the given path and set it as new python environment. Also activate the networking???.


## Usage
Simply type xpy3 ssp in the TopSpin 4.2 command line, while having an HSQC spectra open at the same time. After a few seconds a new pop up window should appear with the predicted secondary structure composition in percent.
pictures should explain it in detail. 

## Further development
This program was part of a master thesis and will not be actively developed by the author in the future. More about the research and programs can be found in the following repository or ...
Since the the NMR databases are growing very fast updating the model could be useful with some time (Last update Nov 2022). Therefore use the workflow described in the mentioned repository (not intended for non programmers).

## More information about the supporting research groups
[![GitHub Logo](https://github.com/Kohulan/DECIMER-Image-to-SMILES/blob/master/assets/CheminfGit.png?raw=true)](https://cheminf.uni-jena.de)
Logo of Peter Bellstedt's group from Zürich University
