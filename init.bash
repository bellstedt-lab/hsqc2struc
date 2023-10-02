## build setup for SSP usage
# build softlink to parent directory (where topspin looks for scripts)
ln -s ../ssp.py ssp.py
# move to parent directory
cd ..
# create a new conda environment (has to be set in topspin as python 3 env)
conda create --name hsqc2struc --file hsqc2struc/requirements.txt
