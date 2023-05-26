import sys
import os
import pandas as pd
import numpy as np
from catboost import CatBoostRegressor

from bruker.api.topspin import Topspin
from bruker.data.nmr import *

class SecStrucPredictor():
    """predicts the three state secondary structure (helix, sheet and coil) of a protein from an N-HSQC spectra"""
    
    def __init__(self):
        """initializes the input matrices and loads the predictor model"""
        self.predictor = CatBoostRegressor()
        self.predictor.load_model("/opt/topspin4.2.0/python/examples/SSP/Model.cbm")
    

    def get_input(self, input_matrices):
        """reads input matrix"""
        self.matrix_20x10 = input_matrices[0]
        self.matrix_26x10 = input_matrices[1]
        self.matrix_10x8 = input_matrices[2]

        return self
    

    def combine_inputs(self):
        """combines and reshapes the input matrices"""
        self.matrix_20x10_1D = self.matrix_20x10.reshape(-1)
        self.matrix_26x10_1D = self.matrix_26x10.reshape(-1)
        self.matrix_10x8_1D = self.matrix_10x8.reshape(-1)

        self.combined_inputs = list(self.matrix_20x10_1D) + list(self.matrix_26x10_1D) + list(self.matrix_10x8_1D)

        return self


    def predict_structure_composition(self):
        self.predictions = self.predictor.predict(self.combined_inputs)

        return self.predictions 




    
def binning(shifts, binsize, shift_min, num_1D_grid):
    """gives the indexes for certain bins by float dividing the shift value by the binsize"""
    bin_indexes = []
    for ppm_value in shifts:
    
        bin_index = int((ppm_value - shift_min) // binsize)
        if bin_index > (num_1D_grid - 1): # to large for set grid limits
            continue # ignore peak
            # bin_index = num_1D_grid - 1 # append to largest quadrant
        elif bin_index < 0: # to small for set grid limits
            continue # ignore peak
            # bin_index = 0 # append to smallest quadrant
        
        bin_indexes.append(bin_index)


    return bin_indexes


def get_shifts(chemical_shifts, atomtype):
    """returns a list of with all shifts of a certain atomtype. X stands for H and Y for N"""
    shifts = []

    for NMR_data in chemical_shifts.values():
        shifts_one_protein = NMR_data[f"{atomtype}_shift"].to_numpy()
        shifts.append(shifts_one_protein)

    return shifts
    

def generate_count_peaks_matrix(binned_H_shifts, binned_N_shifts, H_num_1D_grid, N_num_1D_grid):
    """counts peaks in a certain bin in the defined 2D grid"""
    grid_2D = (N_num_1D_grid, H_num_1D_grid)
    count_peaks_matrix = np.zeros(grid_2D)
    
    for H_shift_bin, N_shift_bin in zip(binned_H_shifts, binned_N_shifts):
        count_peaks_matrix[N_shift_bin, H_shift_bin] += 1

    return count_peaks_matrix
    

if __name__ == "__main__":

    if len(sys.argv) > 1 and sys.argv[1] == "test":
        print("*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+")
        print("Ubiquitin prediction from 5387 backbone spectrum of the BMRB")
        print("*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+")
        peak_list = pd.read_csv("/opt/topspin4.2.0/python/examples/SSP/BMRB_peak_list_ubiquitin.csv")
        N_shifts = peak_list["Y_shift"].to_list()
        H_shifts = peak_list["X_shift"].to_list()

    elif len(sys.argv) > 1 and sys.argv[1].endswith(".csv"):
        print("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ")
        print(f"Spectrum prediction from {os.getcwd() + '/' + sys.argv[1]}")
        print("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ")
        peak_list = pd.read_csv(sys.argv[1])
        N_shifts = peak_list["x(F1) [ppm]"].to_list()
        H_shifts = peak_list["x(F2) [ppm]"].to_list()

    else:
        top = Topspin()
        dp = top.getDataProvider()

        hsqc = dp.getCurrentDataset()
        peak_list = hsqc.getPeakList()

        
        H_shifts = []
        N_shifts = []

        max_intensity = np.max([peak["intensity"] for peak in peak_list])
        
        for peak in peak_list:
            rel_intensity = peak["intensity"] / max_intensity
            if rel_intensity > 0.18: # negative intensity ==> sidechain NH(2) ### here you can adjust sensitivity!!!!!!!!!!!!!!!!!!!!!!
                H_shifts.append(peak["position"][0])
                N_shifts.append(peak["position"][1])


    predictor = SecStrucPredictor() 

    input_matrices = []
    
    H_shift_max = 11
    H_shift_min = 6 

    N_shift_max = 140 
    N_shift_min = 90
    
    for H_num_1D_grid, N_num_1D_grid in [(10,20),(10,26),(8,10)]:
    
        H_binsize = (H_shift_max - H_shift_min) / H_num_1D_grid
        N_binsize = (N_shift_max - N_shift_min) / N_num_1D_grid
    

        binned_H_shifts = binning(H_shifts, H_binsize, H_shift_min, H_num_1D_grid)
        binned_N_shifts = binning(N_shifts, N_binsize, N_shift_min, N_num_1D_grid)

        count_peaks_matrixes = generate_count_peaks_matrix(binned_H_shifts, binned_N_shifts, H_num_1D_grid, N_num_1D_grid)

        input_matrices.append(count_peaks_matrixes)

    
    predictor.get_input(input_matrices)
    predictor.combine_inputs()
    prediction = predictor.predict_structure_composition()
    print(f" Secondary Structure Prediction \n -------------------------------- \n Coil: {round(prediction[0],3)*100}% \n Sheet: {round(prediction[1],3)*100}% \n Helix: {round(prediction[2],3)*100}%")
    
