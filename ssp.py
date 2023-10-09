import sys
import os
import pandas as pd
import numpy as np
import shap

import warnings
warnings.filterwarnings("ignore", message=".*The 'nopython' keyword.*")

from catboost import CatBoostRegressor

from bruker.api.topspin import Topspin
from bruker.data.nmr import *

class SecStrucPredictor():
    """predicts the three state secondary structure (helix, sheet and coil) of a protein from an N-HSQC spectra"""
    
    def __init__(self):
        """initializes the input matrices and loads the predictor model"""
        script_directory = os.path.dirname(os.path.abspath(__file__))
        subdir_path = os.path.join(script_directory, "hsqc2struc")

        self.predictor = CatBoostRegressor()
        self.predictor.load_model(subdir_path+"/"+self.get_model())
    
    def get_model(self):
        #8257 BMRB-PDB match, 3514 proteins, NxH binning: 14x10 (helix), 10x8 (sheet), 16x12 (coil) 
        return 'JD_8257_3979_NxH_a14x10_b10x8_c16x12.cbm'
        
    def get_input(self, input_matrices):
        """reads input matrix"""
        self.matrix_14x10 = input_matrices[0]
        self.matrix_10x8 = input_matrices[1]
        self.matrix_16x12 = input_matrices[2]

        return self
    

    def combine_inputs(self):
        """combines and reshapes the input matrices"""
        self.matrix_14x10_1D = self.matrix_14x10.reshape(-1)
        self.matrix_10x8_1D = self.matrix_10x8.reshape(-1)
        self.matrix_16x12_1D = self.matrix_16x12.reshape(-1)

        self.combined_inputs = list(self.matrix_14x10_1D) + list(self.matrix_10x8_1D) + list(self.matrix_16x12_1D)

        return self


    def predict_structure_composition(self):
        self.predictions = self.predictor.predict(self.combined_inputs)

        return self.predictions 

    
    def calc_shap_values(self):
        """calulates shap values for all quadrants for one sample"""

        explainer = shap.TreeExplainer(self.predictor)
        self.shap_values = explainer.shap_values([self.combined_inputs], [self.predictions])

        return self


    def build_shap_spectra(self):
        """builds spectra based on calculated shap values"""
        
        count_matrix = [self.matrix_14x10_1D, self.matrix_10x8_1D, self.matrix_16x12_1D]
        fig, axs = plt.subplots(figsize=((15,10)), nrows=3, ncols=3)
        fig.suptitle("Shap Values")

        for i, (sec_struc_type, bin_type) in enumerate(zip(["Helix", "Sheet", "Coil"], ["14x10", "10x8", "16x12"])):
            axs[0][i].set_title(bin_type)
            axs[2][i].set_xlabel("H-shift in ppm", fontsize=13)
            axs[i,0].set_ylabel(sec_struc_type + "\n N-Shift in ppm", fontsize=13)
            spectra_14x10 = self.shap_values[i].ravel()[:140].reshape(14,10)
            
            spectra_10x8 = self.shap_values[i].ravel()[140:220].reshape(10,8)
            spectra_16x12 = self.shap_values[i].ravel()[220:].reshape(16,12)
         
            spectra = [spectra_14x10, spectra_10x8, spectra_16x12]

            #for ii, (spec, x, y, H_scale, N_scale) in enumerate(zip(spectra, [10,10,8], [20,26,10], [0.5,0.5,0.625], [2.5,1.9230769230769231,5])):
            for ii, (spec, x, y, H_scale, N_scale) in enumerate(zip(spectra, [10,8,12], [14,10,16], [0.5,0.625,0.4166666666666667], [3.5714285714285716,5,3.125])):

                im = axs[i][ii].imshow(spec, cmap="coolwarm", vmax=0.032, vmin=-0.032)
                cbar = fig.colorbar(im, extend="both")
                axs[i][ii].set_xticks(np.arange(x))
                axs[i][ii].set_yticks(np.arange(y))

                axs[i][ii].set_xlabel("H-Shift in ppm")
                axs[i][ii].set_ylabel("N-Shift in ppm")

                
                axs[i][ii].set_xticks(np.arange(x), [str(round(i,1)) for i in np.arange(6,11,H_scale)], rotation=45)
                axs[i][ii].set_yticks(np.arange(y), [str(round(i,1)) for i in np.arange(90,140,N_scale)])
                
                axs[i][ii].set_xlim(x-0.5,-0.5)

                N_shifts = sorted(list(range(spec.shape[0]))*spec.shape[1])
                H_shifts = list(range(spec.shape[1]))*spec.shape[0]

                axs[i][ii].scatter(x=H_shifts, y=N_shifts, s=count_matrix[ii], color="black", alpha=0.3)
        plt.tight_layout()
        plt.show()

    
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
        script_directory = os.path.dirname(os.path.abspath(__file__))
        subdir_path = os.path.join(script_directory, "hsqc2struc")
        peak_list = pd.read_csv(subdir_path+"/BMRB_peak_list_ubiquitin.csv")
        N_shifts = peak_list["Y_shift"].to_list()
        H_shifts = peak_list["X_shift"].to_list()

    else:
        top = Topspin()
        dp = top.getDataProvider()

        hsqc = dp.getCurrentDataset()

        peak_list = hsqc.getPeakList()

        if len(peak_list) == 0:
            raise ValueError("NO PEAKS SELECTED! Use the pp command or select the peaks manually!")

        
        H_shifts = []
        N_shifts = []

        max_intensity = np.max([peak["intensity"] for peak in peak_list])
        
        for peak in peak_list:
            rel_intensity = peak["intensity"] / max_intensity
            if rel_intensity > 0: # negative intensity ==> sidechain NH(2) ### here you can adjust sensitivity!!!!!!!!!!!!!!!!!!!!!!
                H_shifts.append(peak["position"][0])
                N_shifts.append(peak["position"][1])


    predictor = SecStrucPredictor() 

    input_matrices = []
    
    H_shift_max = 11
    H_shift_min = 6 

    N_shift_max = 140 
    N_shift_min = 90

    # NxH binning: 14x10 (helix), 10x8 (sheet), 16x12 (coil); values are inverted here:
    for H_num_1D_grid, N_num_1D_grid in [(10,14),(8,10),(12,16)]:
    
        H_binsize = (H_shift_max - H_shift_min) / H_num_1D_grid
        N_binsize = (N_shift_max - N_shift_min) / N_num_1D_grid
    

        binned_H_shifts = binning(H_shifts, H_binsize, H_shift_min, H_num_1D_grid)
        binned_N_shifts = binning(N_shifts, N_binsize, N_shift_min, N_num_1D_grid)

        count_peaks_matrixes = generate_count_peaks_matrix(binned_H_shifts, binned_N_shifts, H_num_1D_grid, N_num_1D_grid)

        input_matrices.append(count_peaks_matrixes)

    
    predictor.get_input(input_matrices)
    predictor.combine_inputs()
    prediction = predictor.predict_structure_composition()

    print(f" Secondary Structure Prediction \n -------------------------------- \n Helix: {round(prediction[2],3)*100:.1f}%  \n Sheet: {round(prediction[1],3)*100:.1f}% \n Coil: {round(prediction[0],3)*100:.1f}%")
    
    if "shap" in sys.argv:
        import matplotlib.pyplot as plt
        import shap
        predictor.calc_shap_values()
        predictor.build_shap_spectra()

    if len(sys.argv) > 1 and sys.argv[1] == "test":
        if str(25.5) == str(round(prediction[2],3)*100)[:4] and str(16.9) ==str(round(prediction[1],3)*100)[:4] and str(57.6) == str(round(prediction[0],3)*100)[:4]:
            print("\n \n \t \t <<< PREDICTION TEST PASSED >>> ")
        else:
            print("\n \n \t \t >>> PREDICTION TEST FAILED!!!!!!! <<< ")
            print("Possibly you are currently using a different model?") 
