# Scene representations and behaviour

This code was used for the project "" by Karapetian et al. (2022). The code was written in Matlab 2021a and Python 3.6. 

To clone this repository, use

```
git clone https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making.git
```



## 1. Preprocessing of EEG data

* Preprocess the EEG data according to the standard steps, as described in the manuscript:
[Preprocessing_PDM_full_experiment.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/PREPROCESSING/Preprocessing_PDM_full_experiment.m).

This script relies on functions from the [Fieldtrip toolbox](https://www.fieldtriptoolbox.org/download/). 

The output is a Matlab structure ```timelock.mat``` which is used in subsequent analyses.

## 2. Decoding (MVPA)

### 1. Scene decoding

* Decode individual scene identity from subject-level EEG data: [SVM_object_decoding_full_experiment.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/MVPA/First-level/SVM_object_decoding_full_experiment.m)

* Average the decoding results across individuals and plot the results: [plot_decoding_both_tasks.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/MVPA/Plotting/plot_decoding_both_tasks.m)

* Perform a statistical analysis on the averaged data: [fdr_permutation_cluster_1sample_alld.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/STATS/fdr_permutation_cluster_1sample_alld.m)

### 2. Category decoding

* Decode scene category identity from subject-level EEG date: [dth_pseudotrials_SVM_full_experiment_cross_validated.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/DTH/First-level/dth_pseudotrials_SVM_full_experiment_cross_validated.m)

* Average the decoding results across individuals and plot the results: [plot_decoding_both_tasks.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/MVPA/Plotting/plot_decoding_both_tasks.m)

* Perform a statistical analysis on the averaged data: [fdr_permutation_cluster_1sample_alld.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/STATS/fdr_permutation_cluster_1sample_alld.m)

### 3. Multidimensional scaling (MDS)

* Visualize peak decoding results using MDS: [mds_peak_decoding.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/MVPA/Plotting/mds_peak_decoding.m)

## 3. Distance-to-hyperplane analysis

