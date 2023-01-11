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

### 3. Searchlight analysis 



### 4. Multidimensional scaling (MDS)

* Visualize peak decoding results using MDS: [mds_peak_decoding.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/MVPA/Plotting/mds_peak_decoding.m)

## 3. Distance-to-hyperplane analysis

### 1. Obtain the subject-level distances to the natural/man-made hyperplane

* Using the same script as for category decoding, obtain the distances: [dth_pseudotrials_SVM_full_experiment_cross_validated.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/DTH/First-level/dth_pseudotrials_SVM_full_experiment_cross_validated.m)

### 2. Correlate distances and reaction times 

* Correlate the subject-level distances with median RTs, average over subjects, and plot: [dth_all_distances_median_RT.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/DTH/Average/dth_all_distances_median_RT.m)

### 3. Searchlight analysis

* Perform distance-to-hyperplane analysis in electrode space: 

## 4. Modelling with RCNN 

* Fine-tune the BLNet (from Spoerer et al., original repo [here](https://github.com/cjspoerer/rcnn-sat)) - initially trained on ecoset (Mehrer et al., 2017): 

### 1. Fine-tuning RCNN 
### 2. Feature and RT extraction
### 3. Representational similarity analysis (RSA) with EEG
#### 1. EEG representational dissimilarity matrix (RDM) construction
#### 2. RCNN RDM construction
#### 3. Correlation between RDMs
### 4. Correlation between RTs
### 5. Distance-to-hyperplane analysis between RCNN and human RTs
