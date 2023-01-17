# Brain-behaviour relationship in scene categorization

This code was used for the project "Empirically linking and computationally modelling the brain-behavior relationship for human scene categorization" by Karapetian et al. (2023). The code was written in Matlab 2021a, Python 3.6 and Tensorflow 2.7.  

To clone this repository, use

```
git clone https://github.com/Agnessa14/Perceptual-decision-making.git
```



## 1. Preprocessing of EEG data

* Preprocess the EEG data according to the standard steps, as described in the manuscript:
[Preprocessing_PDM_full_experiment.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/PREPROCESSING/Preprocessing_PDM_full_experiment.m).

This script relies on functions from the [Fieldtrip toolbox](https://www.fieldtriptoolbox.org/download/). 

The output is a Matlab structure ```timelock.mat``` which is used in subsequent analyses.

## 2. Decoding (MVPA)

Decoding and distance-to-hyperplane analyses were performed with the SVM from the [libsvm toolbox:](https://www.csie.ntu.edu.tw/~cjlin/libsvm/).

### 1. Scene decoding

* Decode individual scene identity from subject-level EEG data: [SVM_object_decoding_full_experiment.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/MVPA/First-level/SVM_object_decoding_full_experiment.m) 

* Average the decoding results across individuals and plot the results: [plot_decoding_both_tasks.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/MVPA/Plotting/plot_decoding_both_tasks.m)

### 2. Category decoding

* Decode scene category from subject-level EEG date: [dth_pseudotrials_SVM_full_experiment_cross_validated.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/DTH/First-level/dth_pseudotrials_SVM_full_experiment_cross_validated.m)

* Average the decoding results across individuals and plot the results: [plot_decoding_both_tasks.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/MVPA/Plotting/plot_decoding_both_tasks.m)

### 3. Searchlight analysis 

* Perform decoding in channel space: 
  - subject-level:
    - scene identity: [SVM_object_decoding_full_experiment_searchlight.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/OTHER/SVM_object_decoding_full_experiment_searchlight.m) 
    - scene category: [dth_pseudotrials_SVM_full_experiment_cross_validated_sl.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/OTHER/dth_pseudotrials_SVM_full_experiment_cross_validated_sl.m) 
  - average: [all_subjects_decoding_sl.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/MVPA/Average/all_subjects_decoding_sl.m)
  - plot: [plot_topography_searchlight_decoding.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/OTHER/plot_topography_searchlight_decoding.m)

### 4. Multidimensional scaling (MDS)

* Visualize peak decoding results using MDS: [mds_peak_decoding.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/MVPA/Plotting/mds_peak_decoding.m)

* Visualize peak distance-to-hyperplane results: [mds_peak_dth_correlation.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/OTHER/mds_peak_dth_correlation.m)

## 3. Distance-to-hyperplane analysis

### 1. Obtain the subject-level distances to the natural/man-made hyperplane

* Using the same script as for category decoding, calculate the distances: [dth_pseudotrials_SVM_full_experiment_cross_validated.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/DTH/First-level/dth_pseudotrials_SVM_full_experiment_cross_validated.m)

### 2. Correlate distances and reaction times 

* Correlate the subject-level distances with median RTs, average over subjects, and plot: [dth_all_distances_median_RT.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/DTH/Average/dth_all_distances_median_RT.m)

* Plot results from both tasks together: [plot_dth_both_tasks.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/DTH/Plotting/plot_dth_both_tasks.m)

### 3. Searchlight analysis

* Perform distance-to-hyperplane analysis in channel space: 
  - subject-level: [dth_pseudotrials_SVM_full_experiment_cross_validated_sl.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/OTHER/dth_pseudotrials_SVM_full_experiment_cross_validated_sl.m)
  - average: [dth_all_distances_median_RT_searchlight.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/DTH/Average/dth_all_distances_median_RT_searchlight.m)
  - plot: [plot_topography_searchlight_dth_ak.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/OTHER/plot_topography_searchlight_dth_ak.m)

## 4. Modelling with RCNN 

### 1. Fine-tuning RCNN 

* Fine-tune the BLNet (from Spoerer et al., original repo [here](https://github.com/cjspoerer/rcnn-sat)) - initially trained on ecoset (Mehrer et al., 2017): [rnn_dth_finetune_rnn_larger_dataset.ipynb](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/DNN/rnn_dth_finetune_rnn_larger_dataset.ipynb)

* Load the weights from the fine-tuned model: 

```
import urllib
_, msg = urllib.request.urlretrieve(
    'https://osf.io/', 'model_02.11_2_weights.h5')
print(msg)
```

### 2. Feature and RT extraction

* Extract features from hidden layers and RTs from the prediction layer: [rnn_dth_collect_activations.ipynb](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/DNN/rnn_dth_collect_activations.ipynb)

### 3. Representational similarity analysis (RSA) with EEG
#### 1. EEG representational dissimilarity matrix (RDM) construction

* Create subject-level RDMs at every time point using 1-Pearson's r: [pearson_object_decoding_full_experiment.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/MVPA/First-level/pearson_object_decoding_full_experiment.m)

* Calculate the noise ceiling: [noise_ceiling_rsa.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/OTHER/noise_ceiling_rsa.m)

#### 2. RCNN RDM construction

* Create RDMs using 1-Pearson's r at the selected layers and time steps: [RNN_test_Input_RDM_Correlations.ipynb](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/DNN/RNN_test_Input_RDM_Correlations.ipynb)

#### 3. Correlation between RDMs

* Correlate (Spearman's) subject-level EEG and RCNN RDMs at every time point, selected layer and network time step: [rsa_rnn_eeg_subject_by_subject_averaged_timesteps.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/OTHER/rsa_rnn_eeg_subject_by_subject_averaged_timesteps.m)

### 4. Correlation between RTs

* Correlate (Pearson's) every subject's RTs with RCNN RTs and plot: [plot_RT_rnn_vs_human.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/OTHER/plot_RT_rnn_vs_human.m)

* Calculate the noise ceiling: [noise_ceiling_RT.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/OTHER/noise_ceiling_RT.m)

### 5. Distance-to-hyperplane analysis between EEG and RCNN RTs

* Correlate subject-level EEG distances with RCNN RTs: [rnn_dth_all_distances_median_RT.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/DTH/Average/rnn_dth_all_distances_median_RT.m)

## 5. Statistics examples

* Permutation test and FDR correction: [fdr_permutation_cluster_1sample_alld.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/STATS/fdr_permutation_cluster_1sample_alld.m)

* Bootstrapping: [bootstrap_peak_latency.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/STATS/bootstrap_peak_latency.m)

