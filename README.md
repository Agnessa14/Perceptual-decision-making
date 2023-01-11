# Scene representations and behaviour

This code was used for the project "" by Karapetian et al. (2022). The code was written in Matlab 2021a and Python 3.6. 

To clone this repository, use

```
git clone https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making.git
```



## 1. Preprocessing of EEG data

The first step is to preprocess the EEG data according to the standard steps, as described in the manuscript, using this script:
[Preprocessing_PDM_full_experiment.m](https://github.com/Neural-Dynamics-of-Visual-Cognition-FUB/Perceptual-decision-making/blob/master/ANALYSIS_Full_experiment/PREPROCESSING/Preprocessing_PDM_full_experiment.m).

This script relies on functions from the [Fieldtrip toolbox](https://www.fieldtriptoolbox.org/download/). 

The output is a Matlab structure ```timelock.mat``` which is used in subsequent analyses.

