Imaging neuroscience is currently pushing towards the analysis of the brain from a connectivity perspective, which is unveiling many insights into brain structure and functionality. Neuroscientists are often required to evaluate experimental effects in case-control on specimens with thousands connections. This software is based on an unsupervised machine-learning algorithm, which can capture the multivariate relationships that characterize two distinct groups of connectomes, thus allowing neuroscientists to immediately visualize only the sub-networks that contain information about differences between case and control group. The method exploits recent machine learning techniques which employ sparsity in order to deal with weighted network composed of hundreds of thousands of connections.

Thi code performs the multi-link analysis for case-control comparison. In the example run_examples_MA.m structural connectivity data obtained by the atlasless connectivity tool. This tool can detect the discriminat features (trainTestSLDA function), amd can save the results in trackvis files (writeTracksSel).
With the code you can find the DTI mice data
All source code unless otherwise specified in the source code are copyrighted to Luca Giancardo and  Alessandro Crimi. This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation.

Remember to cite our paper:

@article{Crimi2018multi,
author = {Alessandro Crimi, Luca Giancardo, Fabio Sanbataro, Alessandro Gozzi, Vittorio Murino, Diego Sona},
title = {Multi-Link Analysis:Case-Control Brain Network Comparison via Sparse Connectivity Analysis},
booktitle = { Nature Scientific Reports },
year = {2019 },
month = {January }
}
https://doi.org/10.1038/s41598-018-37300-4

![Mice Differences](https://raw.githubusercontent.com/alecrimi/multi-link/master/Fig2.jpg)

The trainTestSLDA.m scripts has these minimum input and output parameters:
[sldaFeat, sldaWeight, sldaWeightStd] = trainTestSLDA( featMat, classVec );
Input:
     - featMat: a feature matrix containing the highdimensional features for all samples for two different classes
     - classVec: a vector defining the class for each samples (either 0 or 1)

Output:
     - sldaFeat: identified relevant features according to the specified parameters
     - sldaWeight: weights for the identified features 
     - sldaWeightStd: standard deviation for the selected features.

The script also plot an occurence histogram for the relevant fatures
To save the features as Trackvis files, I added the scripts from John Colby (johncolby@ucla.edu) trk_write_new.m
and other external files are in the folder lib

The atlasless script to generate connectivity maps from a fixed grid is in the folder atlasless_connectivity

The Trackvis  mice data are in the folder atlasless_connectivity/dti_res_orig
The unprocessed DTI mice data are on Figshare: https://figshare.com/collections/BTBR_vs_Control_for_MultiLink_Analysis/4350803
