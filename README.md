# multi-link
Multi-Link Analysis:Case-Control Brain Network Comparison via Sparse Connectivity Analysis

Imaging neuroscience is currently pushing towards the analysis of the brain from a connectivity perspective, which is unveiling many insights into brain structure and functionality. Neuroscientists are often required to evaluate experimental effects in case-control on specimens with thousands connections. This software is based on an unsupervised machine-learning algorithm, which can capture the multivariate relationships that characterize two distinct groups of connectomes, thus allowing neuroscientists to immediately visualize only the sub-networks that contain information about differences between case and control group. The method exploits recent machine learning techniques which employ sparsity in order to deal with weighted network composed of hundreds of thousands of connections.

Thi code performs the multi-link analysis for case-control comparison. In the example run_examples_MA.m structural connectivity data obtained by the atlasless connectivity tool. This tool can detect the discriminat features (trainTestSLDA function), amd can save the results in trackvis files (writeTracksSel).
All source code unless otherwise specified in the source code are copyrighted to Luca Giancardo and  Alessandro Crimi. This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation.

Remember to cite our paper:

@article{Crimi2016multi,
author = {Alessandro Crimi, Luca Giancardo, Fabio Sanbataro, Alessandro Gozzi, Vittorio Murino, Diego Sona},
title = {Multi-Link Analysis:Case-Control Brain Network Comparison via Sparse Connectivity Analysis},
booktitle = { },
year = { },
month = { }
}

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
