# DeCom

1. OVERVIEW

We introduce the DeCom toolbox which offers a set of Matlab functions for non-stationary multi-component signal decomposition. It includes routines for visualizing the time-frequency distribution of oscillatory signals, extracting time-frequency ridges for well-separated components, estimating instantaneous frequencies, amplitudes, and waveforms. The routines work for both the (generalized) mode decomposition model and the multiresolution mode decomposition model. These models and decomposition techniques admit both theoretical analysis and efficient numerical implementation.

Applications:

Geophysics:  seismic wave field separation and ground-roll removal.

Materials science:  atomic crystal image analysis, grain boundary and local defects identification, elastic deformation estimation.

Art:  Canvas painting analysis for art forensics, canvas removal for paiting conservation

MechanicalEngineering: fault detection

Astropyhysics: LIGO signal analysis

2. REFERENCE

The folder SynLab is a collection of Matlab and MEX routines which implements 1D and 2D synchrosqueezed transforms proposed in [1]-[4]. It contains numerical examples in [2][5][6].

The folder GeneralModeDecom is a collection of Matlab routines for generalized mode decomposition and multiresolution mode decomposition studied in [1], [9-11].

[1] H. Yang. Synchrosqueezed wave packet transforms and diffeomorphism based spectral analysis for 1d general mode decompositions. Applied and Computational Harmonic Analysis, 39(1):33 – 66, 2015.

[2] H. Yang, J. Lu, and L. Ying. Crystal image analysis using 2D synchrosqueezed transforms. Multiscale Modeling & Simulation, 13(4):1542–1572, 2015.

[3] H. Yang and L. Ying. Synchrosqueezed wave packet transform for 2d mode decompo- sition. SIAM Journal on Imaging Sciences, 6(4):1979–2009, 2013.

[4] H. Yang and L. Ying. Synchrosqueezed curvelet transform for two-dimensional mode decomposition. SIAM Journal on Mathematical Analysis, 46(3):2052–2083, 2014.

[5] H. Yang.  Statistical analysis of synchrosqueezed transforms, To appear, Applied and Computational Harmonic Analysis, 2017.

[6] H. Yang, J. Lu, W. P. Brown, I. Daubechies, and L. Ying, Quantitative Canvas Weave Analysis Using 2D Synchrosqueezed Transforms. IEEE Signal Processing Magazine, Special Issue on Art Investigations, 2015.

[7] J. Lu, B. Wirth and H. Yang. Compbining 2d synchrosqueezed wave packet transforms with optimization for crystal image analysis. Journal of the Mechanics and Physics of Solids, Volume 89, April 2016, Pages 194-210.

[8] H. Yang. Oscillatory data analysis and fast algorithms for integral operators, Ph.D. thesis, Stanford University, 2015.

[9] J. Xu, H. Yang, and I. Daubechies, Recursive Diffeomorphism-Based Regression for Shape Functions. SIAM Journal on Mathematical Analysis, to appear.

[10] H. Yang, Multiresolution Mode Decomposition for Adaptive Time Series Analysis. Submitted.

[11] G. Tang and H. Yang, A Fast Algorithm for Multiresolution Mode Decomposition. Submitted.


3. INSTALLING DeCom

Run the file setpath.m first. It will automatically add all the MATLAB codes to your MATLAB path and compile all MEX files. After this, you can run all demo codes to see how to use this tool box.

4. COPY RIGHT

DeCom is copyright reserved. For further information, please contact
Haizhao Yang at matyh@nus.edu.sg

