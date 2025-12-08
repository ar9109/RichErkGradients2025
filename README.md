# RichErkGradients2025

This repository contains code, source data, and data sets for reproducing the analysis and figures from:  Decaying and expanding Erk gradients process memory of skeletal size during zebrafish fin regeneration
Ashley Rich, Ziqi Lu, Alessandro De Simone, Lucas Garcia, Jacqueline Janssen, Kazunori
Ando, Jianhong Ou, Massimo Vergassola, Kenneth D. Poss, & Stefano Di Talia

The imaging processing code found in 'segmentation_pipeline.m' will perform image processing to quantify Erk activity and Geminin expression from images of regenerating zebrafish rays.

Copyright (C)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

The code requires MATLAB (MathWorks).
The code runs in MATLAB_R2023b. 
The code runs on a Titan Desktop Workstation running Windows 11 Pro (Version 24H2, OS build 26100.6584). Processor: AMD Ryzen Threadripper PRO 5955WX 16-Cores (4.00 GHz). Installed RAM     256 GB.

# INSTALLATION AND USAGE - Imaging Processing Code

A set of example source data files are included within this Github repo in the ‘raw’ folder.

The image processing pipeline uses TGMM (Amat et al. 2014) to segment nuclei, which is too large to be included in this Github repo. The version used for this manuscript can be accessed at the following link: https://www.dropbox.com/scl/fo/lxgj9wf2rdrzsnvi7vy38/AH71h_sksydVsYdhNDrULiY?rlkey=2zwyihyigeddqp2q7it14wukf&st=2olmkwrm&dl=0. 

To run 'segmentation_pipeline.m', the user should download from the Github repo ‘segmentation_pipeline.m’, 'raw’, ‘tgmm_template’, and ‘functions_052924’. The user should also download ‘TGMM_Supplementary_Software_1_0’ from the link above. The ‘.m’ file and the 4 directories should be collected into a single directory, for example, ‘processDirectory’.

The user can run the image processing code by editing line 9 of ‘segmentation_pipeline.m’ 
to set ‘workplace’ to appropriate directory on user’s local machine (i.e. ‘processDirectory’).
MS Windows users have to adapt paths to MS Windows syntax.
 
MATLAB Toolboxes required are:
- Image Processing Toolbox (version 9.5)
- Statistics and Machine Learning Toolbox (version 11.0)
- Financial Toolbox (version 5.8)
- Curve Fitting Toolbox (3.5.4)

The code requires external MATLAB functions:
loadtiff.m, saveastiff.m, weightedcov.m, xml2struct.m, nanconv.m,
brewermap.m
Those functions be downloaded from MathWorks. We attach them to the code together with
their licenses.

The script can run altogether or each section sequentially.
In some instances, user input would be required.

TGMM segmentation step must be run on MS Windows.

Data sample runtime: ~10 minutes.

Expected outcome from each step are located in 'ExpectedResults' directory. Final outputs of imaging processing pipeline are found in 'ExpectedResults/ERK_activity_map' and 'ExpectedResults/results'.

# Installation and Usage - Figure Code
To run code for generating figures, load the appropriate '.m' file. (For example, to generate plots shown in Figure 3, load 'fig3.m'. Relevant data sets for Figure 3 are located in 'Fig3' directory.) '.m' files for each figure may call functions that can be found in 'functions_052924'.
