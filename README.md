This repository contains MATLAB scripts and data files used to generate blank-corrected exometabolite concentrations from six axenic phytoplankton. 

------
First step is to convert Skyline output tables into data tables with nM concentrations.

scripts used: 

Step1_riSkyline_C13.m

Step1_riSkyline_D5.m

------

Next, metabolites known to be problematic were filtered and matrix correction factors were applied.

scripts used:

Step3_filtering_and_matrix_correction.m

------

The corrected data were processed with the following script to choose the SIL-IS type and ion mode:

Step4_combineandsort_Yuting_version.m

-----

Finally, concentration data were corrected with media blanks and the final output table includes nanomolar concentrations, cell-specific concentrations, carbon concentrations, and % contritbution to phytoplankton-excretted DOC.

Step5_normalizedAbundance_two_SILIS.m

-----

The repository is to support a manuscript titled: Characterization of Phytoplankton-Excreted Metabolites Mediating Carbon Flux through the Surface Ocean
Now available in preprint at: https://doi.org/10.1101/2025.11.04.686593

-----
Contact: Yuting Zhu (y2zhu@odu.edu), Assistant Professor, Old Dominion University



