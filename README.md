# wild-binary-segmentation-2.0
Pre-CRAN code to accompany the paper P. Fryzlewicz (2020) "Detecting possibly frequent change-points: Wild Binary Segmentation 2 and steepest-drop model selection", Journal of the Korean Statistical Society (with discussion), to appear.

Files:

- WBS2_SDLL_github_v*.R - code for i.i.d. Gaussian noise

- WBS2_SDLL_robust.R - robust version, for possibly non-Gaussian and/or heterogeneous noise

- utils.R - some functions needed by both of the above

- WBS2_SDLL_comparison_with_other_methods.R - code to reproduce the simulation results from the paper

Instructions: to use the code, source utils.R, and then one or both of WBS2_SDLL_github_v*.R and WBS2_SDLL_robust.R, as desired.
