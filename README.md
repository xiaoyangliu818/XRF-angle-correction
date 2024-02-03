# XRF-angle-correction
#This is the repo of the code for XRF angle correction
Motivation: get elemental concentration from 3D XRF tomography. Commonly, a known concentration standard is measured at an angle same as the sample to calculate concentration (ug/cm2). In 3D XRF, this standard is used for 2D XRF projections collected at different angles, however, this will generate wrong concentration values for XRF taken at different angles. So we want to correct this error. To do this, we calculated the xrf intensity ratio from Sherman equation.

Simulation of XRF spectrum: I used XMI-MSIM to simulate XRF spectrum.
Construct h5 file for XRF spectrum fitting and quantification

1. Theoretic_cal_20231119.py: script for calculation
2. construct_h5_plot_check.py: construct h5 file that can be analyzed by XRF-MAPS and UprobeX
3.image_analysis.py: histogram of selected region, calculated cumulative concentration probability