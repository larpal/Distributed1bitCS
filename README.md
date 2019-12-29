This repo contains the code to generate the figures in "Analysis of Hard-Thresholding for Distributed Compressed Sensing with One-Bit Measurements" by J. Maly and L. Palzer.

Enclosed are the following files:
 - Distributed_1bit_HT.m - main file
 - joint_HT_1Bit_distributed.m - reconstruction algorithm for distributed 1-bit CS - HT_1Bit.m - reconstruction algorithm for 1-bit CS

To run the code, simply run the script "Distributed_1bit_HT.m" in Matlab. The code will generate many instances of 1-bit compressed sensing problems, run the reconstruction algorithms and plot the results as shown in Sec. 5 of the paper. The default parameters are chosen to run fewer simulation as in the paper so that the execution time is within a few minutes minutes on a Laptop (2 cores). The parameters used in the paper are given in the main file. When running the script, a progress bar indicates the current progress of the program.
