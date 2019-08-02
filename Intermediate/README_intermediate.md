# About 

This directory contains intermediate files required to produce the result. 

The four kind of files in this directory for this project is data sets, covariance matrices for the NLME models, data required to create the VPC plots and data files required for creating the simulation plots. 

Below follows a brief description of each kind of file. 

## Data sets

The data set, *Data_whole_tidy_filt.csv* is the combination of data set 1 and data set 2 after the data has been filtered. This is the data set used by all monolix files. 


## VPC-files

The data for the visual predictive checks are created by the matlab-script VPC_check.m, as the reuslt is plotted in R the intermediate files are stored in this folder. 

## Covariance matrices

These files are covariance matrices for the population matrices produced from the parameter estimations in monolix. These matrices are created by a python script (which is required as monolix won't output a correlation matrix), and are required for simulation of data. 

## Simulation-files

These files contains the data produced from simulating the models in Matlab. As the simulation result is plotted in R storage of the intermediate files are stored in this folder. 
