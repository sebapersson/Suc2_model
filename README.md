# About 

This repository contains the code for a project that was done at [Cvijovic Lab](http://www.cvijoviclab.org/) as a part of the course Individual project in mathematics and mathematical statistics (MVE405) at Chalmers University of Technology. The aim of the project was to model recently collected single cell SUC2 expression time series data from *Saccharomyces cerevisiae*, see Figure 1 for the data set. 

To fully replicate the produced result in this project the requirements in the **Requirements for replication of result** section below should be fulfilled. Given this all the tables and figures in the Results-directory can be re-created by running Run_all.sh script in the Code/Shell-directory. When running the Run_all.sh it should be noted that the pathway to Matlab might need to be changed, as the script currently uses absolute path to the Matlab executable. 

## Brief summary of the project 

The motivation behind creating a model to describe the observed data in Figure 1 is that the mechanism behind the observed expression pattern is currently unknown. A powerful tool, and the one used in this project, for testing different possible mechanism is dynamic ODE modelling [1]. If a mechanistic ODE-model, based on a certain hypothesis, fits the data well it's evidence that a proposed mechanism is possible. The ODE-modelling of a single cell data set can be carried by either modelling the average or each individual cell. Modelling of each individual cell should be preferred, as averaging can mask dynamic features as oscillations [2]. To model individual cells two approaches were used and compared, non linear mixed-effects modelling (NLME) and the two stage standard (STS) approach [3, 4]. 

The STS approach and all simulations were implemented in Matlab, while Monolix [5] was used for the NLME approach. To allow the usage of ggplot all the results produced by Matlab and almost all result produced by Monolix were processed and plotted in R. The only Monolix produced text not processed in R was the covariances of the random effects. This is because processing the covariances required string manipulations that were easier to do in Python3 than R. 

## Repository structure

Each directory contains a readme.md file describing the role of that directory. The role of each directory can be summarised as:

**Data**: Contains the original SUC2 data. Upon cloning this project no data files are present in this folder, however the data is available upon request from *niek@chalmers.se*. 
**Docs**: Contains documentation.
**Scratch**: Contains files that can be saftley removed.
**Result**: Contains the figures and files produced by the code.
**Code**: Contains the Python, R, Matlab, Mathematica and Shell-code required to replicate the project. The role of the files is explained in the README_code.md file. 
**Intermediate**: Intermediate files that aren't directly result, but required for producing the final result.

## Requirements for replication of result

This entire repository was created on Ubuntu Linux, and the code should be able to run on any Unix-based system. Below follows a short description of the coding languages and packages used for producing the result. 

R version 3.6.1 was used. The R libraries used were:

* **tidyverse**, version 1.2.1
* **RColorBrewer**, version 1.1.2
* **latex2exp**, version 0.4.0
* **stringr**, version 1.4.0
* **corrplot**, version 0.84
* **reshape2**, version 1.4.3
* **xtable**, version 1.8.3
* **corrr**, version 0.4.0

As only standard functions were used within these libraries the result should be reproducible given other versions. 

Python version 3.7.3 was used. The python packages used were:

* **numpy**, version 1.16.2
* **pandas**, version 0.24.2

As only standard functions were used within these packages the result should be reproducible given other versions. 

Besides R and Python Matlab version 2019R1a and Mathematica version 12.0 were used. The *IdentifiabilityAnalysis*-package used in the Mathematica code is available upon request from *identifiabilityanalysis@fcc.chalmers.se*. 

## References

1. Klipp E, Herwig R, Kowald A, Wierling C, Lehrach H. Systems biology in practice: concepts, implementation and application. John Wiley & Sons; 2008.
2. Cohen-Saidon C, Cohen AA, Sigal A, Liron Y, Alon U. Dynamics and variability of ERK2 response to EGF in individual living cells. Molecular cell. 2009;36(5):885–893.
3. Karlsson M, Janzén DL, Durrieu L, Colman-Lerner A, Kjellsson MC, Cedersund G. Nonlinear mixed-effects modelling for single cell estimation: when, why, and how to use it. BMC systemsbiology. 2015;9(1):52.
4. Almquist J, Bendrioua L, Adiels CB, Goksör M, Hohmann S, Jirstrand M. A nonlinear mixed effects approach for modeling the cell-to-cell variability of Mig1 dynamics in yeast. PloS one. 2015;10(4):e0124050.
5. Monolix version 2019R1. Antony, France: Lixoft SAS; 2019. http://lixoft.com/products/monolix/.
