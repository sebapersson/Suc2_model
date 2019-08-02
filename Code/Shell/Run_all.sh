#!/bin/bash

# Running this script will reproduce all the result in the result folder,
# with other words the aim with this file is to make the entire project
# as reproducible as possible. Since it is impossible to run some of the scripts
# from the command line the program will urge the user to run certain monolix
# files with certain starting values.
# The script has to be run from 

# -------------------------------------------------------------------------
# Start of functions 
# -------------------------------------------------------------------------
# Function that will check if a directory exists. If the directory doesn't
# exist the function will create it. Note that input order matters 
# Args:
#     $1 Directory name 
check_if_dir_exists ()
{
    # Create directory if it doesn't exist
    if [ ! -d $1 ]; then
	mkdir $1
    fi
}

# -------------------------------------------------------------------------
# End of functions 
# -------------------------------------------------------------------------

# Check that the script is run from Shell directory
currentDir=${PWD##*/}

if [ ! $currentDir == "Shell" ]; then
    >&2 echo "The script most be run from Code/Shell directory"
    exit 1
fi 

# -------------------------------------------------------------------------
# Filter the data and create figures of the data set 
# -------------------------------------------------------------------------
# Create the appropriate directories in the result folder 
cd ../../Result
check_if_dir_exists Figures
check_if_dir_exists Files

cd Figures
check_if_dir_exists Data_set

# Run the r-script and create the files and figures 
cd ../../Code/R
echo "Creating figures of the data set"
Rscript ./Explore_data_sets.R 2> /dev/null 

# -------------------------------------------------------------------------
# Run the monolix code 
# -------------------------------------------------------------------------
## Model 1
cd ../Monolix/Model1/

# If the Monolix code hasn't been run prompt the user to run it.
if [ ! -f "Model1/populationParameters.txt" ]; then
    echo "Model 1 hasn't been run in monolix"
    echo "Run Model1.mlxtran file in monolix 2019R1"
    echo "The following initial values should be used:"
    echo "k1 = 143.6"
    echo "k2 = 37.8"
    echo "k3 = 66"
    echo "k4 = 30"
    echo "k5 = 7.1"
    echo "k6 = 143.6"
    echo "k7 = 143.6"
    echo "k8 = 143.6"
    echo "k9 = 143.6"
    echo "k10 = 143.6"
    echo "Suc20 = 4.0"
    echo "Glc0 = 1.6"
    echo "Number of iterations for SAEM should be 2000"
    echo "and a full correlation model should be fitted"
    exit 1
fi

# In case a full model has been fitted copy everything to result folder
check_if_dir_exists Result
cp -r Model1/* Result/

## Model2 
cd ../Model2
# If the Monolix code hasn't been run prompt the user to run it.
if [ ! -f "Model2/populationParameters.txt" ]; then
    echo "Model 2 hasn't been run in monolix"
    echo "Run Model2.mlxtran file in monolix 2019R1"
    echo "The following initial values should be used:"
    echo "k1 = 103.8"
    echo "k2 = 14.1"
    echo "k3 = 45.9"
    echo "k4 = 402.4"
    echo "k5 = 8.5"
    echo "k6 = 80.6"
    echo "k7 = 4.9"
    echo "k8 = 8.8"
    echo "k9 = 1.8"
    echo "k10 = 4"
    echo "Suc20 = 4.1"
    echo "Glc0 = 10.3"
    echo "Number of iterations for SAEM should be 500"
    echo "and all parameters except k10 and Suc20 should"
    echo "be included in the correlation model"
    exit 1
fi

# In case a full model has been fitted copy everything to result folder
check_if_dir_exists Result
cp -r Model2/* Result/

# Move back to head directory
cd ../../..

echo "Both models have been run in Monolix"

# -------------------------------------------------------------------------
# Diagnostics NLME models 
# -------------------------------------------------------------------------
# Running the diagnostics for the NLME models
cd Result/Figures
check_if_dir_exists Model1_nlme
check_if_dir_exists Model2_nlme
cd ../..

# Run the R-script for creating the diagnostic plots (except simulation and VPC)
echo "Creating diagnostic plots for nlme model 1 and 2"
cd Code/R
Rscript Diagnose_nlme.R 2> /dev/null 
echo "Done diagnostic plots for nlme model 1 and 2"

# Move back to root
cd ../..

# -------------------------------------------------------------------------
# Run the VPC analysis
# -------------------------------------------------------------------------
# Create the covariance matrices
cd Code/Python
./Create_cov_mat.py
cd ../..

# Don't run if VPC file already exist (they are expansive to compute)
cd Intermediate/

if [ ! -f VPC_model1.csv ] || [ ! -f VPC_model2.csv ]; then
    echo "Running VPC"
    cd ../Code/Matlab
    # The whole path is used to ensure that R2019a version of matlab is used 
    /usr/local/MATLAB/R2019a/bin/matlab -nodisplay -nosplash -nodesktop -r "run('Perform_VPC.m');exit;" | tail -n +11

    cd ../R/
    Rscript Plot_VPC.R 2> /dev/null
    cd ../../Intermediate/
    
fi 

echo "Done creating VPC plots"
cd ..

# -------------------------------------------------------------------------
# Run the STS estimation 
# -------------------------------------------------------------------------
cd Code/Matlab/STS

# Don't run analysis if result already exist (since it is expansive)
if [ ! -f Model1/Result/Residuals.csv ] || [ ! -f Model2/Result/Residuals.csv ]; then
    # Run STS in matlab
    /usr/local/MATLAB/R2019a/bin/matlab -nodisplay -nosplash -nodesktop -r "run('STS_approach.m');exit;" | tail -n +11
fi

# Move back to project root 
cd ../../..

echo "Done with STS analysis"

# -------------------------------------------------------------------------
# Process the STS result 
# -------------------------------------------------------------------------
echo "Creating diagnosis plots for STS"

cd Result/Figures 
check_if_dir_exists Model1_STS
check_if_dir_exists Model2_STS

cd ../../Code/R/
Rscript Diagnose_STS.R 2> /dev/null

echo "Done creating diagnosis plots for STS"

cd ../..

# -------------------------------------------------------------------------
# Run the simulations 
# -------------------------------------------------------------------------
cd Intermediate

# Don't run if files already exist as this is expansive
if [ ! -Simulation_model1.csv ] || [ ! -f Simulation_model2.csv ]; then
    cd ../Code/Matlab
    echo "Running simulations for nlme models"
    /usr/local/MATLAB/R2019a/bin/matlab -nodisplay -nosplash -nodesktop -r "run('Simulate_cells_nlme.m');exit;" | tail -n +11
    cd ../../Intermediate
fi

# Don't run if files already exist as this is expansive
if [ ! -Simulation_model1_STS.csv ] || [ ! -f Simulation_model2_STS.csv ]; then
    cd ../Code/Matlab
    echo "Running simulations for STS models"
    /usr/local/MATLAB/R2019a/bin/matlab -nodisplay -nosplash -nodesktop -r "run('Simulate_cells_STS.m');exit;" | tail -n +11
    cd ../../Intermediate
fi 

echo "Done simulating cells"
cd ..

# -------------------------------------------------------------------------
# Plot the simulations 
# -------------------------------------------------------------------------
cd Code/R
echo "Plotting the simulations"
Rscript Plot_simulations.R 2> /dev/null
echo "Done plotting the simulations"

cd ../../

# -------------------------------------------------------------------------
# Create the extrapolated data for model 2 
# -------------------------------------------------------------------------
cd Intermediate
# Don't run if files already exist as this is expansive
if [ ! -f Simulation_model2_extraplate.csv ] || [ ! -f Simulation_model2_extra_glc_1_20.csv ] || [ ! -f Simulation_model2_extra_glc_1_2.csv ]; then
    cd ../Code/Matlab
    echo "Running extrapolation for model 2"
    /usr/local/MATLAB/R2019a/bin/matlab -nodisplay -nosplash -nodesktop -r "run('Extrapolate_simulations.m');exit;" | tail -n +11
    cd ../../Intermediate
fi 

cd ..
echo "Done with extrapolation calculations for model 2"

# -------------------------------------------------------------------------
# Plot the extrapolated data 
# -------------------------------------------------------------------------
cd Result/Figures
check_if_dir_exists Model2_extrapolated

cd ../../Code/R

echo "Plotting the extrapolations"
Rscript Plot_extrapolations.R 2> /dev/null
echo "Done plotting the extrapolations"

cd ../../

echo ""
echo "The analysis ran with zero errors. All created figures and files can now be"
echo "found in the result folder."

exit 0 
