#!/usr/bin/env python

import numpy as np
import pandas as pd

'''
Function that will find all indices that match a certain character in a string
Args:
string, the string
character, the character
'''
def find_ch(string, character):
    return [i for i, ltr in enumerate(string) if ltr == character]

'''
Function that will given a path to populationParameters.txt will calculate the 
covariance matrix and write it to a csv-file. 
Args:
path_pop_param, the path, from python folder, to the population parameters
path_save, the path to where the file should be saved 
'''
def create_cor_mat(path_pop_param, path_save):
    # Read the population parameters 
    data_pop = pd.read_csv(path_pop_param)
    data_pop_values = data_pop[["parameter", "value"]]
    # Remove the variance sigma 
    data_pop_values = data_pop_values.iloc[0:-1]
    
    # Get the number of parameters, use the fact that no parameter start on omega
    n_parameters = 0
    for parameter in data_pop_values["parameter"]:
        if parameter[0:5] == "omega":
            break 
        n_parameters += 1
        
        # Allocate the correlation matrix
        cov_mat = np.zeros((n_parameters, n_parameters))
        
        # Associate each parameter with a value
        dict_param = {}
    for i in range(n_parameters):
        # Find the first underline
        parameter = data_pop_values.iloc[i, 0]
        pos_underline = parameter.find("_")
        dict_param[parameter[0:pos_underline]] = i
        
    # Fill the covariance matrix
    for i in range(n_parameters, data_pop_values.shape[0]):
        # Get the parameter
        parameter = data_pop_values.iloc[i, 0]
        # The case for a diagonal entry
        if parameter[0:5] == "omega":
            which_param = parameter[parameter.find("_") + 1:]
            index = dict_param[which_param]
            cov_mat[index, index] = data_pop_values.iloc[i, 1] 
            
            # The case for correlation
        elif parameter[0:4] == "corr":
            # Get all underlines 
            index_underline = find_ch(parameter, "_")
            # Get the numbers for the matrix
            i1 = dict_param[parameter[index_underline[0]+1:index_underline[1]]]
            i2 = dict_param[parameter[index_underline[1]+1:]]
            # Calculate the covariance
            cov_mat[i1, i2] = data_pop_values.iloc[i, 1]*cov_mat[i1, i1]*cov_mat[i2, i2]
            cov_mat[i2, i1] = data_pop_values.iloc[i, 1]*cov_mat[i1, i1]*cov_mat[i2, i2]
    
    # Ensure variance along diagonal
    for i in range(n_parameters):
        cov_mat[i, i] *= cov_mat[i, i]
        
    parameter_names = data_pop_values["parameter"][0:n_parameters]
    data_to_save = pd.DataFrame(cov_mat, columns = parameter_names)
    data_to_save.to_csv(path_save)
    
    return 0

# For Model1
path_pop_param = "../Monolix/Model1/Result/populationParameters.txt"
path_save = "../../Intermediate/Cov_mat_model1_nlme.csv"
create_cor_mat(path_pop_param, path_save)

# For Model2
path_pop_param = "../Monolix/Model2/Result/populationParameters.txt"
path_save = "../../Intermediate/Cov_mat_model2_nlme.csv"
create_cor_mat(path_pop_param, path_save)

