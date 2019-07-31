% ----------------------------------------------------------------------
% Model 1 STS
% ----------------------------------------------------------------------
fprintf("Starting STS for model 1\n")
% Read the whole data set
path_to_data = "../../../Intermediate/Data_whole_tidy_filt.csv";
data = readtable(path_to_data);
cell_data = table2array(data(:, [1, 6, 8]));
clear data 

% Defining the inital values with k1, k2, k3, k4, k5, k6, k7, 
% k8, k9, k10, Suc20 and Glc0, same as for Monolix 
theta_start = [143.6, 37.8, 66, 30, 7.1, 21.3, 3.4, 173.2, ...
    10, 206, 4, 1.6];

% Define the ODE-function 
cd ../Models
ode_fun = @Model1;
cd ../STS

% Function for calculating the cost function 
obj_fun = @calc_cost_STS;

% General parmaeters
n_cells = 125;
n_parameters = 12;
time_int = [0, 1];

% Estimate the parameters using STS 
fprintf("Starting parameter estimation\n");
[param_result, residuals] = perform_STS(cell_data, theta_start, ...
    n_parameters, time_int, n_cells, obj_fun, ode_fun);

% Calculate the predicted values 
fprintf("Calculating the predicted values\n");
start_vec = ones(n_cells, 4); start_vec(:, 4) = 0;
start_vec(:, 3) = param_result(:, end - 2);
start_vec(:, 1) = param_result(:, end - 1);
rate_vec = [param_result(:, 1:end - 3), start_vec(:, 1)];
[ind_fits, pred_fits] = compute_individual_fits(rate_vec, ...
    start_vec, time_int, n_cells, ode_fun);

% Check if results directory exist, else create 
cd Model1
if ~exist("./Result", "dir")
   mkdir("Result")
end

% Save the result to disk 
csvwrite("Result/Parameters_and_cost_func.csv", param_result);
csvwrite("Result/Residuals.csv", residuals);
csvwrite("Result/Predicted_fits.csv", pred_fits);
csvwrite("Result/Individual_fits.csv", ind_fits);
cd ..

% ----------------------------------------------------------------------
% Model 2 STS
% ----------------------------------------------------------------------
fprintf("Starting STS for model 2\n")
% Read the whole data set
path_to_data = "../../../Intermediate/Data_whole_tidy_filt.csv";
data = readtable(path_to_data);
cell_data = table2array(data(:, [1, 6, 8]));
clear data 

% Defining the inital values with k1, k2, k3, k4, k5, k6, k7, 
% k8, k9, k10, Suc20 and Glc0, same as for Monolix 
theta_start = [103.8, 14.1, 45.9, 402.4, 8.5, 80.6, 4.9, 8.8, ...
    1.8, 4.1, 4, 10.3];

% Define the ODE-function 
cd ../Models
ode_fun = @Model2;
cd ../STS

% Function for calculating the cost function 
obj_fun = @calc_cost_STS;

% General parmaeters
n_cells = 125;
n_parameters = 12;
time_int = [0, 1];

% Estimate the parameters using STS 
fprintf("Starting parameter estimation\n");
[param_result, residuals] = perform_STS(cell_data, theta_start, ...
    n_parameters, time_int, n_cells, obj_fun, ode_fun);

% Calculate the predicted values 
fprintf("Calculating the predicted values\n");
start_vec = ones(n_cells, 4); start_vec(:, 4) = 0;
start_vec(:, 3) = param_result(:, end - 2);
start_vec(:, 1) = param_result(:, end - 1);
rate_vec = [param_result(:, 1:end - 3), start_vec(:, 1)];
[ind_fits, pred_fits] = compute_individual_fits(rate_vec, ...
    start_vec, time_int, n_cells, ode_fun);

% Check if results directory exist, else create 
cd Model2
if ~exist("./Result", "dir")
   mkdir("Result")
end

% Save the result to disk 
csvwrite("Result/Parameters_and_cost_func.csv", param_result);
csvwrite("Result/Residuals.csv", residuals);
csvwrite("Result/Predicted_fits.csv", pred_fits);
csvwrite("Result/Individual_fits.csv", ind_fits);
cd ..
