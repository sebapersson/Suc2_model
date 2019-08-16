% -------------------------------------------------------------------
% VPC validation for model 1
% -------------------------------------------------------------------
% Read the effects
data = readtable(...
    "../Monolix/Model1/Result/populationParameters.txt");
parameter_val = table2array(data(1:end, 2));
clear data

% Extract parameters values
fixed_effects = parameter_val(1:10);
Suc20 = parameter_val(11);
Glc0 = parameter_val(12);
sigma_val = parameter_val(91);

% Read the covariance matrix
cov_mat = csvread(...
    "../../Intermediate/Cov_mat_model1_nlme.csv", 1);
cov_mat = cov_mat(:, 2:end);
% Choose the correct model
cd Models
ode_fun = @Model1;
cd ..

% Run the simulations
n_individuals = 125;
simulated_result = calc_VC_quantiles(cov_mat, ode_fun, 1000, ...
    fixed_effects, Suc20, Glc0, sigma_val, n_individuals, 10);

% Save the resulting csv-file
csvwrite("../../Intermediate/VPC_model1.csv", ...
    simulated_result);

% -------------------------------------------------------------------
% VPC validation for model 2
% -------------------------------------------------------------------
% Read the effects
data = readtable(...
    "../Monolix/Model2/Result/populationParameters.txt");
parameter_val = table2array(data(1:end, 2));
clear data

% Extract parameters values
fixed_effects = parameter_val(1:10);
Suc20 = parameter_val(11);
Glc0 = parameter_val(12);
sigma_val = parameter_val(70);

% Read the covariance matrix
cov_mat = csvread(...
    "../../Intermediate/Cov_mat_model2_nlme.csv", 1);
cov_mat = cov_mat(:, 2:end);
% Choose the correct model
cd Models
ode_fun = @Model2;
cd ..

% Run the simulations
n_individuals = 125;
simulated_result2 = calc_VC_quantiles(cov_mat, ode_fun, 1000, ...
    fixed_effects, Suc20, Glc0, sigma_val, n_individuals, 10);

% Save the resulting csv-file
csvwrite("../../Intermediate/VPC_model2.csv", ...
    simulated_result2);
