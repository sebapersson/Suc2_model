% ---------------------------------------------------------------------
% Model 1 STS
% ---------------------------------------------------------------------

% Read the matrix used to generate random effects
param_data = csvread(...
    "../../Result/Files/STS_model1_cov_mean.csv", 1);
cov_mat = param_data(1:end-1, :);
mean_vec = param_data(end, :);
% Choose the correct model 
cd Models
ode_fun = @Model1;
cd ..


% Simulate the data
n_cells_simulate = 50000;
n_states = 4;
time_int = [0, 1];

% Run the simulations
simulated_data = simulate_single_cells_STS(n_cells_simulate, ...
    n_states, mean_vec, cov_mat, time_int, ode_fun);

% Save the result to disk 
csvwrite("../../Intermediate/Simulation_model1_STS.csv", ...
    simulated_data);

% ---------------------------------------------------------------------
% Model 2 STS
% ---------------------------------------------------------------------

% Read the matrix used to generate random effects
param_data = csvread(...
    "../../Result/Files/STS_model2_cov_mean.csv", 1);
cov_mat = param_data(1:end-1, :);
mean_vec = param_data(end, :);

% Choose the correct model 
cd Models
ode_fun = @Model2;
cd ..

% Simulate the data
n_cells_simulate = 50000;
n_states = 4;
time_int = [0, 1];

% Run the simulations, 
simulated_data = simulate_single_cells_STS(n_cells_simulate, ...
    n_states, mean_vec, cov_mat, time_int, ode_fun);

% Save the result to disk 
csvwrite("../../Intermediate/Simulation_model2_STS.csv", ...
    simulated_data);

% ---------------------------------------------------------------------
% Model 2 STS short del 
% ---------------------------------------------------------------------

% Read the matrix used to generate random effects
param_data = csvread(...
    "../../Result/Files/STS_model2_short_del_cov_mean.csv", 1);
cov_mat = param_data(1:end-1, :);
mean_vec = param_data(end, :);
% Choose the correct model 
cd Models
ode_fun = @Model2_short_del;
cd ..

% Simulate the data
n_cells_simulate = 10000;
n_states = 4;
time_int = [0, 1];

% Run the simulations, 
simulated_data = simulate_single_cells_STS(n_cells_simulate, ...
    n_states, mean_vec, cov_mat, time_int, ode_fun);

% Save the result to disk 
csvwrite("../../Intermediate/Simulation_model2_short_del_STS.csv", ...
    simulated_data);

% ---------------------------------------------------------------------
% Model 1 STS alt cov 
% ---------------------------------------------------------------------

% Read the matrix used to generate random effects
param_data = csvread(...
    "../../Result/Files/STS_model1_cov_mean.csv", 1);
cov_mat = param_data(1:end-1, :);
mean_vec = param_data(end, :);
% Choose the correct model 
cd Models
ode_fun = @Model1;
cd ..

% Read the covariance matrix 
cov_mat = csvread(...
    "../../Intermediate/Cov_mat_model1_nlme.csv", 1);
cov_mat = cov_mat(:, 2:end);

% Simulate the data
n_cells_simulate = 50000;
n_states = 4;
time_int = [0, 1];

% Run the simulations
simulated_data = simulate_single_cells_STS(n_cells_simulate, ...
    n_states, mean_vec, cov_mat, time_int, ode_fun);

% Save the result to disk 
csvwrite("../../Intermediate/Simulation_model1_alt_cov_STS.csv", ...
    simulated_data);

% ---------------------------------------------------------------------
% Model 2 STS using the covariance matrix from NLME 
% ---------------------------------------------------------------------

% Read the matrix used to generate random effects
param_data = csvread(...
    "../../Result/Files/STS_model2_cov_mean.csv", 1);
cov_mat = param_data(1:end-1, :);
mean_vec = param_data(end, :);

% Choose the correct model 
cd Models
ode_fun = @Model2;
cd ..

% Read the covariance matrix 
cov_mat = csvread(...
    "../../Intermediate/Cov_mat_model2_nlme.csv", 1);
cov_mat = cov_mat(:, 2:end);

% Simulate the data
n_cells_simulate = 50000;
n_states = 4;
time_int = [0, 1];

% Run the simulations, 
simulated_data = simulate_single_cells_STS(n_cells_simulate, ...
    n_states, mean_vec, cov_mat, time_int, ode_fun);

% Save the result to disk 
csvwrite("../../Intermediate/Simulation_model2_alt_cov_STS.csv", ...
    simulated_data);

