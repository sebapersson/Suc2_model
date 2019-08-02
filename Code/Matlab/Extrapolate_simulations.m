% ----------------------------------------------------------------------- 
% Model 2 increasing the time span
% ----------------------------------------------------------------------- 
% Read the effects 
data = readtable(...
    "../Monolix/Model2/Result/populationParameters.txt");
parameter_val = table2array(data(1:end, 2));
clear data

% Extract parameters values 
fixed_effects = parameter_val(1:10);
Suc20 = parameter_val(11);
Glc0 = parameter_val(12);

% Read the covariance matrix 
cov_mat = csvread(...
    "../../Intermediate/Cov_mat_model2_nlme.csv", 1);
cov_mat = cov_mat(:, 2:end);

% Choose the correct model 
cd Models
ode_fun = @Model2;
cd ..

% Set parameters for cells to simulate 
n_cells_simulate = 50000;
n_states = 4;
n_parameters = 12;
time_int = [0, 2];
n_point_simulate = 300;

% Run the simulations
simulated_data = simulate_single_cells_nlme(n_cells_simulate, n_states, ...
    n_parameters, time_int, fixed_effects, Suc20, Glc0, cov_mat, ...
    ode_fun, n_point_simulate);

csvwrite("../../Intermediate/Simulation_model2_extrapolate.csv", ...
    simulated_data);

% ----------------------------------------------------------------------- 
% Model 2 change glucose in from 1 / 40 to 1/ 20 
% ----------------------------------------------------------------------- 
% Read the effects 
data = readtable(...
    "../Monolix/Model2/Result/populationParameters.txt");
parameter_val = table2array(data(1:end, 2));
clear data

% Extract parameters values 
fixed_effects = parameter_val(1:10);
Suc20 = parameter_val(11);
Glc0 = parameter_val(12);

% Read the covariance matrix 
cov_mat = csvread(...
    "../../Intermediate/Cov_mat_model2_nlme.csv", 1);
cov_mat = cov_mat(:, 2:end);

% Choose the correct model 
cd Models
ode_fun = @Model2;
cd ..

% Set parameters for cells to simulate 
n_cells_simulate = 50000;
n_states = 4;
n_parameters = 12;
time_int = [0, 1];
n_point_simulate = 300;
rate_factor_in = 1 / 20;

% Run the simulations
simulated_data = simulate_single_cells_nlme(n_cells_simulate, n_states, ...
    n_parameters, time_int, fixed_effects, Suc20, Glc0, cov_mat, ...
    ode_fun, n_point_simulate, rate_factor_in);

csvwrite("../../Intermediate/Simulation_model2_extra_glc_1_20.csv", ...
    simulated_data);

% ----------------------------------------------------------------------- 
% Model 2 change glucose in from 1 / 40 to 1/ 2 
% ----------------------------------------------------------------------- 
% Read the effects 
data = readtable(...
    "../Monolix/Model2/Result/populationParameters.txt");
parameter_val = table2array(data(1:end, 2));
clear data

% Extract parameters values 
fixed_effects = parameter_val(1:10);
Suc20 = parameter_val(11);
Glc0 = parameter_val(12);

% Read the covariance matrix 
cov_mat = csvread(...
    "../../Intermediate/Cov_mat_model2_nlme.csv", 1);
cov_mat = cov_mat(:, 2:end);

% Choose the correct model 
cd Models
ode_fun = @Model2;
cd ..

% Set parameters for cells to simulate 
n_cells_simulate = 50000;
n_states = 4;
n_parameters = 12;
time_int = [0, 1];
n_point_simulate = 300;
rate_factor_in = 1 / 2;

% Run the simulations
simulated_data = simulate_single_cells_nlme(n_cells_simulate, n_states, ...
    n_parameters, time_int, fixed_effects, Suc20, Glc0, cov_mat, ...
    ode_fun, n_point_simulate, rate_factor_in);

csvwrite("../../Intermediate/Simulation_model2_extra_glc_1_2.csv", ...
    simulated_data);
