% ----------------------------------------------------------------------- 
% Model 1 
% ----------------------------------------------------------------------- 
% Read the effects 
data = readtable(...
    "../Monolix/Model1/Result/populationParameters.txt");
parameter_val = table2array(data(1:end, 2));
clear data

% Extract parameters values 
fixed_effects = parameter_val(1:10);
Suc20 = parameter_val(11);
Glc0 = parameter_val(12);

% Read the covariance matrix 
cov_mat = csvread(...
    "../../Intermediate/Cov_mat_model1_nlme.csv", 1);
cov_mat = cov_mat(:, 2:end);

% Choose the correct model 
cd Models
ode_fun = @Model1;
cd ..

% Set parameters for cells to simulate 
n_cells_simulate = 50000;
n_states = 4;
n_parameters = 12;
time_int = [0, 1];
n_point_simulate = 100;

% Run the simulations, this model uses the same structure as model 7
simulated_data = simulate_single_cells_nlme(n_cells_simulate, n_states, ...
    n_parameters, time_int, fixed_effects, Suc20, Glc0, cov_mat, ...
    ode_fun, n_point_simulate);

csvwrite("../../Intermediate/Simulation_model1.csv", ...
    simulated_data);

% ----------------------------------------------------------------------- 
% Model 2
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
n_cells_simulate = 5000;
n_states = 4;
n_parameters = 12;
time_int = [0, 1];
n_point_simulate = 100;

% Run the simulations, this model uses the same structure as model 7
simulated_data = simulate_single_cells_nlme(n_cells_simulate, n_states, ...
    n_parameters, time_int, fixed_effects, Suc20, Glc0, cov_mat, ...
    ode_fun, n_point_simulate);

csvwrite("../../Intermediate/Simulation_model2.csv", ...
    simulated_data);

% ----------------------------------------------------------------------- 
% Model 2 short delay model 
% ----------------------------------------------------------------------- 
% Read the effects 
data = readtable(...
    "../Monolix/Model2_test/Result/populationParameters.txt");
parameter_val = table2array(data(1:end, 2));
clear data

% Extract parameters values 
fixed_effects = parameter_val(1:10);
Suc20 = parameter_val(11);
Glc0 = parameter_val(12);

% Read the covariance matrix 
cov_mat = csvread(...
    "../../Intermediate/Cov_mat_model2_short_del_nlme.csv", 1);
cov_mat = cov_mat(:, 2:end);

% Choose the correct model 
cd Models
ode_fun = @Model2_short_del;
cd ..

% Set parameters for cells to simulate 
n_cells_simulate = 50000;
n_states = 4;
n_parameters = 12;
time_int = [0, 1];
n_point_simulate = 100;

% Run the simulations, this model uses the same structure as model 7
simulated_data = simulate_single_cells_nlme(n_cells_simulate, n_states, ...
    n_parameters, time_int, fixed_effects, Suc20, Glc0, cov_mat, ...
    ode_fun, n_point_simulate);

csvwrite("../../Intermediate/Simulation_model2_short_del.csv", ...
    simulated_data);