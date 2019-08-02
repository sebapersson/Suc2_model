function simulation_result = simulate_single_cells_nlme(n_cells_simulate, ...
    n_states, n_parameters, time_int, fixed_effects, suc20, glc0, ...
    cov_mat, ode_fun, n_point_simulate, rate_in_factor)
% Function that simulates the states a n_point_simulate time points for the
% nlme models. Note that this function assumes a model with 12 parameters
% on the form k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, Suc20 and Glc0. 
% Args:
%   n_cells_simulate, the number of cells to simulate 
%   n_states, the number of states in the model
%   n_parameters, the number of parameters in the model 
%   time_int, the time intervall for solving the ode system 
%   fixed_effects, the fixed population parameters 
%   suc20, the inital suc2 value 
%   glc0, the inital glucose value
%   cov_mat, the covariance matrix for the parameters
%   ode_fun, the ode function to simulate cells from
%   n_point_simulate, number of time points to simulate 
%   rate_in_factor, factor that decides the glucose change at time zero. It
%   has a default value of 1/40
% Returns:
%   simulation_result, a matrix with columns time-points, x1, x2, x3, x4
%   and index for which simulated cell it was. 

% See if the rate in factor should have defualt value 
if nargin == 10
   rate_in_factor = 1 / 40; 
end

% Set number of simulated data sets
time_stamps1 = linspace(time_int(1), time_int(2), n_point_simulate);

% Store the data in tidy format, first column is time
simulation_result = zeros(n_point_simulate ...
    * n_cells_simulate, n_states + 2);
% Add index
simulation_result(:, n_states + 2) = ...
    repelem(1:1:n_cells_simulate, n_point_simulate);
% Add time stamps
simulation_result(:, n_states + 1) = ...
    repmat(time_stamps1', n_cells_simulate, 1);

% Draw the random effects
mean_vec = zeros(1, n_parameters);
random_effects = mvnrnd(mean_vec, cov_mat, n_cells_simulate)';
% Fix the parameter values
param_values = ones(n_parameters - 2, n_cells_simulate);
suc2_init_val = zeros(1, n_cells_simulate);
glc_init_val = zeros(1, n_cells_simulate);
for i = 1:1:n_cells_simulate
    param_values(:, i) = ...
        fixed_effects .* exp(random_effects(1:end-2, i));
    suc2_init_val(i) = suc20 * exp(random_effects(end-1, i));
    glc_init_val(i) = glc0 * exp(random_effects(end, i));
end

% Surpress the warnings
warning('off','all')

% Run the simulations
for i = 1:1:n_cells_simulate
    % Allocate matrix for where to store the interpolated result
    result_ode = zeros(n_point_simulate, n_states);
    
    if mod(i, 5000) == 0
        fprintf("Iteration = %d\n", i);
    end
    
    % Solve the ODE-system
    init_val = [glc_init_val(i), 1, suc2_init_val(i), 0];
    theta_vec = [param_values(:, i)', glc_init_val(i)];
    [t, y] = ode45(@(t, y) ode_fun(t, y, theta_vec, rate_in_factor), ...
        time_int, init_val);
    
    % Interpolate the observed measurement points
    for j = 1:1:n_point_simulate
        result_ode(j, :) = interp1(t, y, time_stamps1(j));
    end
    
    % Fix the inital value
    result_ode(1, 3) = suc2_init_val(i);
    
    % Store the result
    min_ind = (i-1) * n_point_simulate + 1;
    max_ind = i * n_point_simulate;
    simulation_result(min_ind:max_ind, 1:n_states) = result_ode;
end

end