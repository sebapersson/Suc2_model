function simulation_result = simulate_single_cells_STS(...
    n_cells_simulate, n_states, mean_vec, cov_mat, time_int, ode_fun)
% Function that will simulate n_cells_simulate cells when using the STS
% approach. This means that the random parameters are drawn from a normal
% distribution, and following this the exp-function is applied. The
% covaraince matrix is estimated from all the cells. 
% Args:
%   n_cells_simulate, the number of cells to simulate 
%   n_states, the number of states in the model 
%   mean_vec, the mean vector computed from the STS estimates 
%   cov_mat, the covariance matrix computed from the STS estimates 
%   time_int, the int in where to solve the ODE 
%   ode_fun, the ode fun for the current model (4 or 6)
% Returns:
%   simulation_result, a matrix where the first n_states column are
%   simulated values for the states, the second last is time index and the
%   least is simulated cell id 

% Read the time-stamps
time_stamps1 = ((0:5:480) / 480)';
% Set number of simulated data sets
n_obs = length(time_stamps1);

% Store the data in tidy format, first column is time
simulation_result = zeros(n_obs * n_cells_simulate, n_states + 2);
% Add index
simulation_result(:, n_states + 2) = ...
    repelem(1:1:n_cells_simulate, n_obs);
% Add time stamps
simulation_result(:, n_states + 1) = ...
    repmat(time_stamps1, n_cells_simulate, 1);

% Draw the random effects
% Note second last column is Suc2, last column is Glc 
rand_effects = mvnrnd(mean_vec, cov_mat, n_cells_simulate)';
rates = exp(rand_effects(1:end-2, :));
rates = [rates; exp(rand_effects(end, :))];
suc2_init_val = exp(rand_effects(end-1, :));
glc_init_val = exp(rand_effects(end, :));

% Surpress the warnings 
warning('off','all')

% Run the simulations
for i = 1:1:n_cells_simulate
    % Allocate matrix for where to store the interpolated result
    result_ode = zeros(n_obs, n_states);
    
    if mod(i, 5000) == 0
       fprintf("Iteration = %d\n", i); 
    end
    
    % Solve the ODE-system
    init_val = [glc_init_val(i), 1, suc2_init_val(i), 0];
    rate_vec = rates(:, i);
    [t, y] = ode45(@(t, y) ode_fun(t, y, rate_vec), time_int, init_val);
    
    % Interpolate the observed measurement points
    for j = 1:1:n_obs
        result_ode(j, :) = interp1(t, y, time_stamps1(j));
    end
    
    % Fix the inital value
    result_ode(1, 3) = suc2_init_val(i);
    
    % Store the result
    min_ind = (i-1) * n_obs + 1;
    max_ind = i * n_obs;
    simulation_result(min_ind:max_ind, 1:n_states) = result_ode;
end

end 