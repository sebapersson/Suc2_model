function [parameters, residuals] = perform_STS(cell_data, theta_start,...
    n_parameters, time_int, n_cells,  obj_fun, ode_fun)
% Function that will perform the STS for n_cells and returns the 
% estimated parameters and residuals for each fit. 
% Args:
%   cell_data, a matrix containing the observed data, with first column
%   being time, second column being observed values and third column being
%   a cell index. The data is assumed to be tidy 
%   theta_start, the starting parameters for each cell
%   n_parameters, the number of parameters to be estiamted in the model
%   time_int, the time interval where the ODE-system is solved
%   n_cells, the number of cells to perform (has to be smaller or equal to
%   the maximum index)
%   ode_fun, the function handle to the ODE-model
%   obj_fun, a function handle to the objective function (file thata
%   calculates the cost)
% Returns:
%   param, a matrix by n_cells x n_parameters + 1, where the last column is
%   the cost function. 
%   res, a matrix of n_cells x n_obs that contains the residuals for each
%   cell 

% Matrix where the result is stored 
parameters = zeros(n_cells, n_parameters + 1);
n_obs = length(find(cell_data(:, 3) == 1));
residuals = zeros(n_cells, n_obs);

% Compuate the time-stamps for scaled time 
time_stamps1 = (0:5:480) / 480;

% Bounds for the optimisation problem 
lb = zeros(n_parameters, 1);
ub = ones(n_parameters, 1) * 2000;

% Turning of notifications 
opts1=  optimset('display','off');

for i = 1:1:n_cells
    % Print progress per tenth cell 
    if mod(i, 10) == 0
        fprintf("i = %d \n", i); 
    end
    
    % Filter out the specific cell 
    i_cell_i = find(cell_data(:, 3) == i);
    data_values = cell_data(i_cell_i, 2);
    
    % Define the cost function 
    my_objective = @(theta) obj_fun(time_stamps1, ... 
        data_values, time_int, ode_fun, theta);
    
    % Solve the optimisation problem 
    [min_res, cost_value, residuals_ode] = ...
        lsqnonlin(my_objective, theta_start, lb, ub, opts1);
    
    % Store result in the result matrix
    parameters(i, 1:end - 1) = min_res;
    parameters(i, end) = cost_value;
    
    % Store the residuals 
    residuals(i, :) = residuals_ode;
end 

end 