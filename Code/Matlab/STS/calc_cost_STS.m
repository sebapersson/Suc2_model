function cost = calc_cost_STS(time_stamps, data_val, ... 
    time_int, ode_fun, theta)
% This function calculates the cost fuction value for the lsqnonlin 
% function for model1 and 2 with all parameters for the STS approach.
% An assumption for this file is that the model contains 12 parameters, 
% in the order k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, Suc20 and Glc0. 
% Args:
%   time_stamps, the observed time stamps for the cells 
%   data_val, the correpsonding data values for the observed time-stamps
%   time_int, the time intervall for solving the ODE-system 
%   ode_fun, the ode-model required to compute the cost function
%   theta, the parameter space for the optimisation 
% Returns:
%   cost, the cost function values requried by lsqnonlin 

% Fix the starting values
inital_glc = theta(end);        % Initial glucose 
initial_Suc2 = theta(end-1);    % Inital Suc2
theta_ode = theta(1:end-2);     % Rate parameters 
initial_val_ode = [inital_glc, 1, initial_Suc2, 0];

% Theta vec also contains glucose for model 7 
theta_vec = [theta_ode, inital_glc];

% Do not display numerical warnings 
warning( 'off' , 'all' )

% Solve the ODE system
[t, y] = ode45(@(t, y) ode_fun(t, y, theta_vec), time_int, initial_val_ode);

% Interpolate to the measured values 
n_data_points = length(data_val);   
y_hat = zeros(n_data_points, 1);
for i = 1:1:n_data_points
   y_hat(i) = interp1(t, y(:, 3), time_stamps(i)); 
end

% Calculate the elements of the cost function 
cost = y_hat - data_val;

