function [ind_fits, pred_fits] = compute_individual_fits(rates, ...
    start_val, time_int, n_cells, ode_fun)
% Function that for a STS model will produce the individual fits using the
% STS estimated parameters. The function will compute predicted fits and
% individuals fits, the later one is better for plotting while the first
% plot is made for an observations vs predictions plot.
% Args:
%   rates, a matrix containing the rates for each cell on the rows
%   start_val, a matrix containing the starting values on each row for each
%   cell
%   n_cells, the number of STS cells (is not allowed to exced the number of
%   rows in rates or start_val)
%   ode_fun, the ode fun used for generating the STS-values
% Returns:
%   pred_val, a n_obs * n_cells x 3 matrix containing all the predicted
%   values, where the first column is time-stamps, the second column is the
%   prediced values and the third column is the an index for each cell.
%   fitted_val, a 300 * n_cells x matrix containing the fitted values. The
%   first column is time-stamps, the second column is the fitted values and
%   the third and last column is an index per cell.


% Read the time-stamps
time_stamps1 = ((0:5:480) / 480)';

% Number of observations in the data set
n_obs = length(time_stamps1);

% Where to store the result
time_stamp_fit = linspace(time_int(1), time_int(2), 300)';
obs_per_fit = length(time_stamp_fit);
% The fitted responses
ind_fits = zeros(obs_per_fit * n_cells, 3);
ind_fits(:, 1) = repmat(time_stamp_fit, n_cells, 1);
ind_fits(:, 3) = repelem(1:1:n_cells, obs_per_fit);
% The predicted responses
pred_fits = zeros(n_obs * n_cells, 3);
pred_fits(:, 3) = repelem(1:1:n_cells, n_obs);
pred_fits(:, 1) = repmat(time_stamps1, n_cells, 1);

% Perform the computations
for i = 1:1:n_cells
    rates_vec = rates(i, :);
    start_vec = start_val(i, :);
    [t, y] = ode45(@(t, y) ode_fun(t, y, rates_vec), ...
        time_int, start_vec);
    
    % Interpolate the fitted values
    for j = 1:1:obs_per_fit
        index = (i - 1) * obs_per_fit + j;
        ind_fits(index, 2) = interp1(t, y(:, 3), time_stamp_fit(j));
    end
    
    % Interpolate the predicted values
    for j = 1:1:n_obs
        index = (i - 1) * n_obs + j;
        pred_fits(index, 2) = interp1(t, y(:, 3), time_stamps1(j));
    end
end

end