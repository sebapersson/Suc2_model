function simulated_result = calc_VC_quantiles(cov_mat, ode_fun, ...
    times_simulate, fixed_effects, suc20, glc0, sigma_val, n_individuals)
% Function that will calculate the median, 0.05, 0.2, 0.3, 0.7, 0.8
% and 0.95 quantile for a model using a covariance matrix for 
% drawing the random effects. Note that this function isn't general, it is
% written to work with 12 parameters models such as model 1 and model 2
% where the last drawn effect is the glucse value. 
% Args:
%   cov_mat, the by monolix estimated covariance matrix 
%   ode_fun, the ode_function corresponding to the model
%   times_simulate, the number of times the experiment is simulated
%   fixed_effects, a vector of the estimated rata constants
%   suc20, the estimated inital suc2 population value 
%   glc0, the estimated inital glucse value
%   sigma_val, the estimated noise for the model 
%   n_individuals, number of cells in the data set 
% Returns:
%   simulated_result, a matrix where the first column is the time indices,
%   the following columns are the 0.05, 0.2, 0.3, 0.7, 0.8, 0.95 quantiles
%   and the last column is an index corresponding to the run. 

% The mean vector for drawing the random effects
mean_vec = zeros(1, 12);

% Get time values
time_stamps1 = 0:5:480 / 480;

% Set number of simulated data sets 
n_obs = length(time_stamps1);
simulated_result = zeros(n_obs * times_simulate, 8);

% Surpress the warnings 
warning('off','all')

% Run the simulations 
for k = 1:1:times_simulate
    
    % Progress-bar
    if mod(k, 50) == 0
        fprintf("Iteration k = %d \n", k);
    end 

    % Draw the random effects 
    random_effects = mvnrnd(mean_vec, cov_mat, n_individuals)';
    % Fix the parameter values
    param_values = ones(10, n_individuals);
    suc2_init_val = ones(1, n_individuals);
    glc_init_val = ones(1, n_individuals);
    for i = 1:1:n_individuals
        % Note that last row of the random effects is Glc, the second last 
        % is Suc2 
        param_values(:, i) = fixed_effects .* exp(random_effects(1:end-2, i));
        suc2_init_val(i) = suc20 * exp(random_effects(end-1, i));
        glc_init_val(i) = glc0 * exp(random_effects(end, i));
    end 

    % Where the result will be stored 
    result_mat = zeros(n_individuals * length(time_stamps1), 3);
    result_mat(:, 1) = repmat(time_stamps1, n_individuals, 1);
    result_mat(:, 3) = repelem(1:1:n_individuals, length(time_stamps1));

    % Solve the system of ODE:s
    time_int = [0, 1];
    for i = 1:n_individuals
        % Solve the ODE-system 
        init_val = [glc_init_val(i), 1, suc2_init_val(i), 0];
        theta_vec = [param_values(:, i)', glc_init_val(i)];
        [t, y] = ode45(@(t, y) ode_fun(t, y, theta_vec), time_int, init_val);
    
        % Interpolate the result 
        y_hat = zeros(n_obs, 1);
        for j = 2:1:n_obs
            y_hat(j) = interp1(t, y(:, 3), time_stamps1(j)); 
        end
        y_hat = y_hat + normrnd(0, sigma_val, n_obs, 1);
        y_hat(1) = suc2_init_val(i);
    
        min_index = (i-1) * n_obs + 1;
        max_index = i * n_obs;
        result_mat(min_index:max_index, 2) = y_hat;
    end

    % Calculate median and the quantiles 
    [~, idx] = sort(result_mat(:, 1));
    sorted_results = result_mat(idx, :);
    result = ones(n_obs, 8);
    prob_vec = [0.05, 0.2, 0.3, 0.7, 0.8, 0.95];
    for i = 1:1:n_obs
        min_ind = (i-1) * n_individuals + 1;
        max_ind = i * n_individuals;
        value_i = sorted_results(min_ind:max_ind, :);
        median_val = median(value_i(:, 2), "omitnan");
        result(i, 1) = value_i(1, 1);
        result(i, 2) = median_val;
        result(i, 3:end) = quantile(value_i(:, 2), prob_vec);
    end
    
    % Store the result
    simulated_result((k-1)*n_obs + 1:k * n_obs, :) = result;

end 

end 