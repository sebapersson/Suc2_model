using CSV
using DifferentialEquations
using Distributions
using Printf
using DataFrames

# Will have to change this later
cd("/home/sebpe/Dropbox/Chalmers/Suc2_model/Code/Julia/")

# ------------------------------------------------------------------------
# Start of functions
# ------------------------------------------------------------------------
# Load the models
include("Models/Models.jl")

# Function that will generate the inital and paramter values required for simulating cells.
# Args
#    n_cells_simulate, the number of cells to simulate
#    n_rates, the number of rates in the model
#    n_init_val, the number of inital values in the model
#    fixed_effects, an array with the fixed effects, the first values are rates
#        k1, ..., kn
#    init_one_list, a tuple with the initial values that should be 1
#    init_rand, a 2-dimensional tuple signalling which random parameters are intial values
#        where for each entry the first element is the position in the inital vector, while
#        the second element is the position in the fixed effects vector.
# Returns
#    rates, the rates to use for simulations
#    inital_values, the inital values to use for simulations, note that
#        individuals are given by columns
function generate_parameters(n_cells_simulate, n_rates, n_init_val, fixed_effects,
    init_one_list, init_rand, cov_mat)

    # Generate the inital value and parameter matrix
    dist_mult_normal = MvNormal(cov_mat)
    rand_effects = rand(dist_mult_normal, n_cells_simulate)
    # For sake of speed apply exponential function once
    rand_effects = exp.(rand_effects)

    # Fixing the parameters
    param_matrix = zeros(n_rates, n_cells_simulate)
    for i in 1:n_rates
        param_matrix[i, :] = fixed_effects[i] .* rand_effects[i, :]
    end

    # Fixing the inital values
    init_matrix = zeros(n_init_val, n_cells_simulate)
    for i in init_one_list
        init_matrix[i, :] = fill(1, n_cells_simulate)
    end
    for i in 1:length(init_rand)
        i_mat, i_fixed_effects = init_rand[i][1], init_rand[i][2]
        init_matrix[i_mat, :] = fixed_effects[i_fixed_effects] .* (rand_effects[i_fixed_effects, :])
    end

    return param_matrix, init_matrix
end


# Function that will simulate n_cells_simulate number of cells using a
# user provided parameter matrix and inital value matrix, plus
# a user provided ODE-function.
# Args:
#    n_cells_simulate, the number of cells to simulate
#    n_states, the number of states in the model
#    time_span_inter, the time span which to interpolate over
#    time_span_solve, the time span to solve ODE-system in
#    param_val, the parameter matrix, where each column corresponds to the
#        parameter values for a certain cell
#    init_val, the inital value matrix, where each column corresponds to the
#        inital values for a cell.
#    ode_model, the ode model, this is a function!
# Returns:
#    simulated_result, a matrix with column 1 being time, column 2 being
#        time  index for a cell and column 3:n_states + 2 being the states
function simulate_cells(n_cells_simulate, n_states, time_span_inter,
    time_span_solve, param_val, init_val, ode_model::Function)

    n_data_points_per_cell = length(time_span_inter)
    # Allocate array to hold the final result
    simulated_result = zeros(  n_data_points_per_cell * n_cells_simulate,
    n_states + 2)
    # First column time-stamps
    simulated_result[:, 1] = repeat(time_span_inter, outer = n_cells_simulate)
    # Second column an index
    simulated_result[:, 2] = repeat(1:n_cells_simulate, inner = n_data_points_per_cell)
    # Column 3-n_states + 2 equals the states

    # Simulate the cells
    counter_fail = 0
    success_symbol = Symbol("Success")
    for i in 1:n_cells_simulate

        if i % 5000 == 0
            @printf("iteration = %d\n", i)
        end

        # Define the ODE-problems
        ode_problem = ODEProblem(ode_model, init_val[:, i], time_span_solve, param_val[:, i])

        # Solve the ODE-problem, it's non stiff
        ode_solution = solve(ode_problem, Tsit5(), verbose = false)

        # Handling integration problems
        if ode_solution.retcode != success_symbol
            counter_fail += 1
            continue
        end

        # Interpolate the solution into the wanted time-slots
        i_min = (i - 1) * n_data_points_per_cell
        for j in 1:n_data_points_per_cell
            simulated_result[i_min + j, 3:n_states+2] = ode_solution(time_span_inter[j])
        end
    end

    @printf("Counter fail = %d\n", counter_fail)
    return simulated_result
end

# ------------------------------------------------------------------------
# End of functions
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Simulate the NLME-models
# ------------------------------------------------------------------------
function model1_nlme()
    # Read csv file
    cov_mat_raw = CSV.read("../../Intermediate/Cov_mat_model1_nlme.csv")
    cov_mat = convert(Matrix, cov_mat_raw[1:end, 2:end])

    # Read the fixed effects
    fixed_effects_raw = CSV.read("../Monolix/Model1/Result/populationParameters.txt")
    fixed_effects = convert(Array, fixed_effects_raw[1:12, 2])

    # Provide the time-span for interpolation
    time_span_inter = range(0, 1, length = 200)

    # Generate the inital values
    n_cells_simulate = 50000
    n_rates = 10; n_init_val = 4
    init_one_list = (2); init_rand = ((1, 12), (3, 11))
    param_val, init_val = generate_parameters(n_cells_simulate, n_rates,
    n_init_val, fixed_effects, init_one_list, init_rand, cov_mat)
    # Add the inital Glc-level as parameter
    param_val = [param_val; init_val[1, :]']

    # Simulate the cells
    n_data_points_per_cell = length(time_span_inter)
    n_states = 4
    time_span_solve = (0.0, 1.0)
    simulated_result = simulate_cells(n_cells_simulate, n_states, time_span_inter,
        time_span_solve, param_val, init_val, model1)

    # Write result to file
    simulated_result_table = convert(DataFrame, simulated_result)
    col_names = ["time", "index", "Glc", "SNF1", "SUC2", "X"]
    names!(simulated_result_table, Symbol.(col_names))
    CSV.write("../../Intermediate/Simulated_model1_nlme_julia.csv",
    simulated_result_table)

    return 0
end

function model2_nlme()
    # Read csv file
    cov_mat_raw = CSV.read("../../Intermediate/Cov_mat_model2_nlme.csv")
    cov_mat = convert(Matrix, cov_mat_raw[1:end, 2:end])

    # Read the fixed effects
    fixed_effects_raw = CSV.read("../Monolix/Model2/Result/populationParameters.txt")
    fixed_effects = convert(Array, fixed_effects_raw[1:12, 2])

    # Provide the time-span for interpolation
    time_span_inter = range(0, 1, length = 200)

    # Generate the inital values
    n_cells_simulate = 50000
    n_rates = 10; n_init_val = 4
    init_one_list = (2); init_rand = ((1, 12), (3, 11))
    param_val, init_val = generate_parameters(n_cells_simulate, n_rates,
      n_init_val, fixed_effects, init_one_list, init_rand, cov_mat)

    # Simulate the cells
    n_data_points_per_cell = length(time_span_inter)
    n_states = 4
    time_span_solve = (0.0, 1.0)
    simulated_result = simulate_cells(n_cells_simulate, n_states, time_span_inter,
        time_span_solve, param_val, init_val, model2)

    # Write result to file
    simulated_result_table = convert(DataFrame, simulated_result)
    col_names = ["time", "index", "Glc", "SNF1", "SUC2", "X"]
    names!(simulated_result_table, Symbol.(col_names))
    CSV.write("../../Intermediate/Simulated_model2_nlme_julia.csv",
    simulated_result_table)

    return 0
end


function model2_short_del_nlme()
    # Read csv file
    cov_mat_raw = CSV.read("../../Intermediate/Cov_mat_model2_short_del_nlme.csv")
    cov_mat = convert(Matrix, cov_mat_raw[1:end, 2:end])

    # Read the fixed effects
    fixed_effects_raw = CSV.read("../Monolix/Model2_test/Result/populationParameters.txt")
    fixed_effects = convert(Array, fixed_effects_raw[1:12, 2])

    # Provide the time-span for interpolation
    time_span_inter = range(0, 1, length = 200)

    # Generate the inital values
    n_cells_simulate = 50000
    n_rates = 10; n_init_val = 4
    init_one_list = (2); init_rand = ((1, 12), (3, 11))
    param_val, init_val = generate_parameters(n_cells_simulate, n_rates,
      n_init_val, fixed_effects, init_one_list, init_rand, cov_mat)

    # Simulate the cells
    n_data_points_per_cell = length(time_span_inter)
    n_states = 4
    time_span_solve = (0.0, 1.0)
    simulated_result = simulate_cells(n_cells_simulate, n_states, time_span_inter,
        time_span_solve, param_val, init_val, model2_short_del)

    # Write result to file
    simulated_result_table = convert(DataFrame, simulated_result)
    col_names = ["time", "index", "Glc", "SNF1", "SUC2", "X"]
    names!(simulated_result_table, Symbol.(col_names))
    CSV.write("../../Intermediate/Simulated_model2_nlme_short_del_julia.csv",
    simulated_result_table)

    return 0
end


function model1_sts()
    # Reading the covariance matrix and fixed effects, last row mean vector
    parameter_data = CSV.read("../../Result/Files/STS_model1_cov_mean.csv")
    cov_mat = convert(Matrix, parameter_data[1:end-1, :])
    fixed_effects = exp.(convert(Vector, parameter_data[end, :]))

    # Generate the inital values
    n_cells_simulate = 50000
    n_rates = 10; n_init_val = 4
    init_one_list = (2); init_rand = ((1, 12), (3, 11))
    param_val, init_val = generate_parameters(n_cells_simulate, n_rates,
      n_init_val, fixed_effects, init_one_list, init_rand, cov_mat)
    param_val = [param_val; init_val[1, :]']

    # Provide the time-span for interpolation
    time_span_inter = range(0, 1, length = 200)

    n_data_points_per_cell = length(time_span_inter)
    n_states = 4
    time_span_solve = (0.0, 1.0)
    simulated_result = simulate_cells(n_cells_simulate, n_states, time_span_inter,
      time_span_solve, param_val, init_val, model1)

    # Write result to file
    simulated_result_table = convert(DataFrame, simulated_result)
    col_names = ["time", "index", "Glc", "SNF1", "SUC2", "X"]
    names!(simulated_result_table, Symbol.(col_names))
    CSV.write("../../Intermediate/Simulated_model1_STS_julia.csv",
      simulated_result_table)

    return 0
end


function model2_sts()
    # Reading the covariance matrix and fixed effects, last row mean vector
    parameter_data = CSV.read("../../Result/Files/STS_model2_cov_mean.csv")
    cov_mat = convert(Matrix, parameter_data[1:end-1, :])
    fixed_effects = exp.(convert(Vector, parameter_data[end, :]))

    # Generate the inital values
    n_cells_simulate = 50000
    n_rates = 10; n_init_val = 4
    init_one_list = (2); init_rand = ((1, 12), (3, 11))
    param_val, init_val = generate_parameters(n_cells_simulate, n_rates,
      n_init_val, fixed_effects, init_one_list, init_rand, cov_mat)

    # Provide the time-span for interpolation
    time_span_inter = range(0, 1, length = 200)

    n_data_points_per_cell = length(time_span_inter)
    n_states = 4
    time_span_solve = (0.0, 1.0)
    simulated_result = simulate_cells(n_cells_simulate, n_states, time_span_inter,
      time_span_solve, param_val, init_val, model2)

    # Write result to file
    simulated_result_table = convert(DataFrame, simulated_result)
    col_names = ["time", "index", "Glc", "SNF1", "SUC2", "X"]
    names!(simulated_result_table, Symbol.(col_names))
    CSV.write("../../Intermediate/Simulated_model2_STS_julia.csv",
      simulated_result_table)

    return 0
end


# Running the code
model1_nlme()
model2_nlme()
model2_short_del_nlme()
model1_sts()
model2_sts()
