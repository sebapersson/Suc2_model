# Function that describes the ODE-system for the model 2 with none short
# delay.
# Args:
#   t, the time
#   u, the different states with:
#       u[1], intralcellular glucose
#       u[2], Snf1p
#       u[3], Suc2 (the outsignal)
#       u[4], X the cellular response to starvation
#   p, paramater vector where p[1] -- p[10] corresponds to the
#   rates k1, ..., k10. p[11] is inital glucose
# Returns:
#   The deriviatives for the system at time t which is used by the
#   ODE-solvers in julia
function model1(du, u, p, t)
    # Defining the rates from parameter vector
    k1 = p[1]
    k2 = p[2]
    k3 = p[3]
    k4 = p[4]
    k5 = p[5]
    k6 = p[6]
    k7 = p[7]
    k8 = p[8]
    k9 = p[9]
    k10 = p[10]
    Glc0 = p[11]

    # Fix the rate in
    if t < 0.01
        rate_in = k1 * 1;
    else
        rate_in = k1 * 1 / 40;
    end

    # The model
    du[1] = rate_in - k2 * u[1] + k3 * u[4]
    du[2] = k4 * (Glc0 - u[1]) - k5 * u[2]
    du[3] = k6 * u[2]^2 - k7 * u[3]
    du[4] = (k10 * u[2] + k7 * u[3]^2) / (u[1] * u[2] + k8) - k9 * u[4]
end


# Function that describes the ODE-system for the model 2 with none short
# delay.
# Args:
#   t, the time
#   u, the different states with:
#       u[1], intralcellular glucose
#       u[2], Snf1p
#       u[3], Suc2 (the outsignal)
#       u[4], X the cellular response to starvation
#   p, paramater vector where p[1] -- p[10] corresponds to the
#   rates k1, ..., k10.
# Returns:
#   The deriviatives for the system at time t which is used by the
#   ODE-solvers in julia
function model2(du, u, p, t)
    # Defining the rates from parameter vector
    k1 = p[1]
    k2 = p[2]
    k3 = p[3]
    k4 = p[4]
    k5 = p[5]
    k6 = p[6]
    k7 = p[7]
    k8 = p[8]
    k9 = p[9]
    k10 = p[10]

    # Fix the rate in
    if t < 0.01
        rate_in = k1 * 1;
    else
        rate_in = k1 * 1 / 40;
    end

    # The model
    du[1] = rate_in - k2 * u[1] + k3 * u[4]
    du[2] = k4 / (k8 + u[1]^2) - k5 * u[2]
    du[3] = k6 * u[2]^2 - k7 * u[3]^2
    du[4] = (k10 * u[2] + k7 * u[3]^2) / (u[1] * u[3] + k8) - k9 * u[4]
end


# Function that describes the ODE-system for the model 2 with short
# delay.
# Args:
#   t, the time
#   u, the different states with:
#       u[1], intralcellular glucose
#       u[2], Snf1p
#       u[3], Suc2 (the outsignal)
#       u[4], X the cellular response to starvation
#   p, paramater vector where p[1] -- p[10] corresponds to the
#   rates k1, ..., k10.
# Returns:
#   The deriviatives for the system at time t which is used by the
#   ODE-solvers in julia
function model2_short_del(du, u, p, t)
    # Defining the rates from parameter vector
    k1 = p[1]
    k2 = p[2]
    k3 = p[3]
    k4 = p[4]
    k5 = p[5]
    k6 = p[6]
    k7 = p[7]
    k8 = p[8]
    k9 = p[9]
    k10 = p[10]

    # Fix the rate in
    if t < 0.0001
        rate_in = k1 * 1;
    else
        rate_in = k1 * 1 / 40;
    end

    # The model
    du[1] = rate_in - k2 * u[1] + k3 * u[4]
    du[2] = k4 / (k8 + u[1]^2) - k5 * u[2]
    du[3] = k6 * u[2]^2 - k7 * u[3]^2
    du[4] = (k10 * u[2] + k7 * u[3]^2) / (u[1] * u[3] + k8) - k9 * u[4]
end
