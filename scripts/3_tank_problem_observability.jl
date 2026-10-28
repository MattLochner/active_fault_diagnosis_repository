using full_simulation
using FiniteDiff
using LinearAlgebra

modelparams = Params(6e-5, 6e-5, 1, 5e-5, 9.81, 0.0154, 0.8, 1.0, 0.6, 0.002, 0)

# Non-degenerate operating point (equal levels cause undefined Jacobian due to sqrt(|xi - xj|))
x0 = [0.6, 0.4, 0.5]

# Wrap RP_F1! as a plain x -> dx function for ForwardDiff
f(x) = begin
    dx = Three_tank_model_process!(x, modelparams, 0.0)
    return dx
end

# A matrix: Jacobian of f at operating point
A = FiniteDiff.finite_difference_jacobian(f, x0)

h(x) = begin
    dx = measurement_prediction_3_tank_system(x, modelparams, 0.0)
    return dx
end

# C matrix: Jacobian of measurement function at operating point
C = FiniteDiff.finite_difference_jacobian(h, x0)

# Observability matrix O = [C; CA; CA²]
O = [C; C*A; C*A^2]

obs_rank = rank(O)
println("A matrix:\n", A)
println("\nC matrix:\n", C)
println("\nObservability matrix O:\n", O)
println("\nRank of O: $obs_rank / 3")
println("System is $(obs_rank == 3 ? "observable" : "NOT observable") at x0 = $x0")
