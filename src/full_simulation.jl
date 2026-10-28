module full_simulation

using LinearAlgebra
using DifferentialEquations

#function related to EKF

include("extended_kalman_filter.jl")
export lyap!
export check_covariance
export ekf
export run_ODE_model #check that this is the same as the one exported from kalman_filter

#functions related to BMR example

include("BMR_gt_generate_measurements.jl")
export BMR_gt_measurement_prediction
export BMR_gt_measurement_prediction_with_samples
export BMR_gt_generate_measurements
export BMR_gt_generate_measurements_with_samples

include("BMR_kf_generate_measurements.jl")
export BMR_kf_measurement_prediction
export BMR_kf_measurement_prediction_with_samples

include("BMR_gt_ODEs.jl")
export calculate_gt_autoclave_odes
export calculate_oxygen_solubility2
export calculate_reaction_rates
export liq_fractions
export solid_fractions
export calculate_model_autoclave_odes

#must still add function to determine solid fraction...

#Functions related to three tank problem

include("3_tank_problem_ODEs.jl")
export Params
export Three_tank_GT_process!
export Three_tank_model_process!

include("3_tank_measurement_prediction.jl")
export measurement_prediction_3_tank_system

include("3_tank_generate_measurements.jl")
export generate_measurements_3_tank_system

include("multiple_model_method.jl")
export multiple_model_method
export flatten_prior
export linearflattening

#Okay let's get this working!

end # module full_simulation
