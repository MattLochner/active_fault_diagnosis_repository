#This file runs the High Pressure ODEs

using DifferentialEquations
using Plots
using full_simulation

include("BMR_gt_params.jl")
using .BMR_ground_truth_params

gtparams = data

include("BMR_initial_conditions.jl")
x0 = initial_conditions
xhat0 = initial_state_estimate  #intial_state estimate for kalman filter

estimationlength = 20      #specify the length of the state estimation period by the number of state estimates obtained at the end of the simulation.
N_meas = 1                  #h, specify the measurement sampling rate
tstart = 0
tend = tstart + estimationlength * N_meas
timesteps = tstart:N_meas:(estimationlength*N_meas)

true_vk_std = measurement_noise#specified in BMR_initial_conditions (measurement noise)

function run_ODEs(realparams, x0, true_vk_std, estimationlength, tstart, tend, timesteps, Measurement_prediction, generate_measurements, ground_truth_ODE) 

    Real_Process = zeros(length(x0),estimationlength+1)
    Real_Process[:,1] = x0                                #store initial values
    Measurements = zeros(length(true_vk_std), estimationlength+1)
    Measurements[:,1] = generate_measurements(x0, realparams, tstart, true_vk_std)              #initial measurement is set equal to the initial state estimate fed to the kalman filter

    #initialise ground truth trajectory

    GT_times = [float(tstart)]          #grows each step
    GT_states = [copy(x0)]              #vector of state vectors
    GT_meas   = [Measurement_prediction(x0, realparams, float(tstart))]


    for i in 2:(estimationlength+1)                           #the first algorithm calculation occurs after the first time step has taken place, i = 1 would call the initial values of the system
        
        #Set up equations for all fault modes
        
        gt_ODE = ODEProblem(ground_truth_ODE, Real_Process[:,i-1], (timesteps[i-1], timesteps[i]), realparams)
        
        #Only solve for fault mode that is selected to be present in the real process.
        
        #Change the ODEProblem called in this line to change which fault mode is being simulated

        ODEsol = solve(gt_ODE, Tsit5(), reltol=1e-6, abstol=1e-8, save_everystep = true, save_start = false, ) #does not save start (start is final value of previous time step)
               
        #Save simulated real process data
        append!(GT_times, ODEsol.t)
        append!(GT_states, ODEsol.u)
        append!(GT_meas, [Measurement_prediction(u, realparams, t) for (u, t) in zip(ODEsol.u, ODEsol.t)])       
        
        Real_Process[:,i] = ODEsol.u[end]

        #Generate measurements

        yk = generate_measurements(Real_Process[:,i], realparams, timesteps[i], true_vk_std)
        Measurements[:,i] = yk

        

    end

    # reshape GT_states from vector-of-vectors → matrix (n_states × n_timepoints)
    GT_states_mat = hcat(GT_states...)
    GT_meas_mat   = reduce(hcat, GT_meas)


return timesteps, Real_Process, Measurements, GT_times, GT_states_mat, GT_meas_mat
end

(timesteps, Real_Process, Measurements, GT_times, GT_states_mat, GT_meas_mat) = run_ODEs(gtparams, x0, true_vk_std, estimationlength, tstart, tend, timesteps, BMR_gt_measurement_prediction, BMR_gt_generate_measurements, calculate_gt_autoclave_odes)


plot(GT_times, GT_states_mat[38,:])

#generate "perfect measurements to compare noisy measurements against"

plot(GT_times, GT_meas_mat[13,:])

#display(scatter!(timesteps, Measurements[17,:]))

plot(GT_times, GT_meas_mat[10,:], lab = "tank 10 temp")
plot(GT_times, GT_meas_mat[11,:], lab = "tank 20 temp")
plot(GT_times, GT_meas_mat[13,:], lab = "V_9")
plot(GT_times, GT_meas_mat[4,:], lab = "temp_AC1")
