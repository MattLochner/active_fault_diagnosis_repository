using DifferentialEquations
using Plots
using LinearAlgebra
using full_simulation

include("BMR_gt_params.jl")
using .BMR_ground_truth_params

include("BMR_kf_params.jl")
using .Kf_model_params

gtparams = data
kfparams = Kf_data

gtparams.faultmode = 0

include("BMR_initial_conditions.jl")
x0 = initial_conditions
xhat0 = initial_state_estimate.*1.05  #initial_state estimate for kalman filter

estimationlength = 200      #specify the length of the state estimation period by the number of state estimates obtained at the end of the simulation.
N_meas = 0.1               #h, specify the measurement sampling rate
#NB: Samples are taken every 1 hour. N_meas must be factor of 1 hour so that the samples and N_meas are taken simultaneously at a specific point.
tstart = 0
tend = tstart + estimationlength * N_meas
timesteps = tstart:N_meas:(estimationlength*N_meas)

#true_vk_std #specified in BMR_initial_conditions (measurement noise)
#true_vk_std = measurement_noise
true_vk_std = meas_noise_3 #much tighter measurements

vk_std = true_vk_std

nx = length(xhat0)
ny = 36
P0 = Matrix(diagm((initial_state_estimate - xhat0).^2))
Q = Matrix(diagm(xhat0.^2).*1e-4)    #diagonal matrix scaled relative to variable size
R = Matrix(diagm(vk_std.^2))

#Bayesian update of fault likelihoods
Priors = [0.333, 0.333, 0.334]

#Run the kalman filter

function run_ODEs(realparams, x0, true_vk_std, Q, R, P0, nx, ny, estimationlength, tstart, tend, timesteps, generate_measurements, GT_measurement_prediction, ground_truth_ODE, xhat0, kalman_ODE, kalmanparams, Priors) 

    Real_Process = zeros(length(x0),estimationlength+1)
    Real_Process[:,1] = x0                                #store initial values
    Measurements = zeros(length(true_vk_std), estimationlength+1)
    Measurements[:,1] = generate_measurements(x0, realparams, tstart, true_vk_std) #initial measurement should maybe be set to xhat0 not x0? 

    kalmanparams.faultmode = 0

    Normal_StateEstimates = zeros(length(xhat0), estimationlength+1)
    Normal_StateEstimates[:,1] = xhat0
    Normal_CovarianceEstimates = fill(zeros(nx,nx), estimationlength+1)
    Normal_CovarianceEstimates[1] = P0
    Normal_Measurement_Estimates = zeros(ny, estimationlength+1)
    Normal_Measurement_Estimates[:,1] = BMR_kf_measurement_prediction(xhat0, kalmanparams, tstart)
    Normal_Measurement_Covariance = fill(zeros(ny, ny), estimationlength+1)

    kalmanparams.faultmode = 4

    F1_StateEstimates = zeros(nx, estimationlength+1)
    F1_StateEstimates[:,1] = xhat0
    F1_CovarianceEstimates = fill(zeros(nx,nx), estimationlength+1)
    F1_CovarianceEstimates[1] = P0
    F1_Measurement_Estimates = zeros(ny, estimationlength+1)
    F1_Measurement_Estimates[:,1] = BMR_kf_measurement_prediction(xhat0, kalmanparams, tstart)
    F1_Measurement_Covariance = fill(zeros(ny, ny), estimationlength+1)

    kalmanparams.faultmode = 5

    F2_StateEstimates = zeros(nx, estimationlength+1)
    F2_StateEstimates[:,1] = xhat0
    F2_CovarianceEstimates = fill(zeros(nx,nx), estimationlength+1)
    F2_CovarianceEstimates[1] = P0
    F2_Measurement_Estimates = zeros(ny, estimationlength+1)
    F2_Measurement_Estimates[:,1] = BMR_kf_measurement_prediction(xhat0, kalmanparams, tstart)
    F2_Measurement_Covariance = fill(zeros(ny, ny), estimationlength+1)   
        
    PosteriorBeliefs = zeros(3, estimationlength+1) #first entry is number of fault modes
    PosteriorBeliefs[:,1] = Priors 
        
    #Set up initial variables that will be modified in the loop:

    N_xhatmin1 = xhat0
    F1_xhatmin1 = xhat0
    F2_xhatmin1 = xhat0
    N_Pposmin1 = P0
    F1_Pposmin1 = P0
    F2_Pposmin1 = P0
    Priors = Priors

    #initialise ground truth trajectory

    GT_times = [float(tstart)]          #grows each step
    GT_states = [copy(x0)]              #vector of state vectors
    GT_meas   = [GT_measurement_prediction(x0, realparams, float(tstart))]

    for i in 2:(estimationlength+1)                           #the first algorithm calculation occurs after the first time step has taken place, i = 1 would call the initial values of the system
        
        #Set up equations for all fault modes
        println("i = $i")
        
        gt_ODE = ODEProblem(ground_truth_ODE, Real_Process[:,i-1], (timesteps[i-1], timesteps[i]), realparams)
        
        #Only solve for fault mode that is selected to be present in the real process.
        
        #Change the ODEProblem called in this line to change which fault mode is being simulated

        ODEsol = solve(gt_ODE, Tsit5(), reltol=1e-6, abstol=1e-8, save_everystep = true, save_start = false, ) #does not save start (start is final value of previous time step)
               
        #Save simulated real process data
        append!(GT_times, ODEsol.t)
        append!(GT_states, ODEsol.u)
        append!(GT_meas, [GT_measurement_prediction(u, realparams, t) for (u, t) in zip(ODEsol.u, ODEsol.t)])       
        
        Real_Process[:,i] = ODEsol.u[end]

        #Generate measurements
        yk = generate_measurements(Real_Process[:,i], realparams, timesteps[i], true_vk_std)
        Measurements[:,i] = yk
        
        #PROPOGATION OF UNCERTAINTY + STATE ESTIMATION WITH KALMAN FILTER
        #****************************************************************
        tStart = timesteps[i-1]
        tEnd = timesteps[i]

        #Normal Mode
            #set fault mode = 0
            kalmanparams.faultmode = 0
        
        (N_yhat, N_Py, N_xhatpos, N_Ppos) = ekf(N_xhatmin1, N_Pposmin1, nx, ny, R, Q, kalmanparams, tStart, tEnd, BMR_kf_measurement_prediction, kalman_ODE, yk)

        Normal_StateEstimates[:,i] = N_xhatpos
        Normal_CovarianceEstimates[i] = N_Ppos
        Normal_Measurement_Estimates[:,i] = N_yhat
        Normal_Measurement_Covariance[i] = N_Py

        N_xhatmin1 = N_xhatpos
        N_Pposmin1 = N_Ppos

        #Faultmode 1
            #set fault mode = 4
            kalmanparams.faultmode = 4

            (F1_yhat, F1_Py, F1_xhatpos, F1_Ppos) = ekf(F1_xhatmin1, F1_Pposmin1, nx, ny, R, Q, kalmanparams, tStart, tEnd, BMR_kf_measurement_prediction, kalman_ODE, yk)

        F1_StateEstimates[:,i] = F1_xhatpos
        F1_CovarianceEstimates[i] = F1_Ppos
        F1_Measurement_Estimates[:,i] = F1_yhat
        F1_Measurement_Covariance[i] = F1_Py

        F1_xhatmin1 = F1_xhatpos
        F1_Pposmin1 = F1_Ppos 

        #Faultmode 2
            #set fault mode = 5
            kalmanparams.faultmode = 5

            (F2_yhat, F2_Py, F2_xhatpos, F2_Ppos) = ekf(F2_xhatmin1, F2_Pposmin1, nx, ny, R, Q, kalmanparams, tStart, tEnd, BMR_kf_measurement_prediction, kalman_ODE, yk)
  
        F2_StateEstimates[:,i] = F2_xhatpos
        F2_CovarianceEstimates[i] = F2_Ppos
        F2_Measurement_Estimates[:,i] = F2_yhat
        F2_Measurement_Covariance[i] = F2_Py

        F2_xhatmin1 = F2_xhatpos
        F2_Pposmin1 = F2_Ppos 

        #BAYESIAN UPDATE OF FAULT LIKELIHOODS

        yhats = [N_yhat'; F1_yhat'; F2_yhat']
        covariances = [N_Py, F1_Py, F2_Py]

        PosteriorBeliefs[:,i] = multiple_model_method(yk, PosteriorBeliefs[:,(i-1)], yhats, covariances)

    end

    # reshape GT_states from vector-of-vectors → matrix (n_states × n_timepoints)
    GT_states_mat = hcat(GT_states...)
    GT_meas_mat   = reduce(hcat, GT_meas)


return timesteps, GT_times, GT_states_mat, GT_meas_mat, Normal_StateEstimates, Normal_CovarianceEstimates, F1_StateEstimates, F1_CovarianceEstimates, F2_StateEstimates, F2_CovarianceEstimates, Real_Process, Measurements, PosteriorBeliefs, Normal_Measurement_Estimates, F1_Measurement_Estimates, F2_Measurement_Estimates
end

(timesteps, GT_times, GT_states_mat, GT_meas_mat, Normal_StateEstimates, Normal_CovarianceEstimates, F1_StateEstimates, F1_CovarianceEstimates, F2_StateEstimates, F2_CovarianceEstimates, Real_Process, Measurements, PosteriorBeliefs, Normal_Measurement_Estimates, F1_Measurement_Estimates, F2_Measurement_Estimates) = run_ODEs(gtparams, x0, true_vk_std, Q, R, P0, nx, ny, estimationlength, tstart, tend, timesteps, BMR_gt_generate_measurements, BMR_gt_measurement_prediction, calculate_gt_autoclave_odes, xhat0, calculate_model_autoclave_odes, kfparams, Priors)

function CalcStdDevs(x, P)
    x_std = fill(0.0, size(x))

    for i in 1:size(x)[1]
        global x_std[i,:] = sqrt.([P[j][i,i] for j in 1:length(P)])
    end
    return x_std
end

#I'm only going to do this for specific states i want to plot, else it's too much...

Normal_x_std = CalcStdDevs(Normal_StateEstimates, Normal_CovarianceEstimates)

#Plot all measurement predictions against measurements and true state to understand Bayesian update of fault likelihood:


# for j in 1:36
# plot(GT_times, GT_meas_mat[j, :], lab = "GT", title = "Measurement $j")
# plot!(timesteps, Normal_Measurement_Estimates[j, :], lab = "Normal")
# plot!(timesteps, F1_Measurement_Estimates[j, :], lab = "F1")
# plot!(timesteps, F2_Measurement_Estimates[j, :], lab = "F2")
# display(scatter!(timesteps, Measurements[j,:]))
# end

#generate "perfect measurements to compare noisy measurements against"
plot(GT_times, GT_states_mat[116,:])
plot!(timesteps, Normal_StateEstimates[116,:], ribbon = (Normal_x_std[116, :], Normal_x_std[116, :]))

scatter(timesteps, Measurements[4, :], lab = "Measurement T_4")
plot!(timesteps, Normal_StateEstimates[100, :], ribbon = (Normal_x_std[112, :], Normal_x_std[112, :]), lab = "T_4 state estimate")
plot!(GT_times, GT_meas_mat[4, :], lab = "Tank 1", xlab = "time (h)", ylab = "temperature (deg C)")

plot(GT_times, GT_meas_mat[4, :], lab = "Tank 10", xlab = "time (h)", ylab = "temperature (deg C)")

plot(GT_times, GT_meas_mat[10, :], lab = "Tank 20 ground truth", xlab = "time (h)", ylab = "temperature (deg C)")
plot!(timesteps, Normal_StateEstimates[21, :], lab = "Normal estimate", xlab = "time (h)", ylab = "temperature (deg C)")
plot!(timesteps, F1_StateEstimates[21, :], lab = "F1 estimate", xlab = "time (h)", ylab = "temperature (deg C)")
plot!(timesteps, F2_StateEstimates[21, :], lab = "F2 estimate", xlab = "time (h)", ylab = "temperature (deg C)")
scatter!(timesteps, Measurements[10, :])

plot(GT_times, GT_meas_mat[13, :], lab = "AC1", xlab = "time (h)", ylab = "temperature (deg C)")
plot!(timesteps, Normal_StateEstimates[31, :], lab = "Normal estimate", xlab = "time (h)", ylab = "temperature (deg C)")
plot!(timesteps, F1_StateEstimates[31, :], lab = "F1 estimate", xlab = "time (h)", ylab = "temperature (deg C)")
plot!(timesteps, F2_StateEstimates[31, :], lab = "F2 estimate", xlab = "time (h)", ylab = "temperature (deg C)")
scatter!(timesteps, Measurements[13, :])

plot(GT_times, GT_meas_mat[15, :], lab = "AC2", xlab = "time (h)", ylab = "temperature (deg C)")
plot!(timesteps, Normal_StateEstimates[41, :], lab = "Normal estimate", xlab = "time (h)", ylab = "temperature (deg C)")
plot!(timesteps, F1_StateEstimates[41, :], lab = "F1 estimate", xlab = "time (h)", ylab = "temperature (deg C)")
plot!(timesteps, F2_StateEstimates[41, :], lab = "F2 estimate", xlab = "time (h)", ylab = "temperature (deg C)")
scatter!(timesteps, Measurements[15, :])

plot(GT_times, GT_meas_mat[17, :], lab = "AC3", xlab = "time (h)", ylab = "temperature (deg C)")
plot(GT_times, GT_meas_mat[32, :], lab = "Tank 50", xlab = "time (h)", ylab = "temperature (deg C)")
plot(GT_times, GT_meas_mat[34, :], lab = "AC4", xlab = "time (h)", ylab = "temperature (deg C)")

plot(GT_times, GT_states_mat[2, :], lab = "Tank 10", xlab = "time (h)", ylab = "solid fraction")
plot(GT_times, GT_meas_mat[12, :], lab = "V_9", xlab = "time (h)", ylab = "flow rate (L/h)")

plot(GT_times, GT_meas_mat[37, :], lab = "Solid fraction stream 1")
plot!(timesteps, Normal_StateEstimates[118, :], lab = "Normal state estimate")
plot!(timesteps, F1_StateEstimates[118, :], lab = "F1 state estimate")
plot!(timesteps, F2_StateEstimates[118, :], lab = "F2 state estimate")

plot(GT_times, GT_states_mat[84, :], lab = "Ground truth Tank 1")
plot!(timesteps, Normal_StateEstimates[84, :], lab = "State estimate Tank 1")

plot(timesteps, PosteriorBeliefs[1, :], lab = "Normal")
plot!(timesteps, PosteriorBeliefs[2, :], lab = "Fault 4")
plot!(timesteps, PosteriorBeliefs[3, :], lab = "Fault 5")

for i in 1:115
check = i
plot(timesteps, Normal_StateEstimates[check, :], lab = "state estimate", ribbon = (Normal_x_std[check, :], Normal_x_std[check, :]))
display(plot!(GT_times, GT_states_mat[check, :], lab = "Ground truth"))
end

check = 60
plot(timesteps, Normal_StateEstimates[check, :], lab = "state estimate", ribbon = (Normal_x_std[check, :], Normal_x_std[check, :]), xlab = "time (h)", ylab = "mass fraction in liquid phase")
display(plot!(GT_times, GT_states_mat[check, :], lab = "Ground truth"))

plot(timesteps, Normal_StateEstimates[14, :], ribbon = (Normal_x_std[84, :], Normal_x_std[84, :]), lab = "state estimate")
plot!(GT_times, GT_states_mat[14,:], lab = "Ground truth")

#display(scatter!(timesteps, Measurements[17,:]))

# plot(GT_times, GT_meas_mat[10,:], lab = "tank 10 temp")
# plot(GT_times, GT_meas_mat[11,:], lab = "tank 20 temp")
# plot(GT_times, GT_meas_mat[13,:], lab = "V_9")
# plot(GT_times, GT_meas_mat[14,:], lab = "temp_AC1")

#I NEED TO LOOK AT THE CHOLESKY WEIGHTING AGAIN!