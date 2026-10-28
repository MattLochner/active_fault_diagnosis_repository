using DifferentialEquations
using Plots
using LinearAlgebra
using full_simulation
using Optim
using FiniteDiff
using Random
using ForwardDiff
#just change to see what happens when i commit

include("BMR_gt_params.jl")
using .BMR_ground_truth_params

include("BMR_kf_params.jl")
using .Kf_model_params

include("BMR_EKF_Input_Design_Functions.jl")

gtparams = data
kfparams = Kf_data

gtparams.faultmode = 5  #selected fault mode
kfparams.faultmode = 5

gtparams.flow.scaling_factor_initial = 0.8  #selected fault mode
gtparams.flow.scaling_factor_rate_change = -0.001  #selected fault mode
gtparams.flow.AC2_cooling_water_max = 10000
kfparams.flow.AC2_cooling_water_max = gtparams.flow.AC2_cooling_water_max

include("BMR_initial_conditions.jl")
#x0 = initial_conditions
#xhat0 = initial_state_estimate.*1.05  #intial_state estimate for kalman filter

#calculate initial steady state value (only applicable to non-runaway faults)

ss_params = kfparams
ss_params.auto_regressive_inputs = constant_inputs

ss_params.faultmode = 0

ss_ODE = ODEProblem(calculate_model_autoclave_odes, initial_state_estimate, (0, 20), ss_params)
ODEsol = solve(ss_ODE, Tsit5(), reltol=1e-6, abstol=1e-8)

N_steady_state = ODEsol.u[end]    #initial conditions for ground truth
x0_N = N_steady_state[1:115]
xhat0_N = N_steady_state*1.01  #intial_state estimate for kalman filter
P0_N = Matrix(diagm((xhat0_N - N_steady_state).^2))

ss_params.faultmode = 8
initial_state_estimate[119] = 0.447
ss_ODE = ODEProblem(calculate_model_autoclave_odes, initial_state_estimate, (0, 20), ss_params)
ODEsol = solve(ss_ODE, Tsit5(), reltol=1e-6, abstol=1e-8)

F8_steady_state = ODEsol.u[end]    #initial conditions for ground truth
x0_faultmode_8 = F8_steady_state[1:115]
xhat0_faultmode_8 = F8_steady_state*1.01  #intial_state estimate for kalman filter
xhat0_faultmode_8[119] = 0.8    #start offset value at a smaller and error and see if it converges to the true error
P0_faultmode_8 = Matrix(diagm((xhat0_faultmode_8 - F8_steady_state).^2))

ss_params.faultmode = 5

initial_state_estimate[119] = 1
ss_ODE = ODEProblem(calculate_model_autoclave_odes, initial_state_estimate, (0, 20), ss_params)
ODEsol = solve(ss_ODE, Tsit5(), reltol=1e-6, abstol=1e-8)

#plot(ODEsol.t, [u[41] for u in ODEsol.u]) #check if it reaches steady state

F5_steady_state = ODEsol.u[end]    #initial conditions for ground truth
x0_faultmode_5 = F5_steady_state[1:115]
xhat0_faultmode_5 = F5_steady_state*1.01  #intial_state estimate for kalman filter

xhat0_faultmode_5[119] = gtparams.flow.scaling_factor_initial #start cooling factor estimation at correct value

P0_faultmode_5 = Matrix(diagm((xhat0_faultmode_5 - F5_steady_state).^2))

#display(ODEsol.u[end])

estimationlength = 500       #specify the length of the state estimation period by the number of state estimates obtained at the end of the simulation.
N_meas = 0.2                 #h, specify the measurement sampling rate (must sum to 1)
tstart = 0
tend = tstart + estimationlength * N_meas
timesteps = tstart:N_meas:(estimationlength*N_meas)

#true_vk_std #specified in BMR_initial_conditions (measurement noise)
#true_vk_std = measurement_noise
true_vk_std = meas_noise_4 #much tighter measurements

#true_vk_std_ext = measurement_noise_ext
true_vk_std_ext = meas_noise_4_ext #much tighter measurements

vk_std = true_vk_std
vk_std_extended = true_vk_std_ext

nx = length(xhat0_N)
ny = 36
n_samples = 13
R = Matrix(diagm(vk_std.^2))
R_ext = Matrix(diagm(vk_std_extended.^2))
#Q[119,119] = Q[119,119] * 100

#Q = Matrix(diagm((xhat0_N).^2).*1e-3)    #diagonal matrix scaled relative to variable size -> could be different for different faults?
#Q_N = Q
#Q_F1 = Q
#Q_F2 = Q

#Measured_states = [1, 11, 21, 31, 41, 51, 52, 62, 81, 97, 98, 99, 100, 101, 102, 103]
Measured_states = [21, 31, 41, 51, 81, 101, 102]
xhat0_N_altered = copy(xhat0_N)
xhat0_N_altered[Measured_states] ./= 1
xhat0_faultmode_8_altered = copy(xhat0_faultmode_8)
xhat0_faultmode_8_altered[Measured_states] ./= 1
xhat0_faultmode_5_altered = copy(xhat0_faultmode_5)
xhat0_faultmode_5_altered[Measured_states] ./= 1

Q_N = Matrix(diagm((xhat0_N_altered*0.01).^2))
Q_F1 = Matrix(diagm((xhat0_faultmode_8_altered*0.01).^2))
Q_F2 = Matrix(diagm((xhat0_faultmode_5_altered*0.01).^2))
Q_F1[119,119] = Q_F1[119,119]*10
#Bayesian update of fault likelihoods
Priors = [0.333, 0.333, 0.334]

#initial ground truth values
# x0_N = x0_N
# x0_faultmode_8 = x0_N
# x0_faultmode_5 = x0_N

struct X0
    N
    F1
    F2
end

x0 = X0(x0_N, x0_faultmode_8, x0_faultmode_5)

#initial state estimates
struct Xhat0
    N
    F1
    F2
end

# xhat0_N = xhat0_N
# xhat0_faultmode_8 = xhat0_N
# xhat0_faultmode_5 = xhat0_N

xhat0 = Xhat0(xhat0_N, xhat0_faultmode_8, xhat0_faultmode_5)

#initial covariance estimates
struct p0
    N
    F1
    F2
end

# P0_N = P0_N
# P0_faultmode_8 = P0_N
# P0_faultmode_5 = P0_N

P0 = p0(P0_N, P0_faultmode_8, P0_faultmode_5)

struct q
    N
    F1
    F2
end

Q = q(Q_N, Q_F1, Q_F2)

#Run the kalman filter

function run_ODEs(realparams, x0, true_vk_std, true_vk_std_ext, Q, R, R_ext, P0, nx, ny, n_samples, estimationlength, tstart, tend, timesteps, generate_measurements, generate_sample_measurements, GT_measurement_prediction, ground_truth_ODE, xhat0, kalman_ODE, kalmanparams, Priors) 
    
    control_horizon    = 1
    prediction_horizon = 3

    afd_windows = [(500, 501)] #One month, with one shift per week with AFD.
    exit_kalmanparams = deepcopy(kalmanparams)
    afd_idx           = 0
    afd_tEnd          = 0.0
    afd_u0_initial    = 0.0
    inloop_u_opt      = Float64[]
    inloop_J          = NaN

    #Update_to_match selected process
    x0_init = x0.F2
    
    Real_Process = zeros(length(x0_init),estimationlength+1)
    Real_Process[:,1] = x0_init                                #store initial values
    yk_all_measurements1 = generate_sample_measurements(x0_init, realparams, tstart, true_vk_std_ext) #initial measurement should maybe be set to xhat0 not x0?
    Measurements = zeros(ny, estimationlength+1)
    Measurements[:,1] = yk_all_measurements1[1:36]              #initial measurement should maybe be set to xhat0 not x0?
    Sample_Measurements = zeros(n_samples, estimationlength+1)
    Sample_Measurements[:,1] = yk_all_measurements1[37:49] 

    kalmanparams.faultmode = 0

    Normal_StateEstimates = zeros(length(xhat0.N), estimationlength+1)
    Normal_StateEstimates[:,1] = xhat0.N
    Normal_CovarianceEstimates = fill(zeros(nx,nx), estimationlength+1)
    Normal_CovarianceEstimates[1] = P0.N
    Normal_Measurement_Estimates = zeros(ny+n_samples, estimationlength+1)
    Normal_Measurement_Estimates[:,1] = BMR_kf_measurement_prediction_with_samples(xhat0.N, kalmanparams, tstart)
    Normal_Measurement_Covariance = fill(zeros(ny+n_samples, ny+n_samples), estimationlength+1)
    Normal_Measurement_Posterior_Estimates = zeros(ny+n_samples, estimationlength+1)
    Normal_Measurement_Posterior_Estimates[:,1] = BMR_kf_measurement_prediction_with_samples(xhat0.N, kalmanparams, tstart)
    Normal_Measurement_Posterior_Covariance = fill(zeros(ny+n_samples, ny+n_samples), estimationlength+1)
        kalmanparams.faultmode = 8

    F1_StateEstimates = zeros(nx, estimationlength+1)
    F1_StateEstimates[:,1] = xhat0.F1
    F1_CovarianceEstimates = fill(zeros(nx,nx), estimationlength+1)
    F1_CovarianceEstimates[1] = P0.F1
    F1_Measurement_Estimates = zeros(ny+n_samples, estimationlength+1)
    F1_Measurement_Estimates[:,1] = BMR_kf_measurement_prediction_with_samples(xhat0.F1, kalmanparams, tstart)
    F1_Measurement_Covariance = fill(zeros(ny+n_samples, ny+n_samples), estimationlength+1)
    F1_Measurement_Posterior_Estimates = zeros(ny+n_samples, estimationlength+1)
    F1_Measurement_Posterior_Estimates[:,1] = BMR_kf_measurement_prediction_with_samples(xhat0.F1, kalmanparams, tstart)
    F1_Measurement_Posterior_Covariance = fill(zeros(ny+n_samples, ny+n_samples), estimationlength+1)

    kalmanparams.faultmode = 5

    F2_StateEstimates = zeros(nx, estimationlength+1)
    F2_StateEstimates[:,1] = xhat0.F2
    F2_CovarianceEstimates = fill(zeros(nx,nx), estimationlength+1)
    F2_CovarianceEstimates[1] = P0.F2
    F2_Measurement_Estimates = zeros(ny+n_samples, estimationlength+1)
    F2_Measurement_Estimates[:,1] = BMR_kf_measurement_prediction_with_samples(xhat0.F2, kalmanparams, tstart)
    F2_Measurement_Covariance = fill(zeros(ny+n_samples, ny+n_samples), estimationlength+1)   
    F2_Measurement_Posterior_Estimates = zeros(ny+n_samples, estimationlength+1)
    F2_Measurement_Posterior_Estimates[:,1] = BMR_kf_measurement_prediction_with_samples(xhat0.F2, kalmanparams, tstart)
    F2_Measurement_Posterior_Covariance = fill(zeros(ny+n_samples, ny+n_samples), estimationlength+1)  

    PosteriorBeliefs = zeros(3, estimationlength+1) #first entry is number of fault modes
    PosteriorBeliefs[:,1] = Priors 
        
    #Set up initial variables that will be modified in the loop:

    #for steady state faults, should i 

    N_xhatmin1 = xhat0.N
    F1_xhatmin1 = xhat0.F1
    F2_xhatmin1 = xhat0.F2
    N_Pposmin1 = P0.N
    F1_Pposmin1 = P0.F1
    F2_Pposmin1 = P0.F2
    Priors = Priors
  
    # preallocate before the loop:
    Applied_Inputs = zeros(2, estimationlength+1)
    Applied_Inputs[1,1] = realparams.control.setpoint[8]
    display("checking Applied_Inputs[1,1]")
    display(Applied_Inputs[1,1])
    Applied_Inputs[2,1] = realparams.control.setpoint[9]

    #initialise ground truth trajectory

    GT_times = [float(tstart)]          #grows each step
    GT_states = [copy(x0_init)]                 #vector of state vectors
    GT_meas   = [GT_measurement_prediction(x0_init, realparams, float(tstart))]

    H2SO4_AC1_cumulative_error = 0.0
    H2SO4_AC4_cumulative_error = 0.0
    supervisory_control_setpoints = zeros(4, estimationlength+1)
    supervisory_control_setpoints[:,1] = [realparams.control.setpoint[21],  # AC1 acid SP
                                        realparams.control.setpoint[13],  # stream 23 flow SP
                                        realparams.control.setpoint[22],  # AC4 acid SP
                                        realparams.control.setpoint[16]]  # stream 20 flow SP    delta_T = 1.0 #timestep between operator setpoint changes
    delta_T = 1.0

    for i in 2:(estimationlength+1)                           #the first algorithm calculation occurs after the first time step has taken place, i = 1 would call the initial values of the system
        
        #Set up equations for all fault modes
        println("i = $i")

        # if i > 500
        #     realparams.faultmode = 8    #introduce change in faultmode at t = 20
        # end
        
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

        yk_all_measurements2 = generate_sample_measurements(Real_Process[:,i], realparams, timesteps[i], true_vk_std_ext)
        yk_online = yk_all_measurements2[1:36]
        Measurements[:,i] = yk_online
                
        if isinteger(timesteps[i])
            yk_sample = yk_all_measurements2[37:49]
            yk = [yk_online; yk_sample]
            Sample_Measurements[:,i] = yk_sample
        else
            yk = yk_online
            Sample_Measurements[:,i] .= NaN   # so it's obvious nothing was sampled
        end
        
        #PROPOGATION OF UNCERTAINTY + STATE ESTIMATION WITH KALMAN FILTER
        #****************************************************************
        tStart = timesteps[i-1]
        tEnd = timesteps[i]
 
        #Normal Mode
            #set fault mode = 0
            kalmanparams.faultmode = 0
        
        if isinteger(timesteps[i])
            n_total = ny + n_samples
            (N_yhat, N_Py, N_xhatpos, N_Ppos) = ekf(N_xhatmin1, N_Pposmin1, nx, n_total, R_ext, Q.N, kalmanparams, tStart, tEnd, BMR_kf_measurement_prediction_with_samples, kalman_ODE, yk)
            N_yhat_store = N_yhat
            N_Py_store   = N_Py
        else
            (N_yhat, N_Py, N_xhatpos, N_Ppos) = ekf(N_xhatmin1, N_Pposmin1, nx, ny, R, Q.N, kalmanparams, tStart, tEnd, BMR_kf_measurement_prediction, kalman_ODE, yk)
            N_yhat_store = vcat(N_yhat, fill(NaN, n_samples))
            N_Py_store   = fill(NaN, ny+n_samples, ny+n_samples)
            N_Py_store[1:ny, 1:ny] .= N_Py
        end
 
        Normal_StateEstimates[:,i] = N_xhatpos
        Normal_CovarianceEstimates[i] = N_Ppos
        Normal_Measurement_Estimates[:,i] = N_yhat_store
        Normal_Measurement_Covariance[i] = N_Py_store
        Normal_Measurement_Posterior_Estimates[:,i] = BMR_kf_measurement_prediction_with_samples(Normal_StateEstimates[:,i], kalmanparams, timesteps[i])
        f(x) = BMR_kf_measurement_prediction_with_samples(x, kalmanparams, timesteps[i])
        J_N = ForwardDiff.jacobian(f, Normal_StateEstimates[:,i]) #linearised measurement function
        Normal_Measurement_Posterior_Covariance[i] = J_N*Normal_CovarianceEstimates[i]*J_N'
 
        N_xhatmin1 = N_xhatpos
        N_Pposmin1 = N_Ppos
 
        #Faultmode 1
            #set fault mode = 4
            kalmanparams.faultmode = 8
        if isinteger(timesteps[i])
            n_total = ny + n_samples
            (F1_yhat, F1_Py, F1_xhatpos, F1_Ppos) = ekf(F1_xhatmin1, F1_Pposmin1, nx, n_total, R_ext, Q.F1, kalmanparams, tStart, tEnd, BMR_kf_measurement_prediction_with_samples, kalman_ODE, yk)
            F1_yhat_store = F1_yhat
            F1_Py_store   = F1_Py
        else
            (F1_yhat, F1_Py, F1_xhatpos, F1_Ppos) = ekf(F1_xhatmin1, F1_Pposmin1, nx, ny, R, Q.F1, kalmanparams, tStart, tEnd, BMR_kf_measurement_prediction, kalman_ODE, yk)
            F1_yhat_store = vcat(F1_yhat, fill(NaN, n_samples))
            F1_Py_store   = fill(NaN, ny+n_samples, ny+n_samples)
            F1_Py_store[1:ny, 1:ny] .= F1_Py        
        end
        F1_StateEstimates[:,i] = F1_xhatpos
        F1_CovarianceEstimates[i] = F1_Ppos
        F1_Measurement_Estimates[:,i] = F1_yhat_store
        F1_Measurement_Covariance[i] = F1_Py_store
        F1_Measurement_Posterior_Estimates[:,i] = BMR_kf_measurement_prediction_with_samples(F1_StateEstimates[:,i], kalmanparams, timesteps[i])
        f(x) = BMR_kf_measurement_prediction_with_samples(x, kalmanparams, timesteps[i])
        J_F1 = ForwardDiff.jacobian(f, F1_StateEstimates[:,i]) #linearised measurement function
        F1_Measurement_Posterior_Covariance[i] = J_F1*F1_CovarianceEstimates[i]*J_F1' 

        F1_xhatmin1 = F1_xhatpos
        F1_Pposmin1 = F1_Ppos
 
        #Faultmode 2
            #set fault mode = 5
            kalmanparams.faultmode = 5
        if isinteger(timesteps[i])
            n_total = ny + n_samples     
            (F2_yhat, F2_Py, F2_xhatpos, F2_Ppos) = ekf(F2_xhatmin1, F2_Pposmin1, nx, n_total, R_ext, Q.F2, kalmanparams, tStart, tEnd, BMR_kf_measurement_prediction_with_samples, kalman_ODE, yk)
            F2_yhat_store = F2_yhat
            F2_Py_store   = F2_Py        
        else
            (F2_yhat, F2_Py, F2_xhatpos, F2_Ppos) = ekf(F2_xhatmin1, F2_Pposmin1, nx, ny, R, Q.F2, kalmanparams, tStart, tEnd, BMR_kf_measurement_prediction, kalman_ODE, yk)
            F2_yhat_store = vcat(F2_yhat, fill(NaN, n_samples))
            F2_Py_store   = fill(NaN, ny+n_samples, ny+n_samples)
            F2_Py_store[1:ny, 1:ny] .= F2_Py         
        end
        F2_StateEstimates[:,i] = F2_xhatpos
        F2_CovarianceEstimates[i] = F2_Ppos
        F2_Measurement_Estimates[:,i] = F2_yhat_store
        F2_Measurement_Covariance[i] = F2_Py_store
        F2_Measurement_Posterior_Estimates[:,i] = BMR_kf_measurement_prediction_with_samples(F2_StateEstimates[:,i], kalmanparams, timesteps[i])
        f(x) = BMR_kf_measurement_prediction_with_samples(x, kalmanparams, timesteps[i])
        J_F2 = ForwardDiff.jacobian(f, F2_StateEstimates[:,i]) #linearised measurement function
        F2_Measurement_Posterior_Covariance[i] = J_F2*F2_CovarianceEstimates[i]*J_F2'

        F2_xhatmin1 = F2_xhatpos
        F2_Pposmin1 = F2_Ppos 

        #BAYESIAN UPDATE OF FAULT LIKELIHOODS

        if i < 250    #allows an initial estimation phase to be specified before model mode probabilities are calculated.      
            PosteriorBeliefs[:,i] = PosteriorBeliefs[:,i-1]       
        else

            # what about adding a "refresh every 8 hours to make sure that faults are continously checked for...
            # if isinteger(timesteps(i)/8)
            #     PosteriorBeliefs[:,i-1] = Priors
            # end

            yhats = [N_yhat'; F1_yhat'; F2_yhat']
            covariances = [N_Py, F1_Py, F2_Py]

            # flattened_prior = 
            #flattened_prior = flatten_prior(PosteriorBeliefs[:,(i-1)], 0.97)
            flattened_prior = linearflattening(PosteriorBeliefs[:,(i-1)], 0.013)

            PosteriorBeliefs[:,i] = multiple_model_method(yk, flattened_prior, yhats, covariances)
        end

        #SUPERVISORY CONTROL

            #sample density of stream 1
                #update ratio of stream 1:(2+3)

            #sample acid in AC1
                #update flow rate of stream 23

            if isinteger(timesteps[i])
            
                    # if timesteps[i] < 24.99
                    #     realparams.control.setpoint[21] = 21.5
                    # elseif timesteps[i] < 49.99
                    #     realparams.control.setpoint[21] = 21.5
                    # else
                    #     realparams.control.setpoint[21] = 21.5
                    # end

                Kc_H2SO4_AC1 = 5.0
                Ti_H2SO4_AC1 = 15.0

                AC1_H2SO4_concentration = Sample_Measurements[1, i]
                AC1_H2SO4_error = realparams.control.setpoint[21] - AC1_H2SO4_concentration
                realparams.control.setpoint[13] = Kc_H2SO4_AC1 * (AC1_H2SO4_error + 1/Ti_H2SO4_AC1 * H2SO4_AC1_cumulative_error) + 40.0   

                if realparams.control.setpoint[13] < 0.0
                    realparams.control.setpoint[13] = 0.0
                    AC1_H2SO4_error = 0.0
                end

                kalmanparams.control.setpoint[13] = realparams.control.setpoint[13]
                
                H2SO4_AC1_cumulative_error = H2SO4_AC1_cumulative_error + AC1_H2SO4_error*delta_T #
            end
            #sample acid in AC4
            if isinteger(timesteps[i])
            
                    # if timesteps[i] < 39.99
                    #     realparams.control.setpoint[22] = 40.0
                    # elseif timesteps[i] < 79.99
                    #     realparams.control.setpoint[22] = 45.0
                    # else
                    #     realparams.control.setpoint[22] = 40.0
                    # end

                Kc_H2SO4_AC4 = 2.5
                Ti_H2SO4_AC4 = 7.6

                AC4_H2SO4_concentration = Sample_Measurements[4, i]
                AC4_H2SO4_error = realparams.control.setpoint[22] - AC4_H2SO4_concentration
                realparams.control.setpoint[16] = Kc_H2SO4_AC4 * (AC4_H2SO4_error + 1/Ti_H2SO4_AC4 * H2SO4_AC4_cumulative_error) + 22.0   

                if realparams.control.setpoint[16] < 0.0
                    realparams.control.setpoint[16] = 0.0
                    AC4_H2SO4_error = 0.0
                end

                kalmanparams.control.setpoint[16] = realparams.control.setpoint[16]
                
                H2SO4_AC4_cumulative_error = H2SO4_AC4_cumulative_error + AC4_H2SO4_error*delta_T #
            end

            #record setpoints at every timestep

            supervisory_control_setpoints[:,i] = [realparams.control.setpoint[21],
                                                  realparams.control.setpoint[13],
                                                  realparams.control.setpoint[22],
                                                  realparams.control.setpoint[16]]
        #INPUT DESIGN SECTION
        #---------------------------------------------------------------------
        #Option 1: Step change 
        #---------------------------------------------------------------------
        # if i < 500
        #     u0 = 130.0
        # else 
        #     u0 = 125.0
        # end

        # realparams.control.setpoint[8] = u0[1, 1]
        # kalmanparams.control.setpoint[8] = u0[1, 1]
        # Applied_Inputs[1,i] = u0[1,1]
        # display(u0[1,1])  
        #-----------------------------------------------------------------------
        #Option 2: Lets try get the input design work work - okay but how are we going to switch between number of input modes easily?
        #-----------------------------------------------------------------------
        run_afd = any(lo <= i < hi for (lo, hi) in afd_windows)
        
        if run_afd
            #commented out code allows to change to having 2 inputs as part of input design 
            Q_modes = [Q.N, Q.F1, Q.F2]

            # --- snapshot for the out-of-loop replication ---
            exit_kalmanparams = deepcopy(kalmanparams)
            afd_idx         = i
            afd_u0_initial    = Applied_Inputs[1,i-1]
            afd_tEnd = tEnd
            # -----------------------------------------------            

            idp = InputDesignParams(tEnd, N_meas, control_horizon, prediction_horizon, [N_xhatpos, F1_xhatpos, F2_xhatpos], [N_Ppos, F1_Ppos, F2_Ppos], nx, ny, R, Q_modes, kalmanparams, BMR_kf_measurement_prediction, kalman_ODE, 3, PosteriorBeliefs[:,i])
            u0_initial = Applied_Inputs[1,i-1]
            display("initial guess")
            display([u0_initial])
            u0_1 = repeat([u0_initial], control_horizon)   # length 2*Hc
            #u0_2 = repeat([125.0], control_horizon)
            lower_1 = fill(125.0, control_horizon)
            upper_1 = fill(135.0, control_horizon)
            #lower_2 = fill(120.0, control_horizon)
            #upper_2 = fill(130.0, control_horizon)

            u0 = u0_1 #vcat(u0_1, u0_2)
            lower = lower_1 #vcat(lower_1, lower_2)
            upper = upper_1 #vcat(upper_1, upper_2)

            res = Optim.optimize(
                u -> Input_Design_Function(u, idp),
                lower, upper, u0,
                Fminbox(LBFGS()),
                Optim.Options(iterations=200);
                autodiff = AutoFiniteDiff()
            )
            #display(res)
            inloop_u_opt = Optim.minimizer(res)
            inloop_J     = Optim.minimum(res)

            realparams.control.setpoint[8]   = inloop_u_opt[1]
            kalmanparams.control.setpoint[8] = inloop_u_opt[1]
            Applied_Inputs[1,i] = inloop_u_opt[1]
        else
            u0  = 130.0
            realparams.control.setpoint[8] = u0[1, 1]
            kalmanparams.control.setpoint[8] = u0[1, 1]
            Applied_Inputs[1,i] = u0[1,1]

        end
        display(Applied_Inputs[1,i])
    end

    # reshape GT_states from vector-of-vectors → matrix (n_states × n_timepoints)
    GT_states_mat = hcat(GT_states...)
    GT_meas_mat   = reduce(hcat, GT_meas)


return timesteps, GT_times, GT_states_mat, GT_meas_mat, Normal_StateEstimates, Normal_CovarianceEstimates, F1_StateEstimates, F1_CovarianceEstimates, F2_StateEstimates, F2_CovarianceEstimates, Real_Process, Measurements, PosteriorBeliefs, Normal_Measurement_Estimates, Normal_Measurement_Covariance, F1_Measurement_Estimates, F1_Measurement_Covariance, F2_Measurement_Estimates, F2_Measurement_Covariance, Applied_Inputs, Normal_Measurement_Posterior_Estimates, Normal_Measurement_Posterior_Covariance, F1_Measurement_Posterior_Estimates, F1_Measurement_Posterior_Covariance, F2_Measurement_Posterior_Estimates, F2_Measurement_Posterior_Covariance, Sample_Measurements, supervisory_control_setpoints, exit_kalmanparams, afd_idx, afd_u0_initial, inloop_u_opt, inloop_J, afd_tEnd
end

struct RUsage
    ru_utime_sec::Clong;  ru_utime_usec::Clong   # user CPU time
    ru_stime_sec::Clong;  ru_stime_usec::Clong   # system CPU time
    ru_maxrss::Clong;     ru_ixrss::Clong
    ru_idrss::Clong;      ru_isrss::Clong
    ru_minflt::Clong;     ru_majflt::Clong
    ru_nswap::Clong;      ru_inblock::Clong
    ru_oublock::Clong;    ru_msgsnd::Clong
    ru_msgrcv::Clong;     ru_nsignals::Clong
    ru_nvcsw::Clong;      ru_nivcsw::Clong
end

function cpu_time()
    ru = Ref{RUsage}()
    ccall(:getrusage, Cint, (Cint, Ref{RUsage}), 0, ru)  # 0 = RUSAGE_SELF
    r = ru[]
    utime = r.ru_utime_sec + r.ru_utime_usec / 1e6
    stime = r.ru_stime_sec + r.ru_stime_usec / 1e6
    return utime + stime   # total CPU seconds
end

function sweep_horizons_mc(snap; control_range = 1:3, prediction_range = 1:5)
    Result = @NamedTuple{control_horizon::Int, prediction_horizon::Int,
                         iterations::Int, f_calls::Int,
                         seconds::Float64, cpu_time::Float64,
                         u_opt::Vector{Float64}}

    results = Result[]

    for Hc in control_range, Hp in prediction_range
        Hc <= Hp || continue

        @info "Running" Hc Hp

        idp = InputDesignParams(snap.afd_tEnd, snap.N_meas, Hc, Hp,
                    [snap.N_x, snap.F1_x, snap.F2_x],
                    [snap.N_P, snap.F1_P, snap.F2_P],
                    snap.nx, snap.ny, snap.R, snap.Q_modes,
                    deepcopy(snap.exit_kalmanparams),
                    BMR_kf_measurement_prediction,
                    calculate_model_autoclave_odes,
                    3, snap.belief)

        u0    = fill(snap.afd_u0_initial, Hc)
        lower = fill(125.0, Hc)
        upper = fill(135.0, Hc)

        #I need to measure the CPU time over this call
        c0 = cpu_time()

        res = Optim.optimize(u -> Input_Design_Function(u, idp),
                    lower, upper, u0, Fminbox(LBFGS()),
                    Optim.Options(iterations=200); autodiff = AutoFiniteDiff())

        cpu_elapsed = 0.0
        cpu_elapsed = cpu_time() - c0
        
        #---------------------------------------------------

        push!(results, (; control_horizon = Hc, prediction_horizon = Hp,
                          iterations = Optim.iterations(res),
                          f_calls    = Optim.f_calls(res),
                          seconds    = res.time_run,
                          cpu_time   = cpu_elapsed,
                          u_opt      = copy(Optim.minimizer(res))))
    end

    return results
end

gt_params_template = gtparams
kalman_params_template = kfparams
BLAS.set_num_threads(1)
n_reps = 10
start = 71
final = 80
all_results = Vector{Any}(undef, n_reps)

for i in 71:80

    rep = i
    #build and replace autoregressive inputs
    gp = deepcopy(gt_params_template)
    kp = deepcopy(kalman_params_template)

    gp.auto_regressive_inputs = build_ar_inputs(2000 + rep)    # reseeds
    Random.seed!(1000 + rep)                                   # measurement seed, after
    # then run_ODEs(gp, ...)

    (ts_out, GT_times, GT_states_mat, GT_meas_mat, Normal_StateEstimates, Normal_CovarianceEstimates, F1_StateEstimates, F1_CovarianceEstimates, F2_StateEstimates, F2_CovarianceEstimates, Real_Process, Measurements, PosteriorBeliefs, Normal_Measurement_Estimates, Normal_Measurement_Covariance, F1_Measurement_Estimates, F1_Measurement_Covariance, F2_Measurement_Estimates, F2_Measurement_Covariance, Applied_Inputs, Normal_Measurement_Posterior_Estimates, Normal_Measurement_Posterior_Covariance, F1_Measurement_Posterior_Estimates, F1_Measurement_Posterior_Covariance, F2_Measurement_Posterior_Estimates, F2_Measurement_Posterior_Covariance, Sample_Measurements, supervisory_control_setpoints, exit_kalmanparams, afd_idx, afd_u0_initial, inloop_u_opt, inloop_J, afd_tEnd) = run_ODEs(gp, x0, true_vk_std, true_vk_std_ext, Q, R, R_ext, P0, nx, ny, n_samples, estimationlength, tstart, tend, timesteps, BMR_gt_generate_measurements, BMR_gt_generate_measurements_with_samples, BMR_gt_measurement_prediction, calculate_gt_autoclave_odes, xhat0, calculate_model_autoclave_odes, kp, Priors)
    fi = afd_idx     # = 500, the AFD step

    snap = (; afd_tEnd, N_meas, nx, ny, R,
             Q_modes = [Q.N, Q.F1, Q.F2],
             exit_kalmanparams, afd_u0_initial,
             N_x  = Normal_StateEstimates[:, fi],
             F1_x = F1_StateEstimates[:, fi],
             F2_x = F2_StateEstimates[:, fi],
             N_P  = Normal_CovarianceEstimates[fi],
             F1_P = F1_CovarianceEstimates[fi],
             F2_P = F2_CovarianceEstimates[fi],
             belief = PosteriorBeliefs[:, fi])

    all_results[i-start+1] = sweep_horizons_mc(snap)
    @info "rep done" i
end

using JLD2

jldsave("Input_design_MonteCarlo_7180.jld2"; all_results)