using DifferentialEquations
using Plots
using Random
using LinearAlgebra
using ForwardDiff
using full_simulation

Random.seed!(4)
#--------------------------
realparams = Params(6e-5, 6e-5, 1, 5e-5, 9.81, 0.0154, 0.8, 1.0, 0.6, 0.002, 1)
modelparams = Params(6e-5, 6e-5, 1, 5e-5, 9.81, 0.0154, 0.8, 1.0, 0.6, 0.002, 0)    #final parameter will be modified in a loop...
#--------------------------
#Import initial values for ground truth and simulation process
x0 = [0.5, 0.5, 0.5];
xhat_0 = [0.3, 0.3, 0.3];
#--------------------------
#Kalman filter parameters

nx = length(xhat_0)             #number of states
ny = 1                          #number of measurements
P0 = Matrix(0.001*I(nx))        #initial covariance matrix
Q = Matrix(10e-6*I(nx))         #process noise covariance matrix
true_vk_std = 0.02              #vk_std used for generating noisy measurements
vk_std = 0.02                   #standard deviation of measurement noise assumed for model
R = vk_std ^2*I(ny)             #measurment noise variance (for system with multiple measured variables, R = vk_std * I)

estimationlength = 50              #specifiy the length of the state estimation period by the number of state estimates obtained.
N_meas = 20                        #s, specify the time difference between samples
tstart = 0
tend = tstart + estimationlength * N_meas
timesteps = tstart:N_meas:(estimationlength*N_meas);
#---------------------------
#Bayesian update of fault likelihoods
Priors = [0.333, 0.333, 0.334]

#Run the Kalman filter

function UpdateBeliefs1Measurement(y_meas, PriorBeliefs, yhats, covariances)              # yhats, std_ys describe the Gaussian probability density function of p(y|x)
    ny = length(yhats)
    p_ycond_x = zeros(ny)
    marginal_likelihood = 0
    PosteriorBeliefs = zeros(ny)                            #This is wrong, yhat and Py are fault mode dependant...
    for i in 1:ny
        p_ycond_x[i] = 1 / (sqrt(covariances[i] * 2*pi)) * exp(- (y_meas - yhats[i]) ^ 2 / (2 * covariances[i]))
        marginal_likelihood = marginal_likelihood + p_ycond_x[i] * PriorBeliefs[i]
    end
    
    #Calculate updated beliefs
    for i in 1:ny
        PosteriorBeliefs[i] = p_ycond_x[i] * PriorBeliefs[i] / marginal_likelihood
    end
    return PosteriorBeliefs
end

function run_kalman_filter(realparams, x0, modelparams, xhat_0, nx, ny, P0, Q, true_vk_std, vk_std, R, estimationlength, N_meas, tstart, tend, timesteps, measurement_prediction_function, generate_measurements, ground_truth_ODE, model_ODE , Priors)

    Real_Process = zeros(nx,estimationlength+1)
    Real_Process[:,1] = x0                                #store initial values
    Measurements = zeros(ny, estimationlength+1)
    Measurements[:,1] = measurement_prediction_function(xhat_0, modelparams, tstart)              #initial measurement is set equal to the initial state estimate fed to the kalman filter

    Normal_StateEstimates = zeros(nx, estimationlength+1)
    Normal_StateEstimates[:,1] = xhat_0
    Normal_CovarianceEstimates = fill(zeros(nx,nx), estimationlength+1)
    Normal_CovarianceEstimates[1] = P0

    F1_StateEstimates = zeros(nx, estimationlength+1)
    F1_StateEstimates[:,1] = xhat_0
    F1_CovarianceEstimates = fill(zeros(nx,nx), estimationlength+1)
    F1_CovarianceEstimates[1] = P0

    F2_StateEstimates = zeros(nx, estimationlength+1)
    F2_StateEstimates[:,1] = xhat_0
    F2_CovarianceEstimates = fill(zeros(nx,nx), estimationlength+1)
    F2_CovarianceEstimates[1] = P0

    PosteriorBeliefs = zeros(3, estimationlength+1) #first entry is number of fault modes
    PosteriorBeliefs[:,1] = Priors

    #Set up initial variables that will be modified in the loop:

    N_xhatmin1 = xhat_0
    F1_xhatmin1 = xhat_0
    F2_xhatmin1 = xhat_0
    N_Pposmin1 = P0
    F1_Pposmin1 = P0
    F2_Pposmin1 = P0
    Priors = Priors

    #store initial input flow rates

    # store_u1_values = [realparams.u1]
    # store_u2_values = [realparams.u2]
    # control_inputs_u1 = fill(realparams.u1, control_horizon)    #initial control input vector
    # control_inputs_u2 = fill(realparams.u2, control_horizon)    #initial control input vector

    for i in 2:(estimationlength+1)                           #the first algorithm calculation occurs after the first time step has taken place, i = 1 would call the initial values of the system
        
        #Set up equations for all fault modes
        
        gt_ODE = ODEProblem(ground_truth_ODE, Real_Process[:,i-1], (timesteps[i-1], timesteps[i]), realparams)
        
        #Only solve for fault mode that is selected to be present in the real process.

        #Change the ODEProblem called in this line to which fault mode is being simulated

        ODEsol = solve(gt_ODE, Tsit5(), save_everystep = false, save_start = false) 

        #Save simulated real process data
        Real_Process[:,i] = ODEsol.u[1]

        #Generate measurements

        yk = generate_measurements(Real_Process, i, true_vk_std)
        Measurements[:,i] = yk

        #PROPOGATION OF UNCERTAINTY + STATE ESTIMATION WITH KF
        #****************************************************

        tStart = timesteps[i-1]
        tEnd = timesteps[i]

        #NormalMode 
            #set fault mode = 0
            modelparams.faultmode = 0

        (N_yhat, N_Py, N_xhatpos, N_Ppos) = ekf(N_xhatmin1, N_Pposmin1, nx, ny, R, Q, modelparams, tStart, tEnd, measurement_prediction_function, model_ODE, yk)
        
        Normal_StateEstimates[:,i] = N_xhatpos
        Normal_CovarianceEstimates[i] = N_Ppos

        N_xhatmin1 = N_xhatpos
        N_Pposmin1 = N_Ppos

        #Faultmode 1
            #set fault mode = 1
            modelparams.faultmode = 1

        (F1_yhat, F1_Py, F1_xhatpos, F1_Ppos) = ekf(F1_xhatmin1, F1_Pposmin1, nx, ny, R, Q, modelparams, tStart, tEnd, measurement_prediction_function, model_ODE, yk)
        F1_StateEstimates[:,i] = F1_xhatpos
        F1_CovarianceEstimates[i] = F1_Ppos

        F1_xhatmin1 = F1_xhatpos
        F1_Pposmin1 = F1_Ppos  

        #Faultmode 2

            #set fault mode = 2
            modelparams.faultmode = 2

        (F2_yhat, F2_Py, F2_xhatpos, F2_Ppos) = ekf(F2_xhatmin1, F2_Pposmin1, nx, ny, R, Q, modelparams, tStart, tEnd, measurement_prediction_function, model_ODE, yk)

        F2_StateEstimates[:,i] = F2_xhatpos
        F2_CovarianceEstimates[i] = F2_Ppos

        F2_xhatmin1 = F2_xhatpos
        F2_Pposmin1 = F2_Ppos  

        #BAYESIAN UPDATE OF FAULT LIKELIHOODS
        #*****************************************************

        yhats = [N_yhat[1,1], F1_yhat[1,1], F2_yhat[1,1]]
        covariances = [N_Py[1,1], F1_Py[1,1], F2_Py[1,1]]   #for one measurement, covariance matrices are required for multiple measurements

        PosteriorBeliefs[:,i] = UpdateBeliefs1Measurement(Measurements[i], PosteriorBeliefs[:,(i-1)], yhats, covariances)   

    end
return timesteps, Normal_StateEstimates, Normal_CovarianceEstimates, F1_StateEstimates, F1_CovarianceEstimates, F2_StateEstimates, F2_CovarianceEstimates, Real_Process, Measurements, PosteriorBeliefs
end

(timesteps, Normal_StateEstimates, Normal_CovarianceEstimates, F1_StateEstimates, F1_CovarianceEstimates, F2_StateEstimates, F2_CovarianceEstimates, Real_Process, Measurements, PosteriorBeliefs) = run_kalman_filter(realparams, x0, modelparams, xhat_0, nx, ny, P0, Q, true_vk_std, vk_std, R, estimationlength, N_meas, tstart, tend, timesteps, measurement_prediction_3_tank_system, generate_measurements_3_tank_system, Three_tank_GT_process!, Three_tank_model_process!, Priors)

function CalcStdDevs(x, P)
    x_std = fill(0.0, size(x))

    for i in 1:size(x)[1]
        global x_std[i,:] = sqrt.([P[j][i,i] for j in 1:length(P)])
    end
    return x_std
end

#I'm only going to do this for specific states i want to plot, else it's too much...

Normal_x_std = CalcStdDevs(Normal_StateEstimates, Normal_CovarianceEstimates)
F1_x_std = CalcStdDevs(F1_StateEstimates, F1_CovarianceEstimates)
F2_x_std = CalcStdDevs(F2_StateEstimates, F2_CovarianceEstimates)

plot(timesteps, Real_Process[3,:])
plot!(timesteps, Normal_StateEstimates[3,:], ribbon = (Normal_x_std[3, :], Normal_x_std[3, :]), lab = "Normal")
plot!(timesteps, F1_StateEstimates[3,:], ribbon = (F1_x_std[3, :], F1_x_std[3, :]), lab = "F1")
plot!(timesteps, F2_StateEstimates[3,:], ribbon = (F2_x_std[3, :], F2_x_std[3, :]), lab = "F2")
display(scatter!(timesteps, Measurements'))

plot(timesteps, Real_Process[2,:], lab = "Ground truth")
plot!(timesteps, Normal_StateEstimates[2,:], ribbon = (Normal_x_std[2, :], Normal_x_std[2, :]), lab = "Normal")
plot!(timesteps, F1_StateEstimates[2,:], ribbon = (F1_x_std[2, :], F1_x_std[2, :]), lab = "F1")
plot!(timesteps, F2_StateEstimates[2,:], ribbon = (F2_x_std[2, :], F2_x_std[2, :]), lab = "F2")

plot(timesteps, Real_Process[1,:], lab = "Ground truth")
plot!(timesteps, Normal_StateEstimates[1,:], ribbon = (Normal_x_std[1, :], Normal_x_std[1, :]), lab = "Normal")
plot!(timesteps, F1_StateEstimates[1,:], ribbon = (F1_x_std[1, :], F1_x_std[1, :]), lab = "F1")
plot!(timesteps, F2_StateEstimates[1,:], ribbon = (F2_x_std[1, :], F2_x_std[1, :]), lab = "F2")

plot(timesteps, PosteriorBeliefs[1,:])
plot!(timesteps, PosteriorBeliefs[2,:])
plot!(timesteps, PosteriorBeliefs[3,:])