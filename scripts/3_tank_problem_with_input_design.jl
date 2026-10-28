using DifferentialEquations
using Plots
using Random
using LinearAlgebra
using ForwardDiff
using full_simulation
using Optim
using FiniteDiff
using Optimization
using OptimizationIpopt

Random.seed!(42)

#--------------------------
realparams = Params(6e-5, 6e-5, 1, 5e-5, 9.81, 0.0154, 0.8, 1.0, 0.6, 0.002, 1)
modelparams = Params(6e-5, 6e-5, 1, 5e-5, 9.81, 0.0154, 0.8, 1.0, 0.6, 0.002, 1)    #final parameter will be modified in a loop...
#--------------------------
#Import initial values for ground truth and simulation process
x0 = [0.5, 0.5, 0.5];
xhat_0 = [0.45, 0.44, 0.46];
#--------------------------
#Kalman filter parameters

nx = length(xhat_0)             #number of states
ny = 1                          #number of measurements
P0 = Matrix(0.001*I(nx))        #initial covariance matrix
Q = Matrix(10e-6*I(nx))         #process noise covariance matrix
true_vk_std = 0.02              #vk_std used for generating noisy measurements
vk_std = 0.02                   #standard deviation of measurement noise assumed for model
R = vk_std ^2*I(ny)             #measurment noise variance (for system with multiple measured variables, R = vk_std * I)

estimationlength = 50                #specifiy the length of the state estimation period by the number of state estimates obtained.
N_meas = 20.0                        #s, specify the time difference between samples
tstart = 0.0                         #s, specify the start time of the simulation
tend = tstart + estimationlength * N_meas
timesteps = tstart:N_meas:(estimationlength*N_meas);
#---------------------------
#Bayesian update of fault probabilities
Priors = [0.333, 0.333, 0.334]        #specify the prior fault mode probabilities

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

function ExtendedStatePrediction(xhatmin1, Pposmin1, nx, ny, R, Q, params, tstart, tend, measurement_pred_function, model_ODE_function)

    #calculate prior state estimate
    xhatneg = run_ODE_model(tstart, tend, xhatmin1, params, model_ODE_function)
    #----------------------------------------------------------------------
    # f = x -> begin
    #     dx = similar(x)          # matches eltype(x), so it accepts Duals
    #     model_ODE_function(dx, x, params, tstart)
    #     return dx
    # end

    # #approximate Jacobian at current state estimate

    # A = ForwardDiff.jacobian(f, xhatneg)
    #----------------------------------------------------------------------
    f(x) = model_ODE_function(x, params, tend)
    A = FiniteDiff.finite_difference_jacobian(f, xhatneg)
    #----------------------------------------------------------------------

    #calculate the prior covariance matrix

    #Option 1: ODE solver
#----------------------------------------------------------------
    # prob = ODEProblem(lyap!, Pposmin1, (tstart, tend), (A, Q))
    # sol = solve(prob, Tsit5(), saveat = tend)

    # Pneg_1 = sol.u[2]

    # Pneg_1 = (Pneg_1 + Pneg_1') / 2

#--------------------------------------------------------------
    #Option 2: "Van Loan" trick
    n = size(A, 1)

    dt = tend-tstart

    M = [A*dt              Q*dt;
         zeros(n, n)     -A'*dt]

    E = exp(M)

    Φ = E[1:n, 1:n]         # = exp(A)
    G = E[1:n, n+1:end]     # = ∫₀¹ exp(A s) Q exp(A' s) ds  (times Φ', see below)
    
    #Pneg_2 = Φ * Pposmin1 * Φ' + Φ * G
    Pneg_2 = Φ * Pposmin1 * Φ' + G * Φ'
#-------------------------------------------------------------------
    # diagnose + symmetrize
    #debugging
    
    #Pneg_2, λmin2, asym2 = check_covariance(Pneg_2; name="Pneg_2")
    #Pneg_1, λmin1, asym1 = check_covariance(Pneg_1; name="Pneg_1")

    # only assert true positive-definiteness on the symmetrized version
    #@assert isposdef(Symmetric(Pneg_1)) "Pneg_1 not PD even after symmetrizing (λmin=$λmin1)"
    #@assert isposdef(Symmetric(Pneg_2)) "Pneg_2 not PD even after symmetrizing (λmin=$λmin2)"

    return xhatneg, Pneg_2         #it must just return the prior state estimate and posterior state estimate so that Bayesian update of probability can work
end

function LinearisedMeasurementPrediction(xhatneg, Pneg, params, R, measurement_pred_function, tend)
    #Measurement prediction
    #------------------------------------------------------
        h(x) = measurement_pred_function(x, params, tend)
        H = FiniteDiff.finite_difference_jacobian(h, xhatneg)
        yhat = measurement_pred_function(xhatneg, params, tend)
    #------------------------------------------------------
        Py = H*Pneg*H'+ R #assumption simplififying MRM to R

    return yhat, Py
end 

struct InputDesignParams{TR, TQ, MF, OF}
    t_start            :: Float64
    N_meas             :: Float64
    control_horizon    :: Int
    prediction_horizon :: Int
    xhat0_modes        :: Vector{Vector{Float64}}
    P0_modes           :: Vector{Matrix{Float64}}
    nx                 :: Int
    ny                 :: Int
    R                  :: TR          # Diagonal{Float64} from vk_std^2*I(ny)
    Q                  :: TQ          # Matrix{Float64}
    propparams         :: Params
    measurement_pred_function :: MF
    model_ODE_function        :: OF
    no_of_modes        :: Int
    beliefs            :: Vector{Float64}
end

Base.copy(q::Params) = Params(q.u1, q.u2, q.c1, q.Sp, q.g, q.A,
                              q.c2, q.c3, q.alpha, q.r, q.faultmode)

function Input_Design_Function(u_vec, p)
    Hc, Hp = p.control_horizon, p.prediction_horizon
    U = reshape(u_vec, 2, Hc)          # column j = [u1(j); u2(j)]

    yhats = [Vector{Vector{Float64}}(undef, p.no_of_modes) for _ in 1:Hp]  # [t][mode] -> mean
    Pys   = [Vector{Matrix{Float64}}(undef, p.no_of_modes) for _ in 1:Hp]  # [t][mode] -> cov
    params = copy(p.propparams)

    for m in 1:p.no_of_modes
        params.faultmode = m - 1             # mode 1,2,3 -> faultmode 0,1,2

        x = copy(p.xhat0_modes[m])           # start from THIS mode's posterior
        P = copy(p.P0_modes[m])

        for j in 1:Hp
            col = j ≤ Hc ? j : Hc            # move over control horizon, then hold
            params.u1 = U[1, col]
            params.u2 = U[2, col]

            tj  = p.t_start + (j-1)*p.N_meas
            tj1 = p.t_start +  j   *p.N_meas
            x, P = ExtendedStatePrediction(x, P, p.nx, p.ny, p.R, p.Q,
                                           params, tj, tj1,
                                           p.measurement_pred_function,
                                           p.model_ODE_function)
            yhats[j][m], Pys[j][m] = LinearisedMeasurementPrediction(
                                         x, P, params, p.R,
                                         p.measurement_pred_function, tj1)
        end
    end
    return Calc_Obj_Function(p.beliefs, yhats, Pys, p.no_of_modes)
end

function bhatt_dist(x, X, y, Y)
    S  = (X + Y) / 2
    Δ  = x - y
    t1 = (Δ' * (S \ Δ)) / 8
    t2 = 0.5 * (logdet(S) - 0.5*(logdet(X) + logdet(Y)))
    return t1 + t2                     # d_ij ≥ 0, NOT exp(-d)
end

function Calc_Obj_Function(beliefs, yhats, Pys, no_of_modes)
    Hp = length(yhats)                 # outer is timesteps now
    s = 0.0
    for i in 1:no_of_modes, j in (i+1):no_of_modes
        dsum = 0.0
        for t in 1:Hp
            dsum += bhatt_dist(yhats[t][i], Pys[t][i], yhats[t][j], Pys[t][j])
        end
        s += sqrt(beliefs[i] * beliefs[j]) * exp(-dsum)
    end
    return s
end

function Input_Design_Constraint!(c, u_vec, p)
    Hc, Hp = p.control_horizon, p.prediction_horizon
    U = reshape(u_vec, 2, Hc)
    params = copy(p.propparams)
    idx = 1
    for m in 1:p.no_of_modes
        params.faultmode = m - 1
        x = copy(p.xhat0_modes[m])
        for j in 1:Hp
            col = j ≤ Hc ? j : Hc
            params.u1 = U[1, col]; params.u2 = U[2, col]
            tj  = p.t_start + (j-1)*p.N_meas
            tj1 = p.t_start +  j   *p.N_meas
            x = run_ODE_model(tj, tj1, x, params, p.model_ODE_function)  # mean only
            for k in 1:p.nx
                c[idx] = x[k]; idx += 1
            end
        end
    end
    return nothing
end

function run_kalman_filter(realparams, x0, modelparams, xhat_0, nx, ny, P0, Q, true_vk_std, vk_std, R, estimationlength, N_meas, tstart, tend, timesteps, measurement_prediction_function, generate_measurements, ground_truth_ODE, model_ODE , Priors)

    control_horizon    = 1
    prediction_horizon = 5

    Real_Process = zeros(nx,estimationlength+1)
    Real_Process[:,1] = x0                                #store initial values
    GT_values = Vector{Vector{Float64}}()
    GT_times = Vector{Float64}()

    #push!(GT_values, x0)
    #push!(GT_times, [tstart])

    Measurements = zeros(ny, estimationlength+1)
    Measurements[:,1] = measurement_prediction_function(xhat_0, modelparams, tstart)              #initial measurement is set equal to the initial state estimate fed to the kalman filter

    Normal_StateEstimates = zeros(nx, estimationlength+1)
    Normal_StateEstimates[:,1] = xhat_0
    Normal_CovarianceEstimates = fill(zeros(nx,nx), estimationlength+1)
    Normal_CovarianceEstimates[1] = P0

    Normal_MeasPred_means = zeros(estimationlength+1)
    Normal_MeasPred_means[1] = xhat_0[3]
    Normal_MeasPred_covariance = zeros(estimationlength+1)
    Normal_MeasPred_covariance[1] = P0[3,3]

    F1_StateEstimates = zeros(nx, estimationlength+1)
    F1_StateEstimates[:,1] = xhat_0
    F1_CovarianceEstimates = fill(zeros(nx,nx), estimationlength+1)
    F1_CovarianceEstimates[1] = P0

    F1_MeasPred_means = zeros(estimationlength+1)
    F1_MeasPred_means[1] = xhat_0[3]
    F1_MeasPred_covariance = zeros(estimationlength+1)
    F1_MeasPred_covariance[1] = P0[3,3]

    F2_StateEstimates = zeros(nx, estimationlength+1)
    F2_StateEstimates[:,1] = xhat_0
    F2_CovarianceEstimates = fill(zeros(nx,nx), estimationlength+1)
    F2_CovarianceEstimates[1] = P0

    F2_MeasPred_means = zeros(estimationlength+1)
    F2_MeasPred_means[1] = xhat_0[3]
    F2_MeasPred_covariance = zeros(estimationlength+1)
    F2_MeasPred_covariance[1] = P0[3,3]

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

    # preallocate before the loop:
    Applied_Inputs = zeros(2, estimationlength+1)
    Applied_Inputs[:,1] = [realparams.u1, realparams.u2]

    # control_inputs_u1 = fill(realparams.u1, control_horizon)    #initial control input vector
    # control_inputs_u2 = fill(realparams.u2, control_horizon)    #initial control input vector

    for i in 2:(estimationlength+1)                           #the first algorithm calculation occurs after the first time step has taken place, i = 1 would call the initial values of the system
        
        println("i = $i")

        #Set up equations for all fault modes
        
        gt_ODE = ODEProblem(ground_truth_ODE, Real_Process[:,i-1], (timesteps[i-1], timesteps[i]), realparams)
        
        #Only solve for fault mode that is selected to be present in the real process.

        #Change the ODEProblem called in this line to which fault mode is being simulated

        ODEsol = solve(gt_ODE, Tsit5(),abstol = 1e-8, reltol = 1e-8)#, save_everystep = false, save_start = false) 

        GT_values = append!(GT_values, ODEsol.u)
        GT_times = append!(GT_times, ODEsol.t)

        #Save simulated real process data
        Real_Process[:,i] = ODEsol.u[end]

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
        Normal_MeasPred_means[i] = N_yhat[1,1]
        Normal_MeasPred_covariance[i] = N_Py[1,1]

        N_xhatmin1 = N_xhatpos
        N_Pposmin1 = N_Ppos

        #Faultmode 1
            #set fault mode = 1
            modelparams.faultmode = 1

        (F1_yhat, F1_Py, F1_xhatpos, F1_Ppos) = ekf(F1_xhatmin1, F1_Pposmin1, nx, ny, R, Q, modelparams, tStart, tEnd, measurement_prediction_function, model_ODE, yk)
        F1_StateEstimates[:,i] = F1_xhatpos
        F1_CovarianceEstimates[i] = F1_Ppos
        F1_MeasPred_means[i] = F1_yhat[1,1]
        F1_MeasPred_covariance[i] = F1_Py[1,1]

        F1_xhatmin1 = F1_xhatpos
        F1_Pposmin1 = F1_Ppos  

        #Faultmode 2

            #set fault mode = 2
            modelparams.faultmode = 2

        (F2_yhat, F2_Py, F2_xhatpos, F2_Ppos) = ekf(F2_xhatmin1, F2_Pposmin1, nx, ny, R, Q, modelparams, tStart, tEnd, measurement_prediction_function, model_ODE, yk)

        F2_StateEstimates[:,i] = F2_xhatpos
        F2_CovarianceEstimates[i] = F2_Ppos
        F2_MeasPred_means[i] = F2_yhat[1,1]
        F2_MeasPred_covariance[i] = F2_Py[1,1]

        F2_xhatmin1 = F2_xhatpos
        F2_Pposmin1 = F2_Ppos  

        #BAYESIAN UPDATE OF FAULT LIKELIHOODS
        #*****************************************************

        yhats = [N_yhat[1,1], F1_yhat[1,1], F2_yhat[1,1]]
        covariances = [N_Py[1,1], F1_Py[1,1], F2_Py[1,1]]   #for one measurement, covariance matrices are required for multiple measurements

        PosteriorBeliefs[:,i] = UpdateBeliefs1Measurement(Measurements[i], PosteriorBeliefs[:,(i-1)], yhats, covariances)

        #INPUT DESIGN
        #-----------------------------------------------------------------------------------------------
        #Option 1: Constant inputs
        #-----------------------------------------------------------------------------------------------
        # u0    = repeat([6e-5, 6e-5], control_horizon)   # length 2*Hc
        # realparams.u1,  realparams.u2  = u0[1,1], u0[2,1]
        # modelparams.u1, modelparams.u2 = u0[1,1], u0[2,1]
        # Applied_Inputs[:,i] = u0[:,1]
        
        #-----------------------------------------------------------------------------------------------
        #Option 2: Optim.optimize, input constraints, no state constraints
        #-----------------------------------------------------------------------------------------------
        # idp = InputDesignParams(tEnd, N_meas, control_horizon, prediction_horizon,
        #         [N_xhatpos, F1_xhatpos, F2_xhatpos],
        #         [N_Ppos,    F1_Ppos,    F2_Ppos],
        #         nx, ny, R, Q, modelparams,
        #         measurement_prediction_function, model_ODE, 3, PosteriorBeliefs[:,i])

        # u0    = repeat([6e-5, 6e-5], control_horizon)   # length 2*Hc
        # lower = 0.8 .* 6e-5
        # upper = 1.2 .* 6e-5
        
        # res = Optim.optimize(
        #     u -> Input_Design_Function(u, idp),
        #     lower, upper, u0,
        #     Fminbox(LBFGS()),
        #     Optim.Options(iterations=200);
        #     autodiff = AutoFiniteDiff()
        # )

        # u_opt = reshape(Optim.minimizer(res), 2, control_horizon)
        # realparams.u1,  realparams.u2  = u_opt[1,1], u_opt[2,1]
        # modelparams.u1, modelparams.u2 = u_opt[1,1], u_opt[2,1]
        # Applied_Inputs[:,i] = u_opt[:,1]
        #----------------------------------------------------------------------------------------------
        #Option 3: Option with state constraints
        #----------------------------------------------------------------------------------------------
        idp = InputDesignParams(tEnd, N_meas, control_horizon, prediction_horizon,
                [N_xhatpos, F1_xhatpos, F2_xhatpos],
                [N_Ppos,    F1_Ppos,    F2_Ppos],
                nx, ny, R, Q, modelparams,
                measurement_prediction_function, model_ODE, 3, PosteriorBeliefs[:,i])

        n_steps = idp.prediction_horizon
        n_con = n_steps * idp.no_of_modes * idp.nx

        initial_u      = repeat(Applied_Inputs[:,i-1], control_horizon)   # length 2*Hc
        u_lowerlimits  = fill(0.5*1e-4, 2*control_horizon)
        u_upperlimits  = fill(1.2*1e-4, 2*control_horizon)

        #per-state box, repeated across (step, mode)
        x_lo = [0.0, 0.0, 0.0]
        x_high = [0.6, 0.6, 0.6]
        lcons = repeat(x_lo, n_steps * idp.no_of_modes)
        ucons = repeat(x_high, n_steps * idp.no_of_modes)
        
        objfun = OptimizationFunction(Input_Design_Function, Optimization.AutoFiniteDiff();
                              cons = Input_Design_Constraint!)

        prob = Optimization.OptimizationProblem(objfun, initial_u, idp;
                                lb = u_lowerlimits, ub = u_upperlimits,
                                lcons = lcons, ucons = ucons)

        opt = IpoptOptimizer(; hessian_approximation = "limited-memory",
                            limited_memory_max_history = 10,
                            acceptable_tol  = 1e-2,
                            additional_options = Dict("acceptable_iter" => 5))

        sol = solve(prob, opt;
                    maxiters = 1000,    # NOT max_iter
                    reltol   = 1e-5,    # NOT tol  (this is Ipopt's `tol`)
                    verbose  = false)   # NOT print_level      

        u_opt = reshape(sol.u, 2, control_horizon)
        realparams.u1,  realparams.u2  = u_opt[1,1], u_opt[2,1]
        modelparams.u1, modelparams.u2 = u_opt[1,1], u_opt[2,1]
        Applied_Inputs[:,i] = u_opt[:,1]
        display(Applied_Inputs[:,i])

        #----------------------------------------------------------------------------------------------

    end
    return timesteps, Normal_StateEstimates, Normal_CovarianceEstimates, F1_StateEstimates, F1_CovarianceEstimates, F2_StateEstimates, F2_CovarianceEstimates, Real_Process, Measurements, PosteriorBeliefs, Applied_Inputs, GT_times, GT_values, Normal_MeasPred_means, Normal_MeasPred_covariance, F1_MeasPred_means, F1_MeasPred_covariance, F2_MeasPred_means, F2_MeasPred_covariance
end

(timesteps, Normal_StateEstimates, Normal_CovarianceEstimates, F1_StateEstimates, F1_CovarianceEstimates, F2_StateEstimates, F2_CovarianceEstimates, Real_Process, Measurements, PosteriorBeliefs, Applied_Inputs, GT_times, GT_values, Normal_MeasPred_means, Normal_MeasPred_covariance, F1_MeasPred_means, F1_MeasPred_covariance, F2_MeasPred_means, F2_MeasPred_covariance) = run_kalman_filter(realparams, x0, modelparams, xhat_0, nx, ny, P0, Q, true_vk_std, vk_std, R, estimationlength, N_meas, tstart, tend, timesteps, measurement_prediction_3_tank_system, generate_measurements_3_tank_system, Three_tank_GT_process!, Three_tank_model_process!, Priors)

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
display(scatter!(timesteps[2:end], Measurements'[2:end]))

plot(timesteps, Real_Process[2,:], lab = "Ground truth")
plot!(timesteps, Normal_StateEstimates[2,:], ribbon = (Normal_x_std[2, :], Normal_x_std[2, :]), lab = "Normal")
plot!(timesteps, F1_StateEstimates[2,:], ribbon = (F1_x_std[2, :], F1_x_std[2, :]), lab = "F1")
plot!(timesteps, F2_StateEstimates[2,:], ribbon = (F2_x_std[2, :], F2_x_std[2, :]), lab = "F2")

plot(timesteps, Real_Process[1,:], lab = "Ground truth")
plot!(timesteps, Normal_StateEstimates[1,:], ribbon = (Normal_x_std[1, :], Normal_x_std[1, :]), lab = "Normal")
plot!(timesteps, F1_StateEstimates[1,:], ribbon = (F1_x_std[1, :], F1_x_std[1, :]), lab = "F1")
plot!(timesteps, F2_StateEstimates[1,:], ribbon = (F2_x_std[1, :], F2_x_std[1, :]), lab = "F2")

plot(timesteps, PosteriorBeliefs[1,:], lab = "Normal")
plot!(timesteps, PosteriorBeliefs[2,:], lab = "Fault mode 1")
plot!(timesteps, PosteriorBeliefs[3,:], lab = "Fault mode 2", xlab = "time (s)", ylab = "probability", legend = :topright)

plot(timesteps, Applied_Inputs[1,:], lab = "U1", ylims = [5e-5, 7e-5])
plot!(timesteps, Applied_Inputs[2,:], lab = "U2", )

plot(timesteps./60, Applied_Inputs[1,:].*60000, xlab = "time (min)", ylab = "flow rate (l/min)", lab = "u1 values")
display(plot!(timesteps./60, Applied_Inputs[2,:]*60000, xlab = "time (min)", ylab = "flow rate (l/min)", lab = "u2 values"))

##PLOTS FOR STATE ESTIMATION RESULTS
    #Select normal operating conditions (fault mode = 0)

    GT_times
    GT_values
    GT_val_matrix = stack(GT_values)

    plot_a = plot(GT_times, GT_val_matrix[3,:], lab = "Ground truth", xlab = "time (s)", ylab = "tank level (m)")
    plot_a = plot!(timesteps, Normal_StateEstimates[3,:], ribbon = (Normal_x_std[3, :], Normal_x_std[3, :]), lab = "State Estimate")
    plot_a = scatter!(timesteps[2:end], Measurements'[2:end], lab = "Measurement")

    plot_b = plot(GT_times, GT_val_matrix[2,:], lab = "Ground truth", xlab = "time (s)", ylab = "tank level (m)")
    plot_b = plot!(timesteps, Normal_StateEstimates[2,:], ribbon = (Normal_x_std[2, :], Normal_x_std[2, :]), lab = "State Estimate")

    plot_c = plot(GT_times, GT_val_matrix[1,:], lab = "Ground truth", xlab = "time (s)", ylab = "tank level (m)")
    plot_c = plot!(timesteps, Normal_StateEstimates[1,:], ribbon = (Normal_x_std[1, :], Normal_x_std[1, :]), lab = "State Estimate")

    display(plot_a)

    # real panels — legend off
    plot_a = plot(GT_times, GT_val_matrix[3,:], lab="Ground truth",
                xlab="time (s)", ylab="tank 3 level (m)", legend=false,
                guidefontsize=16, tickfontsize=13, lw=3)
    plot!(plot_a, timesteps, Normal_StateEstimates[3,:],
        ribbon=(Normal_x_std[3,:], Normal_x_std[3,:]), lw=3)
    plot!(plot_a, timesteps, F1_StateEstimates[3,:],
        ribbon=(Normal_x_std[3,:], F1_x_std[3,:]), lw=3)
    plot!(plot_a, timesteps, F2_StateEstimates[3,:],
        ribbon=(Normal_x_std[3,:], F2_x_std[3,:]), lw=3)
    scatter!(plot_a, timesteps[2:end], Measurements'[2:end])

    plot_b = plot(GT_times, GT_val_matrix[2,:], lab="Ground truth",
                xlab="time (s)", ylab="tank 2 level (m)", legend=false,
                guidefontsize=16, tickfontsize=13, lw=3)
    plot!(plot_b, timesteps, Normal_StateEstimates[2,:],
        ribbon=(Normal_x_std[2,:], Normal_x_std[2,:]), lw=3)

    plot_c = plot(GT_times, GT_val_matrix[1,:], lab="Ground truth",
                xlab="time (s)", ylab="tank 1 level (m)", legend=false,
                guidefontsize=16, tickfontsize=13, lw=3)
    plot!(plot_c, timesteps, Normal_StateEstimates[1,:],
        ribbon=(Normal_x_std[1,:], Normal_x_std[1,:]), lw=3)
    # legend-only panel
    legend_plot = plot(framestyle=:none, legend=:inside, legendfontsize=11)
    plot!(legend_plot, [NaN], [NaN], color=1, lw=2, label="Ground truth")
    plot!(legend_plot, [NaN], [NaN], color=2, lw=2, label="State Estimate")
    scatter!(legend_plot, [NaN], [NaN], color=3, label="Measurement")

    final = plot(plot_a, plot_b, plot_c, legend_plot,
                layout=(2,2), size=(1200, 780), dpi=300,
                left_margin=5Plots.mm, bottom_margin=5Plots.mm,
                legendfontsize=13, titlefontsize=18)

    a = 2
##PLOTS FOR FAULT DIAGNOSIS RESULTS

#Plot 1: 
plot(timesteps, PosteriorBeliefs[1,:], lab = "Normal", ylims = [-0.02, 1.02])
plot!(timesteps, PosteriorBeliefs[2,:], lab = "Fault mode 1")
plot!(timesteps, PosteriorBeliefs[3,:], lab = "Fault mode 2", xlab = "time (s)", ylab = "probability", legend = :topleft)

#Plot 2: 

plot(timesteps, Applied_Inputs[1,:], lab = "U1", ylims = [4.5e-5, 9e-5], xlab = "time (s)", ylab = "input flowrate (m³/s)")
plot!(timesteps, Applied_Inputs[2,:], linestyle = :dash, lab = "U2", lw = 2)

#Plot 3: 

plot_f = plot(GT_times, GT_val_matrix[3,:], lab="Ground truth", xlab = "time (s)", ylab = "tank 3 level (m)", 
                guidefontsize=10, tickfontsize = 8, lw = 2, linestyle=:dash, color=:black, ylims = [0.295, 0.52])
plot!(plot_f, timesteps, Normal_MeasPred_means, ribbon=(sqrt.(Normal_MeasPred_covariance), sqrt.(Normal_MeasPred_covariance)), lw = 2, lab = "Normal")
plot!(plot_f, timesteps, F1_MeasPred_means, ribbon=(sqrt.(F1_MeasPred_covariance), sqrt.(F1_MeasPred_covariance)), lw = 2, lab = "F1")
plot!(plot_f, timesteps, F2_MeasPred_means, ribbon=(sqrt.(F2_MeasPred_covariance), sqrt.(F2_MeasPred_covariance)), lw = 2, lab = "F2")
scatter!(plot_f, timesteps[2:end], Measurements'[2:end], lab = "Measurements")

#Plot 4:

plot_e = plot(GT_times, GT_val_matrix[1,:], lab = "Ground truth", xlab = "time (s)", ylab = "tank 1 level (m)",
                guidefontsize = 10, tickfontsize=8, lw=2, linestyle=:dash, color=:black)
         plot!(plot_e, timesteps, Normal_StateEstimates[1,:],
            ribbon=(Normal_x_std[1,:], Normal_x_std[1,:]), fillalpha=0.25, lw=2, lab = "Normal")
         plot!(plot_e, timesteps, F1_StateEstimates[1,:],
            ribbon=(F1_x_std[1,:], F1_x_std[1,:]), lw=2, lab = "F1")
         plot!(plot_e, timesteps, F2_StateEstimates[1,:],
            ribbon=(F2_x_std[1,:], F2_x_std[1,:]), fillalpha=0.25, lw=3, lab = "F2", linestyle=:dot)

##OLD RANDOM PLOTS...

#normal conditions
c = plot(timesteps, PosteriorBeliefs[1,:], lab="Normal",
         xlab="time (s)", ylab="probability", legend=:topleft,
         guidefontsize=10, tickfontsize=8, legendfontsize=8, lw=2)
plot!(c, timesteps, PosteriorBeliefs[2,:], lab="Fault mode 1", lw=2)
plot!(c, timesteps, PosteriorBeliefs[3,:], lab="Fault mode 2", lw=2)

d = plot(timesteps, Applied_Inputs[1,:].*10^5, lab="U1", ylims=[4.5, 10],
         xlab="time (s)", ylab= "input flow rate (10⁻⁵ m³/s)", legend=:topright,
         guidefontsize=10, tickfontsize=8, legendfontsize=8, lw=2)
plot!(d, timesteps, Applied_Inputs[2,:].*10^5, lab="U2", lw=3, linestyle = :dot)

# n = length(timesteps)
# step = 25                    # spacing — raise/lower to taste
# i1 = 1:step:n                # U1 marker indices
# i2 = (1 + step÷2):step:n     # U2 markers offset halfway → interleaved

# # lines only, no markers
# d = plot(timesteps, Applied_Inputs[1,:].*10^5, lab="U1", lw=3, ylims=[5,7],
#          xlab="time (s)", ylab="input flow rate (10⁻⁵ m³/s)",
#          legend=:topright, guidefontsize=16, tickfontsize=13, legendfontsize=11)
# plot!(d, timesteps, Applied_Inputs[2,:].*10^5, lab="U2", lw=3)

# # markers on subsampled points only
# scatter!(d, timesteps[i1], Applied_Inputs[1,i1].*10^5,
#          marker=:circle, ms=5, color=1, label="")
# scatter!(d, timesteps[i2], Applied_Inputs[2,i2].*10^5,
#          marker=:square, ms=5, color=2, label="")

final = plot(c, d,
             layout=(1,2), size=(1200, 390), dpi=300,
             left_margin=8Plots.mm, right_margin=5Plots.mm,
             top_margin=5Plots.mm, bottom_margin=8Plots.mm,
             titlefontsize=18)

plot_a = plot(GT_times, GT_val_matrix[3,:], lab="Ground truth", xlab="time (s)", ylab= "tank 3 level (m)",
            guidefontsize=10, tickfontsize=8, lw=2, linestyle=:dash, color=:black)
        plot!(plot_a, timesteps, Normal_StateEstimates[3,:],
            ribbon=(Normal_x_std[3,:], Normal_x_std[3,:]), lw=2, lab = "Normal")
        plot!(plot_a, timesteps, F1_StateEstimates[3,:],
            ribbon=(F1_x_std[3,:], F1_x_std[3,:]), lw=2, lab = "F1")
        plot!(plot_a, timesteps, F2_StateEstimates[3,:],
            ribbon=(F2_x_std[3,:], F2_x_std[3,:]), lw=2, lab = "F2")
scatter!(plot_a, timesteps[2:end], Measurements'[2:end], lab = "Measurement")

plot(timesteps, Real_Process[1,:], lab = "Ground truth")
plot!(timesteps, Normal_StateEstimates[1,:], ribbon = (Normal_x_std[1, :], Normal_x_std[1, :]), lab = "Normal")
plot!(timesteps, F1_StateEstimates[1,:], ribbon = (F1_x_std[1, :], F1_x_std[1, :]), lab = "F1")
plot!(timesteps, F2_StateEstimates[1,:], ribbon = (F2_x_std[1, :], F2_x_std[1, :]), lab = "F2")

