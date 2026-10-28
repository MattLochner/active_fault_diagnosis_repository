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
    A = ForwardDiff.jacobian(f, xhatmin1)
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
        H = ForwardDiff.jacobian(h, xhatneg)
        yhat = measurement_pred_function(xhatneg, params, tend)
    #------------------------------------------------------
        Py = H*Pneg*H'+ R #assumption simplififying MRM to R

    return yhat, Py
end 

struct InputDesignParams{TR, TQ, PP, MF, OF}
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
    propparams         :: PP
    measurement_pred_function :: MF
    model_ODE_function        :: OF
    no_of_modes        :: Int
    beliefs            :: Vector{Float64}
end

function Input_Design_Function(u_vec, p)
    Hc, Hp = p.control_horizon, p.prediction_horizon
    U = reshape(u_vec, 1, Hc)         
    #U = reshape(u_vec, 1, 2*Hc)

    #yhats = Vector{Vector{Float64}}(undef, p.no_of_modes) #store a vector of vectors of vectors i.e. each timestep should have its own vector length of prediction horizon
    #Pys   = Vector{Matrix{Float64}}(undef, p.no_of_modes) #store a vector of vector of matrices, length of prediction horizon
    yhats = [Vector{Vector{Float64}}(undef, p.no_of_modes) for _ in 1:Hp]  # [t][mode] -> mean
    Pys   = [Vector{Matrix{Float64}}(undef, p.no_of_modes) for _ in 1:Hp]  # [t][mode] -> cov
    params = deepcopy(p.propparams)               #

    #Normal Operating conditions:
        params.faultmode = 0            # mode index 1,2,3  ->  faultmode 0,1,2
        m = 1
        x = copy(p.xhat0_modes[m])           # start from THIS mode's posterior
        P = copy(p.P0_modes[m])

        for j in 1:Hp
            col = j ≤ Hc ? j : Hc            # change over control horizon, then hold
            params.control.setpoint[8]= U[1, col]
            #params.control.setpoint[9]=U[1,col+Hc]

            tj  = p.t_start + (j-1)*p.N_meas
            tj1 = p.t_start +  j   *p.N_meas
            x, P = ExtendedStatePrediction(x, P, p.nx, p.ny, p.R, p.Q[m],
                                           params, tj, tj1,
                                           p.measurement_pred_function,
                                           p.model_ODE_function)
            yhats[j][m], Pys[j][m] = LinearisedMeasurementPrediction(
                                         x, P, params, p.R,
                                         p.measurement_pred_function, tj1)           
        end

        t_end = p.t_start + Hp*p.N_meas
        # yhats[m], Pys[m] = LinearisedMeasurementPrediction(
        #                        x, P, params, p.R, p.measurement_pred_function, t_end)

    #Fault mode 1
        params.faultmode = 8            # mode index 1,2,3  ->  faultmode 0,1,2
        m = 2
        x = copy(p.xhat0_modes[m])           # start from THIS mode's posterior
        P = copy(p.P0_modes[m])

        for j in 1:Hp
            col = j ≤ Hc ? j : Hc           # change over control horizon, then hold
            params.control.setpoint[8]= U[1, col]
            #params.control.setpoint[9]=U[1,col+Hc]

            tj  = p.t_start + (j-1)*p.N_meas
            tj1 = p.t_start +  j   *p.N_meas
            x, P = ExtendedStatePrediction(x, P, p.nx, p.ny, p.R, p.Q[m],
                                           params, tj, tj1,
                                           p.measurement_pred_function,
                                           p.model_ODE_function)
            yhats[j][m], Pys[j][m] = LinearisedMeasurementPrediction(
                                         x, P, params, p.R,
                                         p.measurement_pred_function, tj1)        
        end

        t_end = p.t_start + Hp*p.N_meas
        # yhats[m], Pys[m] = LinearisedMeasurementPrediction(
        #                        x, P, params, p.R, p.measurement_pred_function, t_end)
    #Faultmode 2
        params.faultmode = 5            # mode index 1,2,3  ->  faultmode 0,1,2
        m = 3
        x = copy(p.xhat0_modes[m])           # start from THIS mode's posterior
        P = copy(p.P0_modes[m])

        for j in 1:Hp
            col = j ≤ Hc ? j : Hc            # change over control horizon, then hold
            params.control.setpoint[8]= U[1, col]
            #params.control.setpoint[9]=U[1,col+Hc]
            
            tj  = p.t_start + (j-1)*p.N_meas
            tj1 = p.t_start +  j   *p.N_meas
            x, P = ExtendedStatePrediction(x, P, p.nx, p.ny, p.R, p.Q[m],
                                           params, tj, tj1,
                                           p.measurement_pred_function,
                                           p.model_ODE_function)
            yhats[j][m], Pys[j][m] = LinearisedMeasurementPrediction(
                                         x, P, params, p.R,
                                         p.measurement_pred_function, tj1)            
        end

        t_end = p.t_start + Hp*p.N_meas
        # yhats[m], Pys[m] = LinearisedMeasurementPrediction(
        #                        x, P, params, p.R, p.measurement_pred_function, t_end)
    
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