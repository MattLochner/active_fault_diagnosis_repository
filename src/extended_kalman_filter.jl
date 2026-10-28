#using ForwardDiff
using FiniteDiff

function ekf(xhatmin1, Pposmin1, nx, ny, R, Q, params, tstart, tend, measurement_pred_function, model_ODE_function, y)

    #calculate prior state estimate mean
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
    A = FiniteDiff.finite_difference_jacobian(f, xhatmin1)
    #----------------------------------------------------------------------

    #calculate the prior covariance matrix

    #Option 1: ODE solver

    # prob = ODEProblem(lyap!, Pposmin1, (tstart, tend), (A, Q))
    # sol = solve(prob, Tsit5(), saveat = tend)

    # Pneg_1 = sol.u[2]

    # Pneg_1 = (Pneg_1 + Pneg_1') / 2

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

    # diagnose + symmetrize
    Pneg_2, λmin2, asym2 = check_covariance(Pneg_2; name="Pneg_2")
    #Pneg_1, λmin1, asym1 = check_covariance(Pneg_1; name="Pneg_1")

    # only assert true positive-definiteness on the symmetrized version
    #@assert isposdef(Symmetric(Pneg_1)) "Pneg_1 not PD even after symmetrizing (λmin=$λmin1)"
    @assert isposdef(Symmetric(Pneg_2)) "Pneg_2 not PD even after symmetrizing (λmin=$λmin2)"
    
    #check1 = cholesky(Symmetric(Pneg_1)) #check if matrix is positive definite
    #check2 = cholesky(Symmetric(Pneg_2)) #check if matrix is positive definite

    #Measurement prediction
    #------------------------------------------------------
        h(x) = measurement_pred_function(x, params, tend)
        H = FiniteDiff.finite_difference_jacobian(h, xhatneg)
        yhat = measurement_pred_function(xhatneg, params, tend)
    #------------------------------------------------------
        # H = ForwardDiff.jacobian(measurement_pred_function, xhatneg)
        # yhat = measurement_pred_function(xhatneg)
    #------------------------------------------------------
    Py = H*Pneg_2*H'+ R #assumption simplififying MRM to R

    #State estimation

    K = Pneg_2*H'*inv(Py)    #kalman gain matrix
    #K = (Py' / (H * Pneg_2'))   #better for numerical stability 

    xhatpos = xhatneg + K*(y .- yhat)

    Ppos = (I - K*H) * Pneg_2 * (I - K*H)' + K*R*K' #assumption simplifying KMRMK to KRK (i.e. M is identity matrix)
  
    return yhat, Py, xhatpos, Ppos         #it must just return the prior state estimate and posterior state estimate so that Bayesian update of probability can work
end


function run_ODE_model(tStart, tEnd, x, p, model_ODE_function)
        prob3 = ODEProblem(model_ODE_function, x, (tStart, tEnd), p)
        ODEsol2 = solve(prob3, Tsit5(), saveat = tEnd)
        return ODEsol2.u[2]
end

function lyap!(dP, P, p, t)
    A, Q = p
    mul!(dP, A, P)              # dP = A*P
    mul!(dP, P, A', 1.0, 1.0)  # dP += P*A'
    dP .+= Q                   # dP += Q
end

function check_covariance(P; name="P", tol=1e-10)
    # 1. Asymmetry magnitude
    asym = maximum(abs.(P - P'))

    # 2. Symmetrize, then look at the true (real) eigenvalues
    Psym = (P + P') / 2
    λ = eigvals(Symmetric(Psym))
    λmin = minimum(λ)

    if λmin < -tol
        @warn "$name has a large negative eigenvalue" λmin asymmetry=asym
    end
    if λmin < 0
        @info "$name has a negative eigenvalue — adding 1e-10 + min(λ) to diagonal" λmin asymmetry=asym
        jitter = abs(minimum(λ))
        Psym += 1e-10 * I + jitter * I

        # Optional: recompute λmin after repair
        λmin = minimum(eigvals(Symmetric(Psym)))
    end

    return Psym, λmin, asym
end