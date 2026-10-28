function multiple_model_method(y_meas, priors, yhats, cov_matrices)
    # Bayesian update of mode probabilities from a single measurement vector.
    # Inputs (d = no. of measured outputs, n = no. of modes):
    #   y_meas       - Vector of length d: the measured outputs at this time step
    #   priors       - Vector of length n: prior probability of each mode (all > 0, should sum to 1)
    #   The measurement predictions are summarised by the matrices, yhats and cov_matrices
    #   yhats        - n×d Matrix: row i is mode i's predicted output (yhats[i,:] must be length d)
    #   cov_matrices - Vector of n d×d matrices: innovation covariance S for each mode
    #                  (must be symmetric positive definite, else cholesky fails)
    # Returns: Vector of length n, the normalised posterior mode probabilities.


    d       = length(y_meas)
    n_modes = length(priors)
    log_p   = zeros(n_modes)

    for i in 1:n_modes
        r        = y_meas - yhats[i,:]
        S        = Symmetric(cov_matrices[i])
        Sc       = cholesky(S)
        mahal2   = dot(r, Sc \ r)

        #Sc \ r          # [0.4545..., 0.1818...]
        #inv(S) * r      # [0.4545..., 0.1818...]   ← iden, i'm just leaving this here as a note. Sc \ R is simply a more numerically efficent way of doing inv(S) * r

        logdetS  = 2 * sum(log, diag(Sc.U))
        log_p[i] = -0.5 * (d*log(2π) + logdetS + mahal2)
    end

    # Combine likelihood and prior in log space:
    log_unnorm = log_p .+ log.(priors)     # log(p(y|xᵢ) · prior_i)

    # Shift so the max is zero, then exponentiate and normalize:
    m = maximum(log_unnorm)
    w = exp.(log_unnorm .- m)                    # unnormalized posterior, max element = 1.0

    posterior = w ./ sum(w)                           # normalize
    return posterior
end

function flatten_prior(vector, factor) #power method -> it shrinks all probabilities towards the uniform
    factored_vector = vector.^factor
    normalisation_constant = sum(factored_vector)
    flattened_posterior = factored_vector./normalisation_constant
    return flattened_posterior 
end

function linearflattening(prior, alpha)
    flattened_prior = (1 - alpha).* prior .+ alpha ./ length(prior)
    return flattened_prior
end

# A = [0.9, 0.05, 0.05]
# B = flatten_posterior(A, 0.9)
#sum(B)
# #testing out fault mode likelihood calculation:

# yhats = vcat([10, 5]',[12,4]')
# yhats[1,:]
# cov_matrices = [[1.0 0.2; 0.2 1],[1.0 0.2; 0.2 1]]

# y_meas = [10.1, 5]

# posterior = UpdateBeliefsMultipleMeasurements2(y_meas, [0.5, 0.5], yhats, cov_matrices)

