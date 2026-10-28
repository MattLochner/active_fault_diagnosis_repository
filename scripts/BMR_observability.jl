using full_simulation
#using ForwardDiff
using FiniteDiff
using LinearAlgebra
using ForwardDiff

include("BMR_kf_params.jl")
using .Kf_model_params

kfparams = Kf_data
kfparams.faultmode = 0

#NOTE: Fault modes are numbered according to Miskin (2016) master's thesis

#0 - normal operating conditions
#8 - fault mode 1 - level measurement positive bias fault
#5 - fault mode 2 - heat exchanger fouling fault

include("BMR_initial_conditions.jl")

xhat0 = initial_state_estimate  #intial_state estimate for kalman filter
xhat1 = xhat0*1.01

t_start = 0.0
#dx_dt = calculate_model_autoclave_odes(xhat0, Kf_data, 0.0)

f(x) = calculate_model_autoclave_odes(x, kfparams, 0.0)
h(x) = BMR_kf_measurement_prediction_with_samples(x, kfparams, 0.0)
#h(x) = BMR_kf_measurement_prediction(x, kfparams, 0.0)

f(xhat0)
j = h(xhat0)

for i in 1:36
    display(j[i])
end
# A matrix: Jacobian of f at operating point
A1 = FiniteDiff.finite_difference_jacobian(f, xhat1)
A = ForwardDiff.jacobian(f, xhat1)

eigs = eigen(A)
e = eigs.values

for i in 1:118
    display(e[i])
end

# C matrix: Jacobian of measurement function at operating point
C1 = FiniteDiff.finite_difference_jacobian(h, xhat1)
C = ForwardDiff.jacobian(h, xhat1)

println("Computing A (130×130)...")
A = ForwardDiff.jacobian(f, xhat1)
println("Computing C (37×130)...")
C = ForwardDiff.jacobian(h, xhat1)

# ── Build observability matrix ────────────────────────────────────────────────
function build_observability(C, A, n_blocks)
    p, n = size(C)
    blocks = Vector{Matrix{Float64}}(undef, n_blocks)
    CAi = copy(C)
    for i in 1:n_blocks
        blocks[i] = CAi
        CAi = CAi * A
    end
    return vcat(blocks...)
end

function check_observability(O, n, label, tol)
    r = rank(O, tol)
    println("\n── $label ────────────────────────────────────────────────────")
    println("  O size:     $(size(O, 1)) × $(size(O, 2))")
    println("  Rank:       $r / $n")
    println("  Observable: $(r == n ? "✓ YES" : "✗ NO — $(n - r) unobservable mode(s)")")
    return r
end

n = size(A, 1)   # 130
p = size(C, 1)   # 37
tol = 1e-6 

O_p = build_observability(C, A, p)          # 37 blocks → 1369 × 130
O_n = build_observability(C, A, n)          # 130 blocks → 4810 × 130

rank(O_p)

dt = 0.2           #length of timestep
Ad = exp(A * dt)   # matrix exponential, NOT element-wise (.*)

# Sanity check — spectral radius must be < 1
ρ = maximum(abs.(eigvals(Ad)))
println("Spectral radius of Ad: $ρ")   # must be < 1 for stable system

# Rebuild observability with Ad
O_p2 = build_observability(C, Ad, p)
O_n2 = build_observability(C, Ad, n)

# Use default tolerance (no manual tol — Julia scales it to O's own singular values)
rank_p = rank(O_p2)
println("rank_p = $rank_p")
rank_n = rank(O_n2)
println("rank_n = $rank_n")
println("\n── Comparison ───────────────────────────────────────────────────")
if rank_p == rank_n
    println("  ✓ Ranks agree — blocks beyond $p added no information")
else
    println("  ! Ranks differ — extra modes found between block $p and $n: $(rank_n - rank_p)")
end

println("\n── SVD analysis of O_n ──────────────────────────────────────────")
U, S, V = svd(O_n2)
println("  Largest 5 singular values:  ", join(round.(S[1:min(5,end)], sigdigits=4), "  "))
println("  Smallest 10 singular values:")
for i in length(S):-1:max(1, length(S)-9)
    println("    σ$i = $(round(S[i], sigdigits=4))")
end

if rank_n < n
    n_unobs = n - rank_n
    println("\n── Unobservable subspace ($n_unobs mode(s)) ──────────────────────────")

    # Last n_unobs columns of V span the unobservable subspace
    V_unobs = V[:, rank_n+1:end]

    # Per-state score: norm of its row in V_unobs (1.0 = fully unobservable, 0.0 = fully observable)
    state_unobs_score = [norm(V_unobs[i, :]) for i in 1:n]
    sorted_idx = sortperm(state_unobs_score, rev=true)

    println("  State  Score   (1.0 = fully unobservable, 0.0 = fully observable)")
    for i in 1:n
        idx = sorted_idx[i]
        score = state_unobs_score[idx]
        score < 0.01 && break
        println("  X[$idx]$(lpad("", max(0, 5-ndigits(idx))))  $(round(score, digits=4))")
    end
else
    println("\n  ✓ System is fully observable.")
end