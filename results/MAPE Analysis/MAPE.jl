using JLD2
using StatsPlots

## Load MAPE vectors (119 states each)

N_mape  = vec(load(joinpath(@__DIR__, "N_mape.jld2"), "N_mape"))
f1_mape = vec(load(joinpath(@__DIR__, "f1_mape.jld2"), "f1_mape"))  #saved under the name "f2_mape" in f1_mape.jld2
f2_mape = vec(load(joinpath(@__DIR__, "f2_mape.jld2"), "f2_mape"))


N_mape_no_integral_states = vcat(N_mape[1:103], N_mape[116:118])
f1_mape_no_integral_states = vcat(f1_mape[1:103], f1_mape[116:118])
f2_mape_no_integral_states = vcat(f2_mape[1:103], f2_mape[116:118])

## Box plot of all 3 models

ytick_vals = [0.01, 0.1, 1, 10, 100, 1000]

boxplot(["N" "F1" "F2"], [N_mape f1_mape f2_mape],
        legend = false, yscale = :log10,
        yticks = (ytick_vals, ["0.01", "0.1", "1", "10", "100", "1000"]),
        xlab = "model", ylab = "MAPE (%)")

plot(N_mape, ylab = "MAPE (%)", xlab = "state", lab = :none)
## Box plot of all 3 models

ytick_vals = [0.01, 0.1, 1, 10, 100, 1000]

boxplot(["N" "F1" "F2"], [N_mape_no_integral_states f1_mape_no_integral_states f2_mape_no_integral_states],
        legend = false, yscale = :log10,
        yticks = (ytick_vals, ["0.01", "0.1", "1", "10", "100", "1000"]),
        xlab = "model", ylab = "MAPE (%)")

plot(N_mape, ylab = "MAPE (%)", xlab = "state", lab = :none)