# Convergence heatmaps for u1, u2, u3 over prediction horizon (PH) x control horizon (CH).
# Pure-Julia version using Plots.jl + GR (no CairoMakie / no Python).
# Each cell = runs that converged correctly (< 126); grey = input absent for that config.
#
# One-time setup (run once in the REPL):  import Pkg; Pkg.add("Plots")
 
using Plots
using Plots.PlotMeasures

gr()
 
# rows = prediction horizon PH 1..5,  cols = control horizon CH 1..3
T = [ 2   NaN  NaN
       6   8    NaN
       10  15   40
       10  25   48
       22  40   61 ]
 
vals = filter(!isnan, vec(T))
vmin, vmax = minimum(vals), maximum(vals)
grad = cgrad(:RdYlGn, rev=true)             # green (low) → yellow → red (high)
nrow, ncol = size(T)

p = heatmap(1:ncol, 1:nrow, T;
    clims = (vmin, vmax), color = grad,
    yflip = true,
    #title = "Temperature", titlefontsize = 11,
    xlabel = "Hc", ylabel = "Hp",
    xticks = 1:ncol, yticks = 1:nrow,
    grid = false,                           # aspect_ratio removed → fills the panel
    colorbar = false, #, colorbar_title = "Temperature (°C)",
    left_margin = 5mm, top_margin = 3mm,
    bottom_margin = 7mm, right_margin = 12mm,
    background_color_inside = "#ededed")

# white cell borders
vline!(p, 0.5:1:(ncol + 0.5); color = :white, lw = 1.5, label = "")
hline!(p, 0.5:1:(nrow + 0.5); color = :white, lw = 1.5, label = "")

# print each value, contrast-aware text colour
for i in 1:nrow, j in 1:ncol
    v = T[i, j]
    isnan(v) && continue
    t = (v - vmin) / (vmax - vmin)
    annotate!(p, j, i, text(string(round(Int, v)), t < 0.32 ? :white : :black, :center, 10))
end

display(p)