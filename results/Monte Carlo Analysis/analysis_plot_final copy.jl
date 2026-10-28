# Convergence heatmaps for u1, u2, u3 over prediction horizon (PH) x control horizon (CH).
# Pure-Julia version using Plots.jl + GR (no CairoMakie / no Python).
# Each cell = runs that converged correctly (< 126); grey = input absent for that config.
#
# One-time setup (run once in the REPL):  import Pkg; Pkg.add("Plots")
 
using Plots
using Plots.PlotMeasures          # brings mm into scope for margins
gr()
 
# rows = prediction horizon PH 1..5,  cols = control horizon CH 1..3
u1 = [ 94  NaN  NaN
       93  NaN  NaN
       94  NaN  NaN
       90  NaN  NaN
       92  NaN  NaN ]
 
u2 = [ NaN  NaN  NaN
       93   93   NaN
       77   92   NaN
       48   93   NaN
       36   92   NaN]
 
u3 = [ NaN  NaN  NaN
       NaN  NaN  NaN
       92   91   65
       81   80   94
       45   55   96 ]
 
panels = [("Control horizon = 1", u1), ("Control horizon = 2", u2), ("Control horizon = 3", u3)]
 
allv = filter(!isnan, vcat(vec(u1), vec(u2), vec(u3)))
vmin, vmax = minimum(allv), maximum(allv)      # shared scale so panels are comparable
grad = cgrad(:RdYlGn)
 
CH, PH = 1:3, 1:5                               # 3 columns, 5 rows per block
 
plts = Any[]
for (ttl, M) in panels
    p = heatmap(CH, PH, M;
        clims = (vmin, vmax), color = grad,
        yflip = true,                          # PH = 1 at the top, like a matrix
        title = ttl, titlefontsize = 11,
        xlabel = "Input sequence",
        ylabel = "Prediction horizon",
        xticks = 1:3, yticks = 1:5,
        aspect_ratio = 1, grid = false,
        colorbar = false,                      # no colour bar -> all three equal size
        background_color_inside = "#ededed")   # colour of empty (NaN) cells
 
    # white cell borders
    vline!(p, 0.5:1:3.5; color = :white, lw = 1.5, label = "")
    hline!(p, 0.5:1:5.5; color = :white, lw = 1.5, label = "")
 
    # print the value inside each populated cell, text colour chosen for contrast
    for i in 1:5, j in 1:3
        v = M[i, j]
        isnan(v) && continue
        t = (v - vmin) / (vmax - vmin)          # 0 = worst (red), 1 = best (green)
        annotate!(p, j, i, text(string(round(Int, v)), t < 0.32 ? :white : :black, :center, 10))
    end
    push!(plts, p)
end
 
fig = plot(plts...; layout = (1, 3), size = (1000, 380),
    left_margin = 4mm, bottom_margin = 5mm, top_margin = 3mm, right_margin = 4mm)