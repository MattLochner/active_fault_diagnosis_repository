using Plots
using Plots.PlotMeasures

Nc  = [1, 1, 1, 1, 1,  2, 2, 2, 2, 3, 3, 3]
Np  = [1, 2, 3, 4, 5,  2, 3, 4, 5, 3, 4, 5]

mean_ = [98.9, 203.4, 364.5, 459.2, 723.8, 319.2, 1585.5, 956.1, 915.5, 817.1, 2416.5, 1607.7]
std_  = [157.7, 249.9, 518.0, 648.9, 1050.1, 281.2, 5360.5, 3021.4, 746.0, 601.4, 7264.7, 4519.0]
mn    = [7.5, 23.4, 28.7, 52.1, 76.9, 68.6, 71.7, 157.7, 248.9, 249.6, 327.4, 332.8]
q1    = [32.2, 71.8, 143.5, 217.0, 289.9, 161.5, 283.1, 389.0, 532.8, 397.5, 591.7, 667.2]
med   = [55.8, 126.5, 209.9, 299.4, 396.5, 216.1, 363.4, 517.6, 680.9, 604.6, 739.4, 841.7]
q3    = [86.9, 209.0, 324.9, 413.9, 605.9, 337.2, 511.1, 757.0, 918.4, 962.9, 1344.7, 1261.9]
mx    = [1046.9, 1552.8, 3387.1, 4926.0, 7000.7, 1672.8, 28779.8, 30582.1, 4601.0, 3276.9, 45868.1, 45326.4]
# ------------------------------------------------------------------
# Lay out x-positions: group by prediction horizon (Np),
# control horizons (Nc) sit next to each other inside each group
# ------------------------------------------------------------------
groups   = sort(unique(Np))
gap      = 0.9
boxwidth = 0.6
palette  = Dict(1 => :blue, 2 => :maroon, 3 => :darkorange)   # colour = Nc

xpos        = zeros(length(Np))
groupcenter = Dict{Int,Float64}()
x = 0.0
for g in groups
    idx = findall(==(g), Np)
    idx = idx[sortperm(Nc[idx])]     # order by control horizon within the group
    xs  = Float64[]
    for i in idx
        x += 1.0
        xpos[i] = x
        push!(xs, x)
    end
    groupcenter[g] = sum(xs) / length(xs)
    x += gap
end

# x-ticks: one tick per Np group, placed at the group centre
tickpos = [groupcenter[g] for g in groups]
ticklab = string.(groups)

# ------------------------------------------------------------------
# Base plot
# ------------------------------------------------------------------
plt = plot(;
    yscale     = :log10,
    ylims      = (1.0, maximum(mx) * 3),
    yticks     = ([10, 100, 720, 1000, 10000], ["10", "100", "720" , "1000", "10000"]),
    xticks     = (tickpos, ticklab),
    xlabel     = "Prediction horizon (Np)",
    ylabel     = "Optimisation solve CPU-time [s]  (log scale)",
    legend       = :topleft,
    legendtitle  = "Control horizon (Nc)",          # <-- legend heading
    legendtitlefontsize = 9,                    # optional: match tick size
    grid       = :y,
    framestyle = :box,
    size       = (950, 560),
    guidefontsize = 10, tickfontsize = 9, titlefontsize = 12,
    left_margin   = 8mm,
    bottom_margin = 8mm,
    top_margin    = 6mm,
    right_margin  = 4mm,
)

# ------------------------------------------------------------------
# Draw each box + whiskers; label one box per Nc for the legend
# ------------------------------------------------------------------
labeled = Set{Int}()
for i in eachindex(Nc)
    xc = xpos[i]; hw = boxwidth / 2; c = palette[Nc[i]]

    lab = ""
    if !(Nc[i] in labeled)
        lab = "Nc = $(Nc[i])"
        push!(labeled, Nc[i])
    end

    plot!(Shape([xc-hw, xc+hw, xc+hw, xc-hw], [q1[i], q1[i], q3[i], q3[i]]);
          fillcolor = c, fillalpha = 0.65, linecolor = :black, lw = 1, label = lab)
    plot!([xc-hw, xc+hw], [med[i], med[i]]; linecolor = :black, lw = 2, label = "")
    plot!([xc, xc], [q3[i], mx[i]]; linecolor = :black, lw = 1, label = "")
    plot!([xc, xc], [mn[i], q1[i]]; linecolor = :black, lw = 1, label = "")
    plot!([xc-hw/2, xc+hw/2], [mx[i], mx[i]]; linecolor = :black, lw = 1, label = "")
    plot!([xc-hw/2, xc+hw/2], [mn[i], mn[i]]; linecolor = :black, lw = 1, label = "")
end

hline!([720]; linecolor = :red, linestyle = :dash, lw = 1.5, label = "measurement interval")
#annotate!(maximum(xpos) + 0.4, 760 * 1.15, text("760 s", 8, :red, :right))

#savefig(plt, "solve_time_boxplot.png")
display(plt)