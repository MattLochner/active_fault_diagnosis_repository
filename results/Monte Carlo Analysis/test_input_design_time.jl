using JLD2
using Plots

all_results = load("Input_design_MonteCarlo_141150.jld2", "all_results")

# for i in 1:12
#     display(all_results[1][i])
# end

save_simulation = all_results[1]

using XLSX

"""
    save_mc_results_wide(filename, seeds, all_results;
                         sheetname="all results", start_row=3)

Wide format: one row per random seed, no header row. First data row
lands at `start_row`. Opens an existing workbook if present
(preserving other sheets), otherwise creates a new one.
"""
function save_mc_results_wide(filename::AbstractString,
                              seeds::AbstractVector,
                              all_results::AbstractVector;
                              sheetname::AbstractString="all results",
                              start_row::Integer=3)

    @assert length(seeds) == length(all_results)
    isempty(seeds) && error("no results to write")

    ncols = 1  # random_seed, then blocks
    for blk in all_results[1]
        ncols += 6 + length(blk.u_opt)      # was 5; +1 for cpu_time
    end

    rows = Vector{Vector{Any}}(undef, length(seeds))
    for (k, res) in enumerate(all_results)
        row = Any[seeds[k]]
        for blk in res
            push!(row, blk.control_horizon, blk.prediction_horizon,
                       blk.iterations, blk.f_calls, blk.seconds, blk.cpu_time)  # cpu after seconds
            append!(row, blk.u_opt)
        end
        length(row) == ncols ||
            error("seed $(seeds[k]) layout mismatch: $(length(row)) vs $ncols cols")
        rows[k] = row
    end

    mode = isfile(filename) ? "rw" : "w"
    XLSX.openxlsx(filename, mode=mode) do xf
        sheet = if sheetname in XLSX.sheetnames(xf)
            xf[sheetname]
        elseif mode == "w"
            s = xf[1]
            XLSX.rename!(s, sheetname)
            s
        else
            XLSX.addsheet!(xf, sheetname)
        end

        for r in 1:length(rows), c in 1:ncols
            sheet[start_row + r - 1, c] = rows[r][c]
        end
    end
    return filename
end

start = 141
final = 150
for i in start:final
    save_simulation = all_results[i-start+1]
    save_mc_results_wide("Book1.xlsx", [1002], [save_simulation]; sheetname="All runs", start_row = i+1)
end