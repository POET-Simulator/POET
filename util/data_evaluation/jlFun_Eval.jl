using RData, DataFrames, Plots

# Load all .rds files in given directory 
function load_data(directory)
    files = readdir(directory)
    # data is dictionary with iteration number as key
    data = Dict{Int,Any}()
    for file in files
        if (endswith(file, ".rds") && startswith(file, "iter"))
            # extract iteration number (iter_XXX.rds)
            iter = parse(Int, split(split(file, "_")[2], ".")[1])
            data[iter] = load(joinpath(directory, file))
        end
    end
    return data
end

function spec_to_mat(in_df::DataFrame, spec::Symbol, cols::Int)
    specvec = in_df[!, spec]
    specmat = transpose(reshape(specvec, cols, :))

    return specmat
end

function plot_field(data::AbstractMatrix, log::Bool=false)
    if log
        heatmap(1:size(data, 2), 1:size(data, 1), log10.(data), c=:viridis)
    else
        heatmap(1:size(data, 2), 1:size(data, 1), data, c=:viridis)
    end
end