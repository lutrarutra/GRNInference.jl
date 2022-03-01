using DataFrames, CSV

function to_vec(m::Matrix{Float64})
    v::Vector{Vector{Float64}} = []

    for i in 1:size(m)[1]
        push!(v, m[i,:])
    end
    v
end


struct Observations
    t::Vector{Float64}
    u::Vector{Vector{Float64}}
    labels::Vector{String}      

    function Observations(filename::String)
        data = DataFrame(CSV.File(filename))
        new(Float64.(data[:,1]), to_vec(Matrix(data[:, 2:6])), names(data)[2:end])
    end
end







