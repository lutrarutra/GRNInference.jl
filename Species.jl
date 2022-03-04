abstract type Species end

mutable struct Gene <: Species
    name::String
    idx::Int64
    u₀::Float64
    targets::Vector{Species}

    num_params::Int64
    param_idx::Vector{Int64}

    Gene(name, idx::Int64) = new(name, idx, 0.05, [], 0, [])
end

idx(g::Gene)::Int64 = g.idx
params(g::Gene, p::Vector{Float64})::Vector{Float64} = p[g.param_idx]
state(g::Gene, u::Vector{Float64})::Float64 = u[g.idx]
# targets_u(g::Gene, u::Vector{Float64})::Vector{Float64} = u[idx.(g.targets)]
# targets_p(g::Gene, p::Vector{Float64})::Vector{Float64} = u[idx.(g.targets)]

function d_dt!(g::Gene, du::Vector{Float64}, u::Vector{Float64}, p::Vector{Float64})
    du[g.idx] = (
        sum(map(gp -> u[gp.idx] * p[gp.param_idx][2], g.targets), init=0.0)
        - u[g.idx] * sum(map(gp -> u[gp.protein.idx] * p[gp.param_idx][1], g.targets), init=0.0)
    )
    nothing
end

mutable struct mRNA <: Species
    name::String
    idx::Int64
    u₀::Float64
    transcriptor::Species
    
    # params
    # transcription_rate::Float64
    # degradation_rate::Float64
    num_params::Int64
    param_idx::Vector{Int64}

    mRNA(name, idx::Int64, transcriptor::Species) =
        new(name, idx, 0.0, transcriptor, 2, [])

end

idx(m::mRNA)::Int64 = m.idx
params(m::mRNA, p::Vector{Float64})::Vector{Float64} = p[m.param_idx]
state(m::mRNA, u::Vector{Float64})::Float64 = u[m.idx]

function d_dt!(m::mRNA, du::Vector{Float64}, u::Vector{Float64}, p::Vector{Float64})
    du[m.idx] = (
        u[m.transcriptor.idx] * p[m.param_idx][1]   # transcription_rate
        - u[m.idx] * p[m.param_idx][2]              # degradation_rate
    )
    nothing
end

mutable struct Protein <: Species
    name::String
    idx::Int64
    u₀::Float64
    targets::Vector{Species}
    translator::mRNA
    
    # params
    # translation_rate::Float64
    num_params::Int64
    param_idx::Vector{Int64}

    Protein(name::String, idx::Int64, translator::mRNA) =
        new(name, idx, 0.0, [], translator, 1, [])
end

idx(protein::Protein)::Int64 = protein.idx
params(protein::Protein, p::Vector{Float64})::Vector{Float64} = p[protein.param_idx]
state(protein::Protein, u::Vector{Float64})::Float64 = u[protein.idx]

function d_dt!(protein::Protein, du::Vector{Float64}, u::Vector{Float64}, p::Vector{Float64})
    du[protein.idx] = (
        u[protein.translator.idx] * p[protein.param_idx][1]
        - u[protein.idx] * sum(map(target -> u[target.idx] * p[target.param_idx][1], protein.targets), init=0.0)
        + sum(map(target -> u[target.idx] * p[target.param_idx][2], protein.targets), init=0.0)
    )
    nothing
end


mutable struct Complex <: Species
    name::String
    idx::Int64
    u₀::Float64
    species1::Species
    species2::Species

    # params
    # binding_rate::Float64
    # release_rate::Float64
    num_params::Int64
    param_idx::Vector{Int64}
    

    function Complex(name, idx::Int64, s1::Species, s2::Species)
        c = new(name, idx, 0.0, s1, s2, 2, [])
        push!(s1.targets, c)
        push!(s2.targets, c)
        c
    end
end

idx(c::Complex)::Int64 = c.idx
sub_idx(c::Complex)::Vector{Int64} = [c.species1.idx, c.species2.idx]

function d_dt!(c::Complex, du::Vector{Float64}, u::Vector{Float64}, p::Vector{Float64})
    du[c.idx] = (
        u[c.protein1.idx] * u[c.protein2.idx] * p[c.param_idx][1]   # binding_rate
        - u[c.idx] * p[c.param_idx][2]                              # release_rate          
    )
    nothing
end


