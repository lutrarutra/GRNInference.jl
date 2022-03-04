using Optim, DifferentialEquations, DiffEqSensitivity, LinearAlgebra

include("Species.jl")
include("Observations.jl")


mutable struct BioNetwork
    ix::Vector{Int64}
    species::Vector{Species}
    num_params::Int64
    obs::Observations
    ode_prob

    BioNetwork() = new([], [], 0)
end

Base.getindex(net::BioNetwork, idx::Int64)::Species = net.species[idx]
Base.getindex(net::BioNetwork, name::String)::Species = net.species[findfirst(s -> lowercase(s.name) == lowercase(name), net.species)]


function add_species!(net::BioNetwork, ttt::DataType, name::String; kwargs...)
    @assert ttt <: Species

    if ttt <: Gene
        push!(net.species, Gene(name, length(net.species)+1))
    elseif ttt <: mRNA
        push!(net.species, mRNA(name, length(net.species)+1, kwargs[:transcriptor]))
    elseif ttt <: Protein
        push!(net.species, Protein(name, length(net.species)+1, kwargs[:translator]))
    elseif ttt <: Complex
        push!(net.species, Complex(name, length(net.species)+1, kwargs[:species1], kwargs[:species2]))
    end

    append!(net.species[end].param_idx, collect(net.num_params+1:net.num_params+net.species[end].num_params))
    net.num_params += net.species[end].num_params
end

function init!(net::BioNetwork, obs::Observations)
    net.obs = obs
    net.ix = id_labels(net)
    for (i, k) in enumerate(net.ix)
        net.species[k].u₀ = net.obs.u[1,i]
    end
    nothing
end


function id_labels(net::BioNetwork)::Vector{Int64}
    ix::Vector{Int64} = []
    for (i,s) in enumerate(net.species)
        for label in net.obs.labels
            if lowercase(s.name) == lowercase(label)
                push!(ix, i)
                break
            end
        end
    end
    ix
end

function labels(net::BioNetwork)::Matrix{String}
    ls::Vector{String} = []
    for i in 1:length(net.species)
        push!(ls, net.species[i].name)
    end
    reshape(ls, 1, length(ls))
end

function d_dt!(net::BioNetwork, du, u, p)
    d_dt!.(net.species, repeat([du], length(net.species)), repeat([u], length(net.species)), repeat([p], length(net.species)))
    nothing
end


function loss(net::BioNetwork, θ)::Float64
    f!(du,u,p,t) = d_dt!(net,du,u,p)
    prob = ODEProblem(f!, state₀(net), (0.0, maximum(net.obs.t)), θ)
    sol = solve(prob)
    u = sol.(net.obs.t)

    L = zeros(Float64, length(net.ix))
    for (t, u_t) in enumerate(u)
        for (i, k) in enumerate(net.ix)
            L[i] += (u_t[k] - net.obs.u[t, i])^2
        end
    end
    sum(L)
end

function state₀(net::BioNetwork)::Vector{Float64}
    u = zeros(Float64, length(net.species))
    for i in 1:length(net.species)
        u[i] = net.species[i].u₀
    end
    u
end


function ∂u(net::BioNetwork, θ)
    N::Int64 = length(θ)
    M::Int64 = size(net.obs.u)[2]

    # _prob = remake(net.ode_prob, u0=net.ode_prob.u0, p=θ)
    u0 = copy(state₀(net))
    f!(du, u, p, t) = d_dt!(net, du, u, p)
    _prob = ODEForwardSensitivityProblem(f!, eltype(θ).(u0), (0.0, maximum(obs.t)), θ)
    sol = solve(_prob, saveat=net.obs.t[2:end])
    x, dp = extract_local_sensitivities(sol)
    
    d = zeros(Float64, (M, length(net.obs.t)-1))

    for (i, k) in enumerate(net.ix)
        d[i,:] = net.obs.u'[i,2:end] - x[k,:]
    end
    
    GD = zeros(Float64, N)
    
    for i in 1:N
        GD[i] = -2.0*dot(dp[i][net.ix,:], d)
    end

    GD
end



# function d_dt!(net::BioNetwork, du, u, p)
#     for (i,s) in enumerate(net.species)
#         du[i] = 0.0
#         if s isa Gene
#             for target in s.targets
#                 du[i] -= u[target.species1.idx] * u[target.species2.idx] * p[target.param_idx][1]  # binding_rate
#                 du[i] += u[target.idx] * p[target.param_idx][2]                                # release_rate
#             end
#         elseif s isa mRNA
#             du[i] = (
#                 u[s.transcriptor.idx] * p[s.param_idx][1]   # transcription_rate
#                 - u[s.idx] * p[s.param_idx][2]              # degradation_rate
#             )
#         elseif s isa Protein
#             du[i] = u[s.translator.idx] * p[s.param_idx][1]
#             for target in s.targets
#                 du[i] -= u[s.idx] * u[target.idx] * p[target.param_idx][1]  # binding_rate
#                 du[i] += u[target.idx] * p[target.param_idx][2]             # release_rate
#             end
#         elseif s isa Complex
#             du[i] = (
#                 u[s.species1.idx] * u[s.species2.idx] * p[s.param_idx][1]   # binding_rate
#                 - u[s.idx] * p[s.param_idx][2]                              # release_rate          
#             )
#         end
#     end
#     nothing
# end

