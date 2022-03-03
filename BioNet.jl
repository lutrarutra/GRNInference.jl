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
    elseif ttt <: GeneProtein
        push!(net.species, GeneProtein(name, length(net.species)+1, kwargs[:gene], kwargs[:protein]))
    elseif ttt <: ProteinComplex
        push!(net.species, ProteinComplex(name, length(net.species)+1, kwargs[:protein1], kwargs[:protein2]))
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

# function make_ode(net::BioNetwork)
#     ODEForwardSensitivityProblem(f!, state₀(net), (0.0, maximum(obs.t)), params(net))
# end

# function system!(du, u, p, t, net::BioNetwork)
#     # du = zeros(Float64, length(network.species))
#     update!(net, state(net), p)
#     return d_dt!(net, du)
# end

# function update!(net::BioNetwork, u, θ)
#     params!(net, θ)
#     for (i,s) in enumerate(net.species)
#         s.u = u[i]
#     end
# end


function d_dt!(du, u, p, t, net::BioNetwork)
    for (i,s) in enumerate(net.species)
        du[i] = 0.0
        if s isa Gene
            for gp in s.gene_proteins
                du[i] -= u[s.idx] * u[gp.protein.idx] * p[gp.param_idx][1]  # binding_rate
                du[i] += u[gp.idx] * p[gp.param_idx][2]                     # release_rate
            end
        elseif s isa mRNA
            du[i] = (
                u[s.transcriptor.idx] * p[s.param_idx][1]   # transcription_rate
                - u[s.idx] * p[s.param_idx][2]              # degradation_rate
            )
        elseif s isa Protein
            du[i] = u[s.translator.idx] * p[s.param_idx][1]
            for target in s.targets
                du[i] -= u[s.idx] * u[target.idx] * p[target.param_idx][1]  # binding_rate
                du[i] += u[target.idx] * p[target.param_idx][2]             # release_rate
            end
        elseif s isa GeneProtein
            du[i] = (
                u[s.protein.idx] * u[s.gene.idx] * p[s.param_idx][1]    # binding_rate
                - u[s.idx] * p[s.param_idx][2]                          # release_rate
            )
        elseif s isa ProteinComplex
            du[i] = (
                u[s.protein1.idx] * u[s.protein2.idx] * p[s.param_idx][1]   # binding_rate
                - u[s.idx] * p[s.param_idx][2]                              # release_rate          
            )
        end
    end
end

# function d_dt!(net::BioNetwork, du)
#     du .= d_dt.(net.species)
# end

function params(net::BioNetwork)::Vector{Float64}
    collect(Iterators.flatten(params.(net.species)))
end

function params!(net::BioNetwork, θ)
    c = 1
    for (i,s) in enumerate(net.species)
        n = length(params(s))
        params!(s, θ[c:c+n-1])
        c+=n
    end
end

function solve_ode(net::BioNetwork, θ; solver)
    f!(du,u,p,t) = d_dt!(du,u,p,t,net)
    prob = ODEProblem(f!, state₀(net), (0.0, maximum(net.obs.t)), θ)
    solve(prob, solver)
end


function loss(net::BioNetwork, θ; solver)::Float64
    f!(du,u,p,t) = d_dt!(du,u,p,t,net)
    prob = ODEProblem(f!, state₀(net), (0.0, maximum(net.obs.t)), θ)
    sol = solve(prob, solver)
    u = sol.(net.obs.t)

    L = zeros(Float64, length(net.ix))
    for (t, u_t) in enumerate(u)
        for (i, k) in enumerate(net.ix)
            L[i] += (u_t[k] - net.obs.u[t, i])^2
        end
    end
    sum(L)
end

function state₀(net::BioNetwork)
    u = zeros(Float64, length(net.species))
    for i in 1:length(net.species)
        u[i] = net.species[i].u₀
    end
    u
end


function ∂u(net::BioNetwork, θ; solver)
    N::Int64 = length(θ)
    M::Int64 = size(net.obs.u)[2]

    # _prob = remake(net.ode_prob, u0=net.ode_prob.u0, p=θ)
    f!(du,u,p,t) = d_dt!(du,u,p,t,net)
    _prob = ODEForwardSensitivityProblem(f!, copy(state₀(net)), (0.0, maximum(obs.t)), θ)
    sol = solve(_prob, solver, saveat=net.obs.t[2:end])
    x, dp = extract_local_sensitivities(sol)
    
    d = zeros(Float64, (M, length(net.obs.t)-1))

    for (i, k) in enumerate(net.ix)
        d[i,:] = net.obs.u'[i,2:end] - x[k,:]
    end
    
    GD = zeros(Float64, N)
    
    for i in 1:N
        GD[i] = -2.0*dot(dp[i][net.ix,:], d)
    end
    # dp1 = dp[1][:,2:end]
    # dp2 = dp[2][:,2:end]
    # [-dot(dp1, d)*2, -dot(dp2, d)*2]
    GD
end

function GD(net::BioNetwork; learning_rate=0.001, max_iter=1000)
    θ₀ = params(net)
    path = zeros(Float64, (max_iter+1, length(θ₀)))
    path[1,:] = θ₀
    for i in 1:max_iter
        k = path[i,:]
        g = ∂u(net, k)
        
        k = max.(0.0, k - learning_rate * g)
        # println(learning_rate*g)
        println(loss(net, k))
        path[i+1,:] = k
    end
    path
end

