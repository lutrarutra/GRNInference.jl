using Optim, DifferentialEquations

include("Observations.jl")

include("Species.jl")


mutable struct BioNetwork
    species::Vector{Species}
    obs::Observations
    ix::Vector{Int64}
    idx::Dict{Species, Int64}

    BioNetwork() = new([], nothing, [])
    function BioNetwork(species::Vector{Species}, obs::Observations)
        net = new(species, obs)
        net.ix = id_labels(net)
        for (i, k) in enumerate(net.ix)
            net.species[k].u₀ = obs.u[1][i]
            net.species[k].u = obs.u[1][i]
        end
        net
    end
end

function d_dt!(net::BioNetwork, du::Vector{Float64})::Vector{Float64}
    du .= d_dt.(net.species)
    du
end

function θ(net::BioNetwork)
    collect(Iterators.flatten(θ.(net.species)))
end

function θ!(net::BioNetwork, vals::Vector{Float64})
    c = 1
    for (i,s) in enumerate(net.species)
        n = length(θ(s))
        θ!(s, vals[c:c+n-1])
        c+=n
    end
end

function loss!(net::BioNetwork, θ::Vector{Float64})::Float64
    θ!(net, θ)
    prob = ODEProblem(system!, state₀(net), extrema(net.obs.t), net)
    sol = solve(prob)

    u = sol.(obs.t)

    L = zeros(Float64, length(net.ix))
    for (t, u_t) in enumerate(u)
        for (i, k) in enumerate(net.ix)
            L[i] += (u_t[k] - net.obs.u[t][i])^2
        end
    end
    sum(L)
end

function fit!(net::BioNetwork)
    θ₀ = θ(net)
    fun(θ) = loss!(net, θ)
    optimize(fun, θ₀)
end

function update!(net::BioNetwork, u::Vector{Float64})
    for (i,s) in enumerate(net.species)
        s.u = u[i]
    end
end

function system!(du::Vector{Float64}, u::Vector{Float64}, net::BioNetwork, t)::Vector{Float64}
    # du = zeros(Float64, length(network.species))
    update!(net, u)
    return d_dt!(net, du)
end

function id_labels(net::BioNetwork)::Vector{Int64}
    ix::Vector{Int64} = []
    for (i,s) in enumerate(net.species)
        for label in obs.labels
            if lowercase(s.name) == lowercase(label)
                push!(ix, i)
                break
            end
        end
    end
    ix
end


function state(net::BioNetwork)
    u = zeros(Float64, length(net.species))
    for i in 1:length(net.species)
        u[i] = net.species[i].u
    end
    u
end

function state₀(net::BioNetwork)
    u = zeros(Float64, length(net.species))
    for i in 1:length(net.species)
        u[i] = net.species[i].u₀
    end
    u
end


function labels(net::BioNetwork)::Matrix{String}
    ls::Vector{String} = []
    for i in 1:length(net.species)
        push!(ls, net.species[i].name)
    end
    reshape(ls, 1, length(ls))
end


