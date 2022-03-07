using Optim, DifferentialEquations, DiffEqSensitivity, LinearAlgebra, JSON

include("Species.jl")
include("Observations.jl")


mutable struct BioNetwork
    ix::Vector{Int64}
    species::Vector{Species}
    num_params::Int64
    N::Int64
    obs::Observations
    ode_prob
    θ::Vector{Float64}
    desc::String

    BioNetwork() = new([], [], 0)
end

Base.getindex(net::BioNetwork, idx::Int64)::Species = net.species[idx]
Base.getindex(net::BioNetwork, name::String)::Species = net.species[findfirst(s -> lowercase(s.name) == lowercase(name), net.species)]


function add_species!(net::BioNetwork, ttt::DataType, name::String)
    @assert ttt <: Species
    push!(net.species, ttt(name, length(net.species)+1))
    # if ttt <: Gene
    # elseif ttt <: mRNA
    #     push!(net.species, mRNA(name, length(net.species)+1))
    # elseif ttt <: Protein
    #     push!(net.species, Protein(name, length(net.species)+1))
    # elseif ttt <: Complex
    #     push!(net.species, Complex(name, length(net.species)+1))
    # end
    net.num_params += net.species[end].num_params
end

function init!(net::BioNetwork, obs::Observations)
    net.obs = obs
    net.ix = id_labels(net)
    for (i, k) in enumerate(net.ix)
        net.species[k].u₀ = net.obs.u[1,i]
    end

    net.N = length(net.species)

    f!(du, u, p, t) = d_dt!(net, du, u, p)
    net.ode_prob = ODEForwardSensitivityProblem(f!, state₀(net), (0.0, maximum(obs.t)), zeros(net.num_params))
    nothing
end

function init_params(net::BioNetwork, θ)
    θ_i = 1
    transcription = zeros(eltype(θ), net.N, net.N)
    translation = zeros(eltype(θ), net.N, net.N)
    degradation = zeros(eltype(θ), net.N)
    release = zeros(eltype(θ), net.N, net.N)
    binding = zeros(eltype(θ), net.N, net.N, net.N)

    for s in net.species
        if s isa mRNA
            transcription[s.transcriptor.idx, s.idx] = θ[θ_i]
            degradation[s.idx] = -θ[θ_i+1]
        elseif s isa Protein
            translation[s.translator.idx, s.idx] = θ[θ_i]
        elseif s isa Complex
            # s_i ⋅ s_j -> c_k <-> s_j ⋅ s_i -> c_k
            binding[s.species1.idx, s.species2.idx, s.idx] = θ[θ_i]
            binding[s.species2.idx, s.species1.idx, s.idx] = θ[θ_i]

            binding[s.species1.idx, s.species1.idx, s.species1.idx] = -θ[θ_i]
            binding[s.species1.idx, s.species2.idx, s.species2.idx] = -θ[θ_i]
            binding[s.species2.idx, s.species1.idx, s.species1.idx] = -θ[θ_i]
            binding[s.species2.idx, s.species2.idx, s.species2.idx] = -θ[θ_i]
            # binding[s.species2.idx, s.species1.idx, s.species1.idx] = -θ[θ_i]
            # binding[s.species2.idx, s.species1.idx, s.species2.idx] = -θ[θ_i]
            # binding[s.species2.idx, s.species2.idx, s.species1.idx] = -θ[θ_i]
            # binding[s.species2.idx, s.species2.idx, s.species2.idx] = -θ[θ_i]

            release[s.idx, s.species1.idx] = θ[θ_i+1]
            release[s.idx, s.species2.idx] = θ[θ_i+1]
            release[s.idx, s.idx] = -θ[θ_i+1]
        end
        θ_i += s.num_params
    end
    # println("crip: $(transcription[:,1])")
    # println("lat: $(translation[:,1])")
    # println("rel: $(release[:,1])")
    # println("deg: $(degradation[1])")
    (transcription, translation, degradation, release, binding)
end

function d_dt!(net, du, u, p)
    transcription, translation, degradation, release, binding = init_params(net, p)
    du .= zeros(Float64, net.N)
    # uu = u'.*u

    for i in 1:net.N
        # 0th order
        du[i] += u[i] * degradation[i]
        # 1st order
        du .= du + u[i] * (transcription[i,:] + translation[i,:] + release[i,:])
        # du .= du + u[i] * transcription[i,:]
        # du .= du + u[i] * translation[i,:]
        # du .= du + u[i] * release[i,:]
        # if du[1] > 0
        #     println("crip: $(u[i] * transcription[i,:])")
        #     println("lat: $(u[i] * translation[i,:])")
        #     println("rel: $(u[i] * release[i,:])")
        # end

        # 2nd order
        for s1 in 1:net.N
            for s2 in 1:s1
                du[i] += u[s1] * u[s2] * binding[s1, s2, i]
            end
        end
    end
    # println(du[1:5])
    nothing
end


function id_labels(net::BioNetwork)::Vector{Int64}
    ix::Vector{Int64} = []
    for (j,label) in enumerate(net.obs.labels)
        for (i,s) in enumerate(net.species)
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

function BIC(net::BioNetwork, L)
    σ = mean(obs.u) * 0.1    # noise magnitude 10% of measurement
    t1 = - sum(log(σ * sqrt(2.0 * π)))
    t2 = - 1.0 / (2.0 * σ^2)*L
    log_likelihood = t1 + t2
    d = net.num_params
    log_likelihood - 0.5 * d * log(length(net.obs.t))
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

    _prob = remake(net.ode_prob, u0=net.ode_prob.u0, p=θ)
    sol = solve(_prob, BS3(), saveat=net.obs.t[2:end]; abstol=1e-6, reltol=1e-6, maxiters=1000)
    x, dp = extract_local_sensitivities(sol)
    
    d = zeros(Float64, (M, length(net.obs.t)-1))

    for (i, k) in enumerate(net.ix)
        d[i,:] = net.obs.u'[i,2:end] - x[k,:]
    end
    
    GD = zeros(Float64, N)
    
    for i in 1:N
        GD[i] = -2.0*dot(dp[i][net.ix,:], d)
    end

    GD, sum(d.^2)
end

function save(net::BioNetwork, θ::Vector{Float64}, file::String)
    open(file, "w") do io
       for s in net.species
            write(io, csv(s)*"\n")
       end
       write(io, JSON.json(θ))
    end
    nothing
end

function load(file::String, obs::Observations)
    net = BioNetwork()
    net.desc = string(split(file, ".")[1])
    lines = readlines(file)
    θ::Vector{Float64} = []

    transcriptors::Vector{Tuple{String, String}} = []
    translators::Vector{Tuple{String, String}} = []
    reactants::Vector{Tuple{String, String, String}} = []

    for line in lines
        data = split(line, ",")
        if data[1] == "Gene"
            add_species!(net, Gene, string(data[3]))
        elseif data[1] == "mRNA"
            add_species!(net, mRNA, string(data[3]))
            m = net.species[end]
            push!(transcriptors, (m.name, string(data[5])))
        elseif data[1] == "Protein"
            add_species!(net, Protein, string(data[3]))
            p = net.species[end]
            push!(translators, (p.name, string(data[5])))
        elseif data[1] == "Complex"
            add_species!(net, Complex, string(data[3]))
            c = net.species[end]
            push!(reactants, (c.name, string(data[5]), string(data[6])))
        else
            θ = JSON.parse(line)
        end
    end

    for (m, g) in transcriptors
        add_transcriptor!(net[m], net[g])
    end
    for (p, m) in translators
        add_translator!(net[p], net[m])
    end
    for (s, s1, s2) in reactants
        add_reactants!(net[s], net[s1], net[s2])
    end
    init!(net, obs)
    (net, θ)
end