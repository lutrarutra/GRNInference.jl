using Revise, Plots, Statistics
plotlyjs()

include("BioNetwork.jl")
include("Observations.jl")

obs = Observations("data.csv")


function plot_net(net, θ; all=false, kwargs...)
    f!(du, u, p, t) = d_dt!(net, du, u, p)
    prob = ODEProblem(f!, state₀(net), extrema(net.obs.t), θ)
    sol = solve(prob, saveat=minimum(net.obs.t):1:maximum(net.obs.t))
    
    m = M(sol.u)
    # plot(sol, labels=labels(net), color=[1,2,3,4])

    plot(;title=kwargs[:title], titlefontsize=10)
    for (i, k) in enumerate(net.ix)
        plot!(sol.t, m[:, k], labels=labels(net)[k], color=i)
        scatter!(net.obs.t, net.obs.u[:, i], label="true " * net.obs.labels[i], color=i)
    end
    if all
        for i in 1:net.N
            if i ∉ net.ix
                plot!(sol.t, m[:, i], labels=labels(net)[i], color=i+length(net.ix))
            end
        end
    end
    plot!()
end
    
function plot_net_from_file(file::String)
    net, θ = load(file, obs)
    L = loss(net, θ)
    title="$(net.desc) BIC: $(round(BIC(net, L), digits=1)), SE: $(round(L, digits=5))"
    plot_net(net, θ; title=title)
end


function summarise_files(path::String)
    species_count = Dict{String, Int64}()
    for (i,file) in enumerate(readdir(path; join=true))
        net, p = load(file, obs)
        L = loss(net, p)
        bic = BIC(net, L)
        if bic > -300
            for s in net.species
                if s isa Complex
                    if !haskey(species_count, s.name)
                        species_count[s.name] = 0
                    end
                    species_count[s.name] += 1
                end
            end
        end
        println("Random model ($i) & $(length(net.species)) & $(net.num_params) & $(round(L, digits=5)) & $(round(bic,digits=1)) \\\\")
    end
    counts = sort(collect(species_count), by=last)
    for (key, value) in counts
        println("$key : $value")
    end
end


summarise_files("networks/")


plot_net_from_file("networks/random_network_4_9.grn")
plot_net_from_file("networks/simplified_1.grn")
plot_net_from_file("networks/simplified_2.grn")
plot_net_from_file("networks/random_network_8_23.grn")
plot_net_from_file("networks/random_network_8_21.grn")

