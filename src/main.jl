using Revise, DifferentialEquations, Random

println("Using $(Threads.nthreads()) threads")

include("BioNetwork.jl")
include("Builder.jl")

obs = Observations("data.csv")
# net = true_network()
# init!(net, obs)


function global_fit(net::BioNetwork, epochs::Int64, max_iter=200)
    # params::Vector{Tuple{Vector{Float64}, Float64}} = [] 
    params = Array{Tuple{Vector{Float64}, Float64}}(undef, epochs)

    Threads.@threads for epoch in 1:epochs

        learning_rate = 0.05
        θ = rand(net.num_params) * 0.1
        g = zeros(net.num_params)

        L = 10.0
        L_best = 10.0
        θ_best = θ
        L_prev::Float64 = Inf
        for i in 1:max_iter
            try
                g, L = ∂u(net, θ)
                θ = max.(0.0, θ - learning_rate * g)
            catch e
                println(e)
                println("Premature GD stop due to error...")
                break
            end
            if L < L_best
                L_best = L
                θ_best = θ
            end
            if L_prev < L
                learning_rate *= 0.9
                println("Reduced learning rate by 20% to: $(learning_rate)")
            end
            L_prev = L
            println("Epoch: $epoch, iteration: $i, loss: $L_prev @$(net.desc)")
            flush(stdout)
        end
        L = loss(net, θ_best)
        println("Finished epoch $epoch with loss: $L")

        params[epoch] = (θ, L)
    end

    min_L = params[1][2]
    best_p = params[1][1]
    for i in 1:length(params)
        if params[i][2] < min_L
            min_L = params[i][2]
            best_p = params[i][1]
        end
    end
    println("Best loss: $min_L")
    best_p, BIC(net, min_L)
end

# @time θ, θ_bic = global_fit(net, 12)

function test_networks(nets::Vector{BioNetwork})
    for (net_i, net) in enumerate(nets)
        println("Fitting network structure $net_i")
        flush(stdout)
        init!(net, obs)
        θ, θ_bic = global_fit(net, 12)
        save(net, θ, "../networks/$(net.desc)_$net_i.grn")
    end
end

nets = Vector{BioNetwork}()

append!(nets, random_network.([9,10,10,10]))

test_networks(nets)

print("DONE!")




