using Revise, DifferentialEquations, Plots, BenchmarkTools, Random
plotlyjs()

include("BioNet.jl")
include("Templates.jl")

obs = Observations("data.csv")
net = true_network()
init!(net, obs)

function LHS(bins, n::Int64)::Vector{Float64}
    is = randperm(MersenneTwister(2), length(bins))
    bins[is[1:n]]
end

θ = LHS(0:0.005:0.3, net.num_params)
println("loss: $(loss(net, θ))")
@btime ∂u(net, θ)

function desc()
    learning_rate = 0.0005
    θ = LHS(0:0.005:0.3, net.num_params)
    for i in 1:100
        g = ∂u(net, θ)
        θ = max.(0.0, θ - learning_rate * g)
        # println("loss: $(loss(net, θ))")
    end
end


@btime desc()

function global_fit(net::BioNetwork, epochs::Int64, max_iter=100)::Vector{Float64}
    min_loss::Float64 = Inf
    best_θ::Vector{Float64} = zeros(net.num_params)

    learning_rate = 0.01

    for epoch in 1:epochs    
        θ = LHS(0:0.005:0.3, net.num_params)
        for i in 1:max_iter
            g = ∂u(net, θ)
            θ = max.(0.0, θ - learning_rate * g)
            println(i)
        end
        L = loss(net, θ)
        println("Epoch $epoch, loss: $L")
        if L < min_loss
            best_θ = θ
            min_loss = L
        end
    end
    best_θ
end

global_fit(net, 1)


s = net[1]

toM(sub_u.(s.targets, repeat([u], length(s.targets))))

u