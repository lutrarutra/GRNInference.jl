using Revise, DifferentialEquations, Plots, BenchmarkTools, Random
plotlyjs()

include("BioNet.jl")
include("Templates.jl")

obs = Observations("data.csv")
net = true_network()
init!(net, obs)


function LHS(bins, n::Int64)::Vector{Float64}
    is = randperm(MersenneTwister(), length(bins))
    bins[is[1:n]]
end

θ = LHS(0:0.005:0.3, net.num_params)
@time ∂u(net, θ)

function global_fit(net::BioNetwork, epochs::Int64, max_iter=100)::Vector{Float64}
    min_loss::Float64 = Inf
    best_θ::Vector{Float64} = zeros(net.num_params)

    solver=BS3()
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

p = LHS(0:0.005:0.2, net.num_params)

u = state₀(net)
s = net["gCbf1"]

du = s.u₀
for gp in s.targets
    du -= u[s.idx] * u[gp.protein.idx] * p[gp.param_idx][1]  # binding_rate
    du += u[gp.idx] * p[gp.param_idx][2]                     # release_rate
end
du


state(s, u) * dot(u[idx.(s.targets)], p[idx.(s.targets)])

s.idx

@which u[]
getindex.([u,u], sub_idx.(s.targets))

s.targets

u[idx.(s.targets)]
state(s, u)

