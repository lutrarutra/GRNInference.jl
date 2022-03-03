using Revise, DifferentialEquations, Plots, BenchmarkTools, Distributions
plotlyjs()

NORMAL = Normal(0.01, 0.005)

include("BioNet.jl")
include("Templates.jl")

obs = Observations("data.csv")
net = true_network()
init!(net, obs)

solver=Vern7()
learning_rate = 0.05

θ₀ = rand(NORMAL, net.num_params)
@time ∂u(net, θ₀, solver=solver)

path = zeros(Float64, (100+1, length(θ₀)))
path[1,:] = θ₀

@time for i in 1:100
    k = path[i,:]
    g = ∂u(net, k, solver=solver)
    
    k = max.(0.0, k - learning_rate * g)
    # k = k - learning_rate * g

    # println(learning_rate*g)
    println(loss(net, k, solver=solver))
    path[i+1,:] = k
end


loss(net, path[end,:], solver=solver)

sol = solve_ode(net, path[end, :], solver=solver)
sol.u

plot(sol, labels=labels(net))

scatter!(obs.t, obs.u, labels=reshape(obs.labels, 1, length(obs.labels)))
labels(net)

labels(net)

