using Revise, DifferentialEquations, Plots
plotlyjs()

include("BioNet.jl")
include("Templates.jl")
include("Observations.jl")

obs = Observations("data.csv")
net = BioNetwork(true_network(), obs)


plot(solve_ODE(net, obs.t), labels=labels(net))


fun(θ) = loss!(net, θ)

params = θ(net)

@time fun(params)

optimize(fun, params, SimulatedAnnealing())