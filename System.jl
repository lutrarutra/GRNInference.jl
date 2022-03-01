using Revise, DifferentialEquations, Plots
plotlyjs()

include("BioNet.jl")
include("Templates.jl")
include("Observations.jl")

obs = Observations("data.csv")
net = BioNetwork(true_network(), obs)

θ(net)

fit!(net)


tt = 0:0.1:maximum(obs.t)
plot(sol(tt), label=labels(net))

loss(net, sol.(obs.t))

state₀(net)


[1 2 3 4 5 6 7 8][2:2]
[1 2 3 4 5 6 7 8][4:6]
