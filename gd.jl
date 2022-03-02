using Revise, DifferentialEquations, Plots, Distributions
using Optim, DiffEqSensitivity, Calculus, ForwardDiff, Tracker
using BenchmarkTools, LinearAlgebra

M(x) = mapreduce(permutedims, vcat, x)

t_data = range(1, 10, 20)
Y_data = [
    [2.46453992, 3.60858988]
    [1.36147358, 2.94974486]
    [0.83659004, 2.42189492]
    [0.63479081, 1.21896289]
    [1.96550665, 0.79799601]
    [2.63884836, 1.05561561]
    [3.9317634, 1.85744158]
    [2.37791467, 3.43918307]
    [1.63918609, 3.56459394]
    [1.13540439, 2.03678334]
    [1.30275284, 1.52561854]
    [1.62360963, 1.17760926]
    [2.68444394, 0.98298436]
    [3.83887552, 1.88928201]
    [2.23165204, 3.46577548]
    [1.07424108, 2.4358289 ]
    [0.57951303, 1.70852535]
    [1.38550887, 1.13571176]
    [1.91259666, 1.06761892]
    [2.90236636, 1.2190146 ]]

y_data = collect(transpose(reshape(Y_data, 2,20)))

function d_dt!(du, u, p, t)
    du[1] = 2.0 * u[1] - p[1] * u[1] * u[2]
    du[2] = p[1]*u[1]*u[2] - p[2]*u[2]
    nothing
end


function test_f(p)
    prob = ODEProblem(d_dt!, u0, extrema(t_data), p)
    solve(prob)[end]
end

u0 = [2.5, 1.0]
p = [1.0, 1.2]


ts = 1:0.5:10

function target_function(params::Vector{Float64})

    prob = ODEProblem(d_dt!, u0, extrema(t_data), params)
    sol = solve(prob)
    u = sol(t_data).u
    s = 0.0
    for (i,v) in enumerate(u)
        s += (y_data[i, 1] - v[1])^2 + (y_data[i, 2] - v[2])^2
    end
    s
end

function extended!(du, u, k, t)
    # Jx[1 3 2 4]
    Jx = [[2 - k[1] * u[2], k[1] * u[2]] [-k[1] * u[1], k[1] * u[1] - k[2]]]
    Jk = [[-u[1] * u[2], u[1] * u[2]] [0, -u[2]]]
    S = [[u[3], u[5]] [u[4], u[6]]]

    dS = Jx*S + Jk
    # du = zeros(Float64, 2)
    d_dt!(du, u[1:2], k, nothing)

    du .= [du[1], du[2], dS[1,1], dS[1,2], dS[2,1], dS[2,2]]
end


function target_grad(k)
    N = length(t_data)
    s0 = zeros(Float64, 4)
    u = vcat(u0, s0)
    t_ext = vcat(0.0, t_data)
    # ext(θ) = extended(θ, k)
    prob = ODEProblem(extended!, u, extrema(t_ext), k)
    sol = solve(prob)
    X_tilde = M(sol.(t_ext))[2:end,:]
    
    X = X_tilde[:,1:2]
    S = X_tilde[:,3:6]
    
    d = y_data - X
    g = zeros(Float64, (2,1))
    
    for i in 1:N
        S_i = [[S[i,1], S[i,2]] [S[i,3], S[i,4]]]
        g += -2*(S_i*d[i,:])
    end
    g
end



function GD(k0, n_iter, r)
    path = zeros(Float64, (n_iter+1, length(k0)))
    path[1,:] = k0
    for i in 1:n_iter
        k = path[i,:]
        g = G(k)
        if sqrt(sum(g.^2)) < 0.0001
            return path[begin:i,:]
        end
        k -= r * g
        path[i+1,:] = k
    end
    path
end

@btime GD([2.0, 1.5], 50, 0.001)

prob  = ODEForwardSensitivityProblem(d_dt!, u0, extrema(t_ext), p)

function G(p)
    t_ext = vcat(0.0, t_data)
    _prob = remake(prob, u0=convert.(eltype(p),prob.u0), p=p)
    sol = solve(_prob, saveat=t_ext)
    x, dp = extract_local_sensitivities(sol)

    dp1 = dp[1][:,2:end]
    dp2 = dp[2][:,2:end]
    d = y_data' - x[:,2:end]

    [-dot(dp1, d)*2, -dot(dp2, d)*2]
end


@btime G(p)

