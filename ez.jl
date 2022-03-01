using Revise, DifferentialEquations, Plots, Distributions
using Optim

plotlyjs()

t_data = [0.25, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12] # time points
Y_data = [
    [0.43363251, 0.87767177, 0.52257947]
    [0.59829549, 0.68515400, 0.85511583]
    [0.57771761, 0.44965359, 1.32524242]
    [0.35806804, 0.34666734, 1.44662417]
    [0.14424849, 0.33451828, 1.54588553]
    [0.39095813, 0.13332767, 1.34309046]
    [0.55361681, 0.13633700, 1.07001529]
    [1.00481754, 0.12037373, 0.53871530]
    [1.26352803, 0.19847235, 0.49535966]
    [1.17064081, 0.10360132, 0.31091261]]

y_data = transpose(reshape(Y_data, 3,10))

plot(t_data, transpose(reshape(Y_data, 3,10)))


function create_A_b(M, theta)
    A = zeros(Float64, (3,3))
    b = theta[1:3]
    c = 4
    for i in 1:3
        for j in 1:3
            if M[i,j] == 1
                A[i,j] = theta[c]
                c += 1
            end
        end
    end
    return A, b
end

m1 = [[1,0,1] [0, 1, 0] [1, 0, 1]]
m2 = [[1,0,0] [0, 1, 0] [0, 0, 1]]
theta = [1, 2, 3, 4, 5, 6, 8, 9]
theta2 = [1, 2, 3, 4, 5, 6]
A, b =create_A_b(m1, theta)

x_init = [0, 1.0, 0]

function f!(du, u, p, t)
    du[1] = p[1][1,1] * u[1] + p[1][1,2] * u[2] + p[1][1,3] * u[3] + p[2][1]
    du[2] = p[1][2,1] * u[1] + p[1][2,2] * u[2] + p[1][2,3] * u[3] + p[2][2]
    du[3] = p[1][3,1] * u[1] + p[1][3,2] * u[2] + p[1][3,3] * u[3] + p[2][3]
end

function to_m(v::Vector{Vector{Float64}})
    m = zeros(Float64, length(v), length(v[1]))
    for (i,col) in enumerate(v)
        m[i,:] = v[i]
    end
    m
end

function solve_system(M, theta, t)
    prob = ODEProblem(f!, x_init, (0.0, maximum(t_data)), create_A_b(M, theta))
    sol = solve(prob)
    to_m(sol.(t))
end

m = solve_system(m1, theta, t_data)
solve_system(m2, theta2, t_data)

function log_likelihood(M, theta)
    σ = 0.15
    mu = solve_system(M, theta, t_data)
    t1 = - log(σ * sqrt(2.0 * π))
    t2 = - 1.0 / (2.0 * σ^2)*((y_data - mu).^2)
    terms = t1 .+ t2
    sum(terms)
end

log_likelihood(m2, theta2)

function residuals(M, theta)
    mu = solve_system(M, theta, t_data)
    sum((y_data - mu).^2)
end


function fit_model(M)
    nrml = Normal(0.05, 0.01)
    d = sum(M)+3
    θ₀ = rand(nrml, d)
    fun(theta) = residuals(M, theta)
    optimize(fun, θ₀)
end

m3 = [[1,1,1] [1,1,1] [1,1,1]]
m4 = [[1,0,1] [1,1,1] [1,1,1]]

tt = 0:0.1:maximum(t_data)

function test_model(M)
    plot(tt, solve_system(M, Optim.minimizer(fit_model(M)), tt))
    scatter!(t_data, y_data)
end

test_model(m4)
