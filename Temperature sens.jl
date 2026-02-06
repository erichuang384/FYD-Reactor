using DifferentialEquations
using Plots
using Roots   # lightweight root-finding package

# Consecutive first-order reactions: A → B → C
# dA/dt = -k1 A
# dB/dt =  k1 A - k2 B
# dC/dt =  k2 B

function two_step!(du, u, p, t)
    A, B, C = u
    k1, k2 = p
    du[1] = -k1 * A
    du[2] =  k1 * A - k2 * B
    du[3] =  k2 * B
end

# Time span (minutes) – adjust as needed
tspan = (0.0, 10.0)

# Initial condition: start with pure UNH, [A]_0 = 1 (fraction), B=C=0
u0 = [1.0, 0.0, 0.0]

# Example: use parameters close to the given Arrhenius equations
# Let's take T = 300 °C = 573 K (within both ranges)
R = 1.987   # cal/mol·K 
T = 350.0 + 273.15
n_points = 500
T_range = LinRange(240,450,n_points) .+ 273.15
Ea1 = 11200.0     # cal/mol
A1  = 2.96e5      # min⁻¹
k1  = A1 .* exp.(-Ea1 ./ (R .* T_range))

Ea2 = 30800.0
A2  = 2.08e11
k2  = A2 .* exp.(-Ea2 ./ (R .* T_range))

C_target = 0.995

# Closed-form expression for C(t)
t_target = zeros(n_points)
for i in 1:n_points
    C_exact(t) = 1 - (k2[i]*exp(-k1[i]*t) - k1[i]*exp(-k2[i]*t)) / (k2[i] - k1[i])

# Algebraic equation to solve
    f(t) = C_exact(t) - C_target

# Solve for t (provide a reasonable bracket in minutes)
    t_target[i] = find_zero(f, (0.0, 1000.0))
end

plot1 = plot(T_range,t_target)

# higher temperature - smaller vessel, higher heat duty probably, 

# 