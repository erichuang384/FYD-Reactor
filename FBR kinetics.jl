using DifferentialEquations
using Plots

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
T = 313.0 + 273.15
Ea1 = 11200.0     # cal/mol
A1  = 2.96e5      # min⁻¹
k1  = A1 * exp(-Ea1 / (R * T))

Ea2 = 30800.0
A2  = 2.08e11
k2  = A2 * exp(-Ea2 / (R * T))

# Solve
p = [k1, k2]
prob = ODEProblem(two_step!, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)



# Plot fractions vs time
p1 = plot(sol, vars=(1), label="UNH (A)",           lw=2.5, color=:darkblue)
plot!(sol, vars=(2), label="Monohydrate (B)", lw=2.5, color=:orange)
plot!(sol, vars=(3), label="UO₃ (C)",             lw=2.5, color=:green3)

plot!(xlabel="Time (min)", ylabel="Fraction", 
      title="Two-step decomposition of UNH → UO₃ at 300°C",
      ylims=(0,1.05), legend=:outertopright, grid=true, dpi=120)

