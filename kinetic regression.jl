using Plots

T = [336, 345, 377, 395, 440, 454]
kg = [7.33, 7.46, 10.48, 15.17, 19.01, 21.99] # cm/min
De = [7.407, 7.457, 7.576, 7.675, 7.752, 7.825] # cm2/min
ks = [2.47, 5.24, 9.37, 15.63, 25.38, 30.40] # cm/min

plot_kine_params = scatter(1 ./T, log.(De))
scatter!(1 ./T, log.(ks))
scatter!(1 ./T, log.(kg))


plot_kg = scatter(T .^ 1.75 + T .^ 0.88,kg)


# First-order decay model (Thein, 2000)
function firstorder!(du, u, p, t)
    k = p
    du[1] = -k * u[1]
end

# Parameters
temp_list = [600, 610, 620, 650] .+ 273.15 # K
k = [2.53e-3, 6.29e-3, 17e-3, 90e-3]        # min^-1 (600, 610, 620, 650 °C)

plot_UO3 = scatter(1 ./ temp_list, -log.(k))


using LinearAlgebra
using Statistics
using Plots

# Data from Table 2.9
T_C   = [600, 610, 620, 650]
T_K   = T_C .+ 273.15
k     = [2.53e-3, 6.29e-3, 17.0e-3, 90.0e-3]

x = 1 ./ T_K           # 1/T (K⁻¹)
y = log.(k)            # ln(k)

# Linear regression: y = m x + b
n = length(x)
m = (n*sum(x.*y) - sum(x)*sum(y)) / (n*sum(x.^2) - sum(x)^2)   # slope
b = (sum(y) - m*sum(x)) / n                                    # intercept

E_a = -m * 8.314 / 1000          # kJ/mol
A   = exp(b)                     # min⁻¹

println("Slope (m)       = ", round(m, sigdigits=6))
println("Intercept (b)   = ", round(b, sigdigits=6))
println("Activation energy Eₐ = ", round(E_a, sigdigits=4), " kJ/mol")
println("Pre-exponential A   = ", round(A, sigdigits=4), " min⁻¹")
println("R² (coefficient of determination):")

# R²
y_mean = mean(y)
SST = sum((y .- y_mean).^2)
SSR = sum((y .- (m.*x .+ b)).^2)
R² = 1 - SSR/SST
println("R² = ", round(R², sigdigits=5))

# Plot
p = scatter(x, y,
    label="data",
    xlabel="1/T  (K⁻¹)",
    ylabel="ln(k)  (min⁻¹)",
    title="Arrhenius plot – Type III UO₃ → U₃O₈",
    legend=:topright)

plot!(x, m.*x .+ b, label="fit", lw=2, color=:red, ls=:dash)
display(p)

# Optional: predicted k at 630 °C for example
T_test = 630 + 273.15
k_pred = A * exp(-E_a*1000 / (8.314 * T_test))
println("Predicted k at 630 °C = ", round(k_pred, sigdigits=4), " min⁻¹")