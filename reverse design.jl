"""
    fluidized_bed_conversion(;
        L_fluidized::Float64,     # fluidized bed height (m)
        u_b::Float64,             # bubble rise velocity (m/s)
        U::Float64,               # superficial gas velocity (m/s)
        U_mf::Float64,            # minimum fluidization velocity (m/s)
        d_b::Float64,             # bubble diameter (m)
        D::Float64,               # gas diffusivity (m²/s)
        ε_mf::Float64,            # voidage at minimum fluidization
        k::Float64,               # first-order rate constant (1/s)
        γ_b::Float64,             # volume of solids dispersed in bubbles (-), typically 0.001–0.01
        α::Float64 = 0.6          # cloud-wake parameter, typically 0.25–1.0
    ) -> Float64

Computes the reactant conversion X_A for a bubbling fluidized bed
using the three-phase (bubble-cloud-emulsion) model with first-order kinetics.

Returns: X_A ∈ [0,1]
"""
function fluidized_bed_conversion(;
    L_fluidized::Float64,     # bed height
#    u_b::Float64,             # bubble velocity
    U::Float64,               # superficial velocity
#    U_mf::Float64,            # min. fluidization velocity
#    d_b::Float64,             # equivalent bubble diameter
    D::Float64,               # molecular diffusivity of gas
    ε_mf::Float64,            # void fraction at mf
#    k::Float64,               # rate constant (s⁻¹)
    γ_b::Float64,             # solids fraction in bubble phase
    α::Float64,          # cloud-wake interchange parameter
    T::Float64                # K bed temperature
)::Float64

    density_solid = 7300 # kg/m3
    d_p = 200 * 10^(-6)
    U_mf = u_mf(d_p, density_solid, T)

    R = 1.987   # cal/mol·K 
    Ea2 = 30800.0
    A2  = 2.08e11
    k  = A2 * exp(-Ea2 / (R * T)) / 60 # s -1

    g = 9.81                  # gravity (m/s²)

    # Bubble rise velocity relative to emulsion
    d_b = 0.54 * (U - U_mf) ^ 0.4 * L_fluidized ^ 0.8 / g ^ 0.2

    u_br = 0.711 * √(g * d_b)

    u_b  = U - U_mf + u_br  # bubble velocity


    # γ_c ─ solids in cloud-wake phase
    γ_c = (1 - ε_mf) * (3 * (U_mf / ε_mf) / (u_br - U_mf / ε_mf)+ α)

    # δ ─ bubble fraction in bed (approximate form most commonly used)
    δ = (U - U_mf) / u_b
    # More exact form (uncomment if needed):
    # δ = (U - (1 - δ - α*δ)*U_mf) / u_b   → requires solver

    # γ_e ─ solids in emulsion phase
    γ_e = (1 - ε_mf) * (1 - δ) / δ - γ_b - γ_c

    # K_bc ─ bubble ↔ cloud exchange coefficient
    K_bc = 4.5 * (U_mf / d_b) + 5.85 * (D^0.5 * g^0.25) / (d_b^1.25)

    # K_ce ─ cloud ↔ emulsion exchange coefficient
    K_ce = 6.78 * ((ε_mf * D * u_b / d_b^3)^0.5)

    # χ₃ = 1/K_ce + 1/(γ_e k)
    χ₃ = 1/K_ce + 1/(γ_e * k)

    # χ₂ = γ_c k + 1/χ₃
    χ₂ = γ_c * k + 1/χ₃

    # χ₁ = 1/K_bc + 1/χ₂
    χ₁ = 1/K_bc + 1/χ₂

    # Overall term inside the exponential
    overall = (γ_b * k + 1/χ₁) * (L_fluidized / u_b)

    # Conversion for plug flow of bubbles (first-order reaction)
    X_A = 1 - exp(-overall)

    return X_A   # keep in physically meaningful range
end

# Example parameter set
params = Dict(
    :L_fluidized => 0.38,        # m
    :U           => 0.45,       # m/s
    :D           => 1.8e-4,     # m²/s  (typical for gases at ~300–500 K)
    :ε_mf        => 0.45,
    :γ_b         => 0.01,       # 0.001 to 0.01      
    :α           => 0.5,       # 0.25  to 1.0
    :T           => 450 + 273.15       
)

X = fluidized_bed_conversion(; params...)

println("Conversion X_A = ", round(X, digits=4))   
