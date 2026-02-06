# Pilot plant data
#
using Clapeyron, Roots, Polynomials
dia_pilot =   [0.15, 0.25, 1.17]  # m pilot plant diameter
len_pilot =   [0.38, 1.52, 3.7]            # m pilot reactor length
nozzle_location_pilot = [0.15, 1.07, 3.25] # was 3'6'', in m

T_pilot   = [450, 430, 288] .+ 273.15        # K pilot temperature
T_com = T_pilot

pilot_prod_rate_area = [730, 1470, 1420] # kg/(hr*m^2) first is 420-430

fluidizing_velocity = [0.45, 1.5, 2.0]    # m/s first is 0.4-0.55, second is 0.3-2.8

#### Scaling relations
#=
inertia/gravity force
u0^2/gL

particle inertia/gas viscous
ρ_s * u0 * dp / μ 

gas inertia/gas viscous
ρ_f * u0 * L / μ 

solid inertia/gas inertia
ρ_s/ ρ_f

bed heigh/bed diameter
L/D
=#

# Required production rate
com_prod_rate   = 0.810411 * (238.03 + 3 * 15.99) * 10 ^(-3) * 60 ^ 2  # kg/h of UNH

com_dia         = com_prod_rate./pilot_prod_rate_area

####### Apply scaling laws
# bed height/bed diameter
com_len = com_dia .* len_pilot ./ dia_pilot

com_nozzle = com_len/len_pilot * nozzle_location_pilot #  m nozzle location of full scale

# inertia / gravity forces for superficial gas velocity
#""The average superficial velocity of the fluidizing gas in the heater zone is 1.7 ft/sec and drop to a value of 1.25 ft/eec in the spray zone.""

pilot_u0 = fluidizing_velocity # m/s in AUSTRALIAN ATOMIC ENERGY COMMISSION also has u0 of 0.5 m/s

com_u0   = sqrt.(pilot_u0 .^ 2 ./ dia_pilot .* com_dia) # m/s

# particle inertia/gas viscous u0 * dp
# assume particle diameter of 0.01 inches, graph/table hard to figure out
pil_dp = 0.01 * 0.0254  # m
com_dp = pil_dp .* pilot_u0 ./ com_u0 # within graph


#### Check it is above minimium fluidization velocity
# Archimedes number

function u_mf(dp, density_solid, T)
    """
    minimum fluidization velocity with: dp - m, density_solid/gas - kg/m^3, viscosity_gas - m^2 / s, T - K
    Note it is kinetmatic viscosity of gas
    """
    p = 1e5
    air_model   = PR(["Nitrogen","Oxygen"];idealmodel = WalkerIdeal)
    z_air = [0.785, 0.215]  # mole fractions

    density_gas = mass_density(air_model,p,T,z_air) # kg/m3
    dynamic_viscosity_gas = 1.458 * 10 ^ (-6) * T ^ 1.5 * (1/ (T + 110.4)) # Pa s
    viscosity_gas = dynamic_viscosity_gas/density_gas
    Ar    = 9.81 * dp^3 * (density_solid - density_gas) / (viscosity_gas ^ 2 * density_gas)
    Re_mf = Ar / (1400 + 5.22 * sqrt(Ar))
    umf   = Re_mf * viscosity_gas / dp

    return umf
end

# assume solid density 4000 kg/m^3, density of gas based on PR air 
umf = u_mf.(com_dp, 4000,T_com)
println(com_u0>umf)


### elutriation velocity

function u_elu_Ar(dp, density_solid, T)
    """
    elutriation/terminal velocity from an Ar - Re correlation
    """
    p = 1e5
    air_model   = PR(["Nitrogen","Oxygen"];idealmodel = WalkerIdeal)
    z_air = [0.785, 0.215]  # mole fractions

    density_gas = mass_density(air_model,p,T,z_air) # kg/m3
    dynamic_viscosity_gas = 1.458 * 10 ^ (-6) * T ^ 1.5 * (1/ (T + 110.4)) # Pa s
    viscosity_gas = dynamic_viscosity_gas/density_gas
    Ar    = 9.81 * dp^3 * (density_solid - density_gas) / (viscosity_gas ^ 2 * density_gas)
    Re_elu = 1.74 * Ar / (31.3 + sqrt(Ar))
    uelu   = Re_elu * viscosity_gas / dp

    return uelu
end
#=
function u_elu_RE2(dp, density_solid, T)
    """
    elutriation/terminal velocity from RE2
    """
    p = 1e5
    air_model   = PR(["Nitrogen","Oxygen"];idealmodel = WalkerIdeal)
    z_air = [0.785, 0.215]  # mole fractions

    density_gas = mass_density(air_model,p,T,z_air) # kg/m3
    dynamic_viscosity_gas = 1.458 * 10 ^ (-6) * T ^ 1.5 * (1/ (T + 110.4)) # Pa s
    viscosity_gas = dynamic_viscosity_gas/density_gas

    #Re_gas  = uelu * dp / viscosity_gas
    Re_gas  = u_elu_RE2(dp, density_solid, T) * dp / visocisty_gas
    C_D     = exp(-5.5 + 69.43/ (log(Re_gas) + 7.99))
    uelu   = sqrt(4 * dp * (density_solid - density_gas)/ (3 * C_D * density_gas))
    println(uelu)

    return uelu
end
=#
using Roots   # Make sure to add Roots.jl: ]add Roots

using Roots   # ]add Roots if not already

function u_elu_2(dp, density_solid, T)
    """
    Terminal / elutriation velocity – non-iterative via quartic in z = √Re
    ξ(Re) = 21/Re + 6/√Re + 0.28
    """
    g = 9.81
    p = 1e5

    air_model = PR(["Nitrogen", "Oxygen"]; idealmodel = WalkerIdeal)
    z_air = [0.785, 0.215]

    ρg = mass_density(air_model, p, T, z_air)          # kg/m³
    μ  = 1.458e-6 * T^1.5 / (T + 110.4)               # Pa·s

    Δρ = density_solid - ρg
    if Δρ ≤ 0
        error("density_solid ($density_solid) ≤ ρg ($ρg) → no settling velocity")
    end

    Ar = ρg * Δρ * g * dp^3 / (μ^2)

    # Define the function we want f(z) = 0
    f(z) = 0.28 * z^4 + 6.0 * z^3 + 21.0 * z^2 - (4.0/3) * Ar

    # Choose bracket depending on Ar magnitude (helps convergence a lot)
    if Ar < 1e-2
        # Very small Ar → very small Re → small z
        bracket = (1e-8, 0.1)
    elseif Ar < 1e3
        bracket = (1e-6, 100.0)
    else
        # Typical to large particles
        bracket = (1e-5, 1e5)
    end

    local z  # Declare so it's visible after try

    try
        # Use a robust method; Bisection is very safe here
        z = find_zero(f, bracket, Roots.Bisection(), rtol=1e-10, atol=1e-12)
        
        if z ≤ 0
            error("Solver returned non-positive z = $z (invalid)")
        end
    catch e
        println("Root-finding failed for Ar = $Ar (dp = $dp m, T = $T K)")
        println("Bracket was: $bracket")
        println("Error was: ", e)
        # Optional fallback: approximate for intermediate regime
        # z ≈ sqrt( (4/3 * Ar) / 0.28 )   # Newton regime
        rethrow(e)  # or return NaN / approximate value
    end

    Re   = z^2
    u_elu = Re * μ / (ρg * dp)

    # Sanity checks
    if u_elu ≤ 0 || isnan(u_elu)
        @warn "Invalid velocity u_elu = $u_elu m/s (Ar = $Ar, z = $z)"
    end

    return u_elu
end

uelu = u_elu_Ar.(com_dp, 4000, T_com)

uelu_2 = u_elu_2.(com_dp, 4000, T_com)

# Reynolds
function Re_gas(dp, T, u0)
    """
    minimum fluidization velocity with: dp - m, density_solid/gas - kg/m^3, viscosity_gas - m^2 / s, T - K
    Note it is kinetmatic viscosity of gas
    """
    p = 1e5
    air_model   = PR(["Nitrogen","Oxygen"];idealmodel = WalkerIdeal)
    z_air = [0.785, 0.215]  # mole fractions

    density_gas = mass_density(air_model,p,T,z_air) # kg/m3
    dynamic_viscosity_gas = 1.458 * 10 ^ (-6) * T ^ 1.5 * (1/ (T + 110.4)) # Pa s
    viscosity_gas = dynamic_viscosity_gas/density_gas

    Re = dp * u0 / viscosity_gas
    return Re
end

Re_g = Re_gas.(com_dp, T_com, com_u0)