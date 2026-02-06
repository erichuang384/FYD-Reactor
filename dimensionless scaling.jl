# Pilot plant data
# https://www.osti.gov/servlets/purl/4193860
dia_pilot = 10 * 2.54 * 10 ^ (-2)   # m pilot plant diameter
len_pilot = 6 * 0.3048              # m pilot reactor length
nozzle_location_pilot = 3 * 0.3048 + 6 * 0.0254 # was 3'6'', in m

T_pilot   = (595-32) * 5/9  + 273.15        # K pilot temperature
T_com = T_pilot
pilot_prod_rate = 140 * 0.45359237        # kg/hr UO3
pilot_prod_rate_area = pilot_prod_rate/ (pi * (dia_pilot/2) ^ 2) # kg/(hr*m^2)

fluidizing_air_flow = 20 * 1.699    # m3/hr air, it is in standard CFM need to check conversion

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

com_dia         = com_prod_rate/pilot_prod_rate_area    

####### Apply scaling laws
# bed height/bed diameter
com_len = com_dia * len_pilot / dia_pilot

com_nozzle = com_len/len_pilot * nozzle_location_pilot #  m nozzle location of full scale

# inertia / gravity forces for superficial gas velocity
#""The average superficial velocity of the fluidizing gas in the heater zone is 1.7 ft/sec and drop to a value of 1.25 ft/eec in the spray zone.""

pilot_u0 = 1.7 * 0.3048 # m/s in AUSTRALIAN ATOMIC ENERGY COMMISSION also has u0 of 0.5 m/s

com_u0   = sqrt(pilot_u0 ^ 2 / dia_pilot * com_dia) # m/s

# particle inertia/gas viscous u0 * dp
# assume particle diameter of 0.01 inches, graph/table hard to figure out
pil_dp = 0.01 * 0.0254  # m
com_dp = pil_dp * pilot_u0 / com_u0 # within graph


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
umf = u_mf(com_dp, 4000,T_com)
println(com_u0>umf)
function u_elu(dp, density_solid, T)
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
    Re_elu = 1.74 * Ar / (31.3 + sqrt(Ar))
    uelu   = Re_elu * viscosity_gas / dp

    return uelu
end
uelu = u_elu(com_dp, 4000, T_com)
uelu2 = u_elu_2(com_dp, 4000, T_com)

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

Re_g = Re_gas(com_dp, T_com, com_u0)