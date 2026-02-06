using DifferentialEquations

### To convert liquid feed and concentration

using Plots, Clapeyron

function U_concentration_to_molar(U_feed_conc)
    """
    converts g(U)/L of UNH to molar of UNH and molar of H2O
    """
    water_model = SAFTgammaMie(["Water"])
 
    UN_MW  = 238.03 + 16*2 + 14*2 + 16*6       

    U_MW   = 238.03                 # g/mol   
    p_bed  = 1e5
    t_nozzle = 80 + 273.15

    water_mass_density  = mass_density(water_model, p_bed, t_nozzle) # kg/m3

    UN_feed_concentration = U_feed_conc / U_MW  * UN_MW # g(UNH)/L = kg/m3 

    # 1 / UNH feed concentration = mass or mole fraction UNH / density UNH + mass or molefraction water / density water 
    # at a certain temperature
    density_UN = 2.807 * 1000             # kg/m3, assume constant

    mass_x_H2O = water_mass_density * (density_UN - UN_feed_concentration) / (UN_feed_concentration *( UN_feed_concentration - water_mass_density)) # mass fraction of water
    mole_x_H2O = mass_x_H2O / 18.015 / (mass_x_H2O/18.015 + (1 - mass_x_H2O) / UN_MW)
    return mole_x_H2O
end

function U_concentration_to_mole_fraction(U_feed_conc)
    """
    Converts g(U)/L to mole fractions of UNH and H2O
    assuming ideal volume mixing
    """

    water_model = SAFTgammaMie(["Water"])

    UN_MW = 394.03        # g/mol
    U_MW   = 238.03          # g/mol
    H2O_MW = 18.015          # g/mol

    p_bed    = 1e5
    t_nozzle = 80 + 273.15

    ρ_H2O = mass_density(water_model, p_bed, t_nozzle)   # kg/m³
    ρ_UN = 2.807e3                                      # kg/m³

    # g(U)/L → kg(UNH)/m³
    c_UN = U_feed_conc / U_MW * UN_MW                  # kg/m³

    # Ideal volume mixing density
    ρ_mix = ρ_H2O * (1 - c_UN / ρ_UN) + c_UN

    # Mass fractions
    w_UN = c_UN / ρ_mix
    w_H2O = 1.0 - w_UN

    # Mole fractions
    x_H2O = (w_H2O / H2O_MW) /
            (w_H2O / H2O_MW + w_UN / UN_MW)


    return x_H2O
end



U_concentration_to_mole_fraction(1000)


U_concentration_to_molar(1000)

N = 100

conc_range = LinRange(900,1600,N)

U_conc_range = U_concentration_to_molar.(conc_range)

plot(conc_range, U_conc_range)




######## Assumed values

initial_dp = 200e-6     # particle diameter, m

t_bed = 350 + 273.15    # K

p_bed = 1e5             # Pa

air_to_liquid_ratio = 485

U = 265                         # W/m2/K

U_feed_concentration = 1000       # g(U)/L

t_air = t_bed                   # assume air preheated to bed temperature

t_nozzle = 80 + 273.15 # will be slightly less than evaporator and less than 100C

t_jacket = 100 + t_bed

fluidization_velocity = 2 # m/s


###### convert feed_concetration to UNH and water flowrates
# UNH flowrate is set value
UNH_molar_flowrate = 0.810411   # mol/s

UNH_MW = 502.13                 # g/mol          238.03 + 16*(8 + 6) + 14 * 2 + 12
U_MW   = 238.03                  # g/mol   

water_mass_density  = mass_density(water_model, p_bed, t_nozzle) # kg/m3
water_molar_density = molar_density(water_model, p_bed, t_nozzle) # mol/m3

UNH_mass_flowrate = UNH_molar_flowrate * UNH_MW # g/s

UNH_feed_concentration = U_feed_concentration / U_MW  * UNH_MW # g(UNH)/L = kg/m3 

# 1 / UNH feed concentration = mass or mole fraction UNH / density UNH + mass or molefraction water / density water 
# at a certain temperature
density_UNH = 2.807 * 1000             # kg/m3, assume constant

mass_x_H2O = water_mass_density * (density_UNH - UNH_feed_concentration) / (UNH_feed_concentration *( UNH_feed_concentration - water_mass_density)) # mass fraction of water
mole_x_H2O = mass_x_H2O / 18.015 / (mass_x_H2O/18.015 + (1 - mass_x_H2O) / UNH_MW)
mole_x_UNH = 1- mole_x_H2O

###### water flowrate 
H2O_molar_flowrate = UNH_molar_flowrate * mole_x_H2O / ( 1 - mole_x_H2O) # mol/s

H2O_density        = molar_density(water_model,p_bed, t_nozzle)

air_density        = molar_density(air_model, p_bed, t_nozzle, z_air) # mol/m3

volume_nozzle_air_flowrate = H2O_molar_flowrate / H2O_density * air_to_liquid_ratio    # m3/s

air_nozzle_molar_flowrate        = volume_nozzle_air_flowrate * air_density           # mol/s



function simplified_dp_ode!(ddp, dp, p, t)
    # unpack parameters
    M_bed, UO3_mass_rate = p

   # total_area = M_bed / solid_density * (6 / dp)

    # ODE: d(dp)/dt
    ddp[1] = (2 * UO3_mass_rate / (6 * M_bed )) * dp[1]
end

function simplified_dp_ode_2!(ddp, dp, p, t)
    # unpack parameters
    M_bed, UO3_mass_rate = p

   # total_area = M_bed / solid_density * (6 / dp)
   #total_area = 8 *  M_bed / solid_density * (d_p)
   dp_0 = 200e-6 * 0.999

    # ODE: d(dp)/dt
    ddp[1] = (2 * UO3_mass_rate / (8 * M_bed * (dp[1] ^ 3 - dp_0 ^ 3) / (dp[1] ^ 4 - dp_0 ^ 4) ))
end


UNH_molar_flowrate = 0.810411        # mol/s
MW_UO3 = 238 + 16*3                  # kg/kmol
UO3_mass_rate = UNH_molar_flowrate * MW_UO3   # g/s

M_bed = 10000.0                       # kg
dp_0  = 200e-6                       # m

p = (M_bed, UO3_mass_rate)
dp0 = [dp_0]

tspan = (0.0, 3600)                  # seconds

prob = ODEProblem(simplified_dp_ode_2!, dp0, tspan, p)
sol = solve(prob)
dp_final = sol(60.0)[1] * 10^6
println("Final particle diameter = $(dp_final) m")

mass = pi * (0.666/2) ^ 2 * 4.8 * 7300