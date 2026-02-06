using Plots, Clapeyron

water_model = SAFTgammaMie(["Water"];idealmodel = WalkerIdeal)
air_model   = PR(["Nitrogen","Oxygen"];idealmodel = WalkerIdeal)
z_air = [0.785, 0.215]  # mole fractions

water_mass_density  = mass_density(water_model,1e5, 80 + 273.15) # kg/m3
water_molar_density = molar_density(water_model,1e5, 80 + 273.15) # mol/m3


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

Heat_capacity_water = isobaric_heat_capacity(water_model, p_bed, t_nozzle)      # J/mol/K

Heat_capacity_UNH = 120 * 4.184# J/mol/K THIS IS VALUE AT 300 K https://pubs.acs.org/doi/pdf/10.1021/ja01867a061?ref=article_openPDF

Heat_capacity_feed = Heat_capacity_UNH * mole_x_UNH + Heat_capacity_water * mole_x_H2O

Heat_capacity_air = isobaric_heat_capacity(air_model,p_bed,t_air,z_air)         # J/mol/K

##### Water enthalpy of vaporization 

enthalpy_vap_water = enthalpy_vap(water_model, t_bed)   # J/mol

enthalpy_reaction_UNH_to_UO3 = 136 * 4.184 * 1000 # J/mol, this is value at 300C https://onlinelibrary.wiley.com/doi/epdf/10.1002/jctb.5010140907?saml_referrer=

area = (UNH_molar_flowrate * enthalpy_reaction_UNH_to_UO3 + H2O_molar_flowrate * enthalpy_vap_water 
- (UNH_molar_flowrate * Heat_capacity_UNH + H2O_molar_flowrate * Heat_capacity_water + air_nozzle_molar_flowrate * Heat_capacity_air) * (t_nozzle - t_bed)) / (U * (t_jacket - t_bed)) # mol/s * J/mol    +   

Heating_duty = area * U * (t_jacket-t_bed)

diameter = 1.2 # m

length = area/(pi * diameter)

volume_fbr = (diameter/2)^2 * pi * length

volume_fbr_design = (diameter/2)^2 * pi * 10


# Fluidizing air required

volume_fluidizing_air = fluidization_velocity * pi * (diameter/2) ^ 2