using DifferentialEquations, Plots, Clapeyron

# Sizing for a 99.5% conversion

# Air equations of state for CAg
function air_CAg(p,T)
    """
    Input temperature in Celsius and pressure in Pa
    """
    eos_air = PR(["Nitrogen","Oxygen"])
    x_air = [0.785, 0.215]  # mole fractions
    air_molar_density = molar_density(eos_air, p, T+273.15, x_air) * 0.215# * 10^6
    return air_molar_density
end

# 3UO2 + O2 -> U3O8  (JEONG, 2006). BASED ON SHRINKING CORE MODEL

# Table 1 data
params = Dict(
    336 => (kg = 7.33,  De = 7.407, ks = 2.47),
    345 => (kg = 7.46,  De = 7.457, ks = 5.24),
    377 => (kg = 10.48, De = 7.576, ks = 9.37),
    395 => (kg = 15.17, De = 7.675, ks = 15.63),
    440 => (kg = 19.01, De = 7.752, ks = 25.38),
    454 => (kg = 21.99, De = 7.825, ks = 30.40)
)

function t_SCM_model(x; T, b, n, ρs, r, CAg)

    p = params[T]
    kg, De, ks = p.kg, p.De, p.ks

    term1 = (ρs * r) / (2 * b * kg * CAg) * x
    term2 = (ρs * r^2) / (4 * b * De * CAg) *
            (x + (1 - x) * log(1 - x))
    term3 = (ρs * r) / (b * ks * CAg) *
            (-log(1 - x))^(1 / n)

    return term1 + term2 + term3
end

# constants
b   = 3.0
n   = 4.0
ρs  = (10.97/270.03)/10^-6 # mol/m3, is a function of T,p but assumed constant here since it's a solid. g/cm3/(g/mol) 
r   = 7.14e-3* 10^2 # cm 
r = 500e-6 * 10^2 / 2

# Choose temperature
T = 454
p = 1e5 + 45e3

# use PR for this
CAg = air_CAg(p,T)

xi = 0.99
tval = t_SCM_model(xi; T=T, b=b, n=n, ρs=ρs, r=r, CAg=CAg)

###### Weight of bed
# Requires inlet flowrate
U_MW        = 238.02891 # g/mol
U3O8_MW     = U_MW * 3 + 15.9994 * 8 # g/mol
UO2_MW      = U_MW + 15.9994 * 2 # g/mol  

U_prod      = 5000 * 10^3 / (1.0 * 365 * 24 *3600)      # kg/s, replace operating days
U_prod_mol  = U_prod / U_MW                             # kmol/s
U3O8_prod   = U_prod_mol / 3                            # kmol/s

F_in_UO2 = U3O8_prod * 3 / xi                           # kmol/s

# Residence time
#t_res =  ρs * r / (2 * b * params[T].kg * CAg) + ρs * r ^ 2 / (4 * b * params[T].De * CAg) +  ρs * r / (b * params[T].ks * CAg)

# Bed weight
bed_weight = F_in_UO2 * UO2_MW * tval * 60 # kg/s * s

# Volume
bed_volume = bed_weight / ((10.97 * 1000)*(1-0.4)) # m3
