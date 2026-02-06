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