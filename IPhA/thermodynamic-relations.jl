include("IPhA/physconst.jl")



function isothermal_sound_speed_ideal_gas(
    T::Float64, # Temperature in Kelvin
    mu::Float64; # Mean molecular weight (dimensionless)
    mParticle::Float64=mProton, # Mean particle mass in grams (default is proton mass from physconst.jl)
    kB::Float64=kBoltzmann, # Boltzmann constant in erg/K (default is kBoltzmann from physconst.jl)
    )
    """
    Computes the isothermal sound speed for an ideal gas at temperature T and mean molecular weight mu.
    """
    
    return sqrt(kB * T / (mu * mParticle))
end