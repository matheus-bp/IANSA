function blackbody_lambda(
    λ::Float64, 
    T::Float64;
    c::Float64=physconst.cLight,
    h::Float64=physconst.hPlanck,
    kB::Float64=physconst.kBoltzmann
    )
    """
    Computes the blackbody intensity (Planck's formula) at a certain wavelength λ for a temperature T
    """
    numerator = 2 * h * c^2 / λ^5
    exponent = (h * c) / (λ * kB * T)
    denominator = exp(exponent) - 1
    
    return numerator / denominator
end