# https://ui.adsabs.harvard.edu/abs/1970MNRAS.149..111C/abstract

include("IPhA/physconst.jl")


function Nl(
    T::Float64,
    N_e::Float64, 
    N_i::Float64, 
    gl::Float64, 
    bl::Float64, 
    chil::Float64; 
    m_p::Float64=physconst.mProton, 
    kB::Float64=physconst.kBoltzmann, 
    h::Float64=physconst.hPlanck 
    )
    """
    Computes the number density of the lower level of a two-level atom in thermal equilibrium.
    """    
    return gl * (h^3 / (2*(2*pi*m_p*kB*T)^(3/2))) * N_e*N_i*bl*exp(chil/kB*T)
end

###########################################

function tau_0(
    r::Float64,
    v::Float64,
    nu0::Float64,
    Nl::Float64,
    Nu::Float64;
    gl::Float64=1.0,
    gu::Float64=1.0,
    gf::Float64=1.0,
    c::Float64=physconst.cLight,
    m_p::Float64=physconst.mProton,
    elec::Float64=physconst.eCharge,
    )

    """
    Computes the optical depth at the line center for a two-level atom in a stellar wind. 
    """

     return ((pi * elec^2) / (m_p * c)) * gf * ( ((Nl/gl) - (Nu/gu))/(nu0 * v / c) ) * r
end

############################################

# function planck_lambda_vector(λs::AbstractArray{<:Float64}, T::Float64)
#     # Fundamental physical constants
#     hc = h * c
#     hc2 = h * c^2
#     kBT = k * T
    
#     B = similar(λs)
#     @inbounds @simd for (i, λ) in enumerate(λs)
#     # @simd: Enables SIMD (Single Instruction Multiple Data) vectorization
#     # @inbounds: Disables array bounds checking
#         exponent = hc / (λ * kBT)
#         B[i] = (2 * hc2 / λ^5) / (exp(exponent) - 1)
#     end
#     return B
# end

#################################################

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