# https://ui.adsabs.harvard.edu/abs/1970MNRAS.149..111C/abstract

include("IPhA/physconst.jl")


function Nl(
    T::Real,
    N_e::Real, 
    N_i::Real, 
    gl::Real, 
    bl::Real, 
    chil::Real; 
    m_p::Real=physconst.mProton, 
    kB::Real=physconst.kBoltzmann, 
    h::Real=physconst.hPlanck 
    )
    """
    Computes the number density of the lower level of a two-level atom in thermal equilibrium.
    """    
    return gl * (h^3 / (2*(2*pi*m_p*kB*T)^(3/2))) * N_e*N_i*bl*exp(chil/kB*T)
end

function tau_0(r, v, Nl, )
     return ((pi * elec^2) / (m_p * c_cgs)) * gf * ( ((Nl/gl) - (Nu/gu))/(nu_0 * v / c_cgs) ) * r
end

function planck_lambda(λs::AbstractArray{<:Real}, T::Real)
    # Fundamental physical constants
    hc = h * c
    hc2 = h * c^2
    kBT = k * T
    
    B = similar(λs)
    @inbounds @simd for (i, λ) in enumerate(λs)
    # @simd: Enables SIMD (Single Instruction Multiple Data) vectorization
    # @inbounds: Disables array bounds checking
        exponent = hc / (λ * kBT)
        B[i] = (2 * hc2 / λ^5) / (exp(exponent) - 1)
    end
    return B
end