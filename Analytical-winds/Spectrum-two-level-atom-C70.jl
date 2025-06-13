# https://ui.adsabs.harvard.edu/abs/1970MNRAS.149..111C/abstract


Nl = gl * (h^3 / (2*(2*pi*m_proton*k_B*T)^(3/2))) * N_e*N_i*bl*exp(chil/k_b*T)

tau_0(r, v, Nl, ) = ((pi * elec^2) / (m_p * c_cgs)) * gf * ( ((Nl/gl) - (Nu/gu))/(nu_0 * v / c_cgs) ) * r

function planck_lambda(λs::AbstractArray{<:Real}, T::Real)
    # Fundamental physical constants
    hc = h * c
    hc2 = h * c^2
    kBT = k * T
    
    B = similar(λs)
    @inbounds @simd for (i, λ) in enumerate(λs)
        exponent = hc / (λ * kBT)
        B[i] = (2 * hc2 / λ^5) / (exp(exponent) - 1)
    end
    return B
end