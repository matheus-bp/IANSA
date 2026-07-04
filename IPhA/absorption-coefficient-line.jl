include("physconst.jl")


function absorption_coefficient(transition_frequency, occupation_number_i, occupation_number_j, einstein_B_ji, einstein_B_ij, line_profile; planck_h = planck_h)
    return ( (planck_h * transition_frequency) / 4*pi ) * 
    (occupation_number_j * einstein_B_ji - occupation_number_i * einstein_B_ij) * line_profile
end