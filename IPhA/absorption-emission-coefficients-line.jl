include("physconst.jl")

function emission_coefficien(transition_frequency_ij, occupation_number_i, einstein_A_ij, line_profile; planck_h = planck_h)
    return ( (planck_h * transition_frequency_ij) / 4*pi ) * occupation_number_i * einstein_A_ij * line_profile
end

function absorption_coefficient(transition_frequency_ij, occupation_number_i, occupation_number_j, einstein_B_ji, einstein_B_ij, line_profile; planck_h = planck_h)
    return ( (planck_h * transition_frequency_ij) / 4*pi ) * 
    (occupation_number_j * einstein_B_ji - occupation_number_i * einstein_B_ij) * line_profile
end