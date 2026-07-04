include("physconst.jl")

function line_profile_doppler(frequency, projected_velocity; c = c_light)
    return line_profile(frequency * (1 - (1/c) * projected_velocity))
end