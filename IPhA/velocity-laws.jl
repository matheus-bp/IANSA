function beta_law(
    r::Float64,
    vinf::Float64, 
    beta::Float64; 
    vcore::Float64=0.0, 
    )
    # Calculate the velocity at a given radius using the beta law
    return vcore + (vinf - vcore) * (1 - (1 / r)^(beta)) 
end