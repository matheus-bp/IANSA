
function W_geom(
    Rstar :: Float64,
    r     :: Float64
)
"""
    Function that computes the geometrical dillution factor, which represents the solid angle covered by the star seen from the position "r"
"""    
    
    return (1/2) * ( 1 - (1-(Rstar/r)^2)^0.5 )

end

########################

function v_crit(
    r_crit    :: Float64,
    Mstar     :: Float64, 
    c_sound   :: Float64, 
    alpha_CAK :: Float64, 
    gamma_e   :: Float64; 
    G_grav    :: Float64
)
"""
    Computes the velocity at the critical point based on the conditions of an isothermal line-driven wind with non-negligible pressure term.
"""

    return ( a_sound^2 + (1/(1-alpha_CAK)) * ( (G_grav * Mstar * (1-gamma_e) / r_crit) - 2*c_sound^2 ) )^0.5

end

