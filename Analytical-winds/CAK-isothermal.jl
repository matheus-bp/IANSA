include("IEManDJa/root-finding-methods.jl")
include("IEManDJa/numerical-ode-solution-methods.jl")
include("IPhA/physconst.jl")   

# Following the Chapter 8 of Stellar Winds

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
    gamma_e   :: Float64,
    G_grav    :: Float64,
    )
"""
    Computes the velocity at the critical point based on the conditions of an isothermal line-driven wind with non-negligible pressure term.
"""

    return ( a_sound^2 + (1/(1-alpha_CAK)) * ( (G_grav * Mstar * (1-gamma_e) / r_crit) - 2*c_sound^2 ) )^0.5

end

########################



###########################################################

# The solution of the isothermal wind equation is essentially the solution of an ODE. However, the dv/dr is not easily isolatable.
# Therefore One needs to first modify the equation so that the ODE can be solved numerically (e.g. using Euler, RK4, etc).
# 1. Modify variables: z = r^2 * v
# 2. Find the root of the EOM (e.g., using Newton-Raphson, Brent, Secant, etc)
# 3. Recover the dv/dr from the root (z_root) of the EOM
# 4. Solve the ODE using a numerical method (e.g., RK4, Euler, etc)

function F(
    r, # radial distance
    v, # velocity
    Mstar, # stellar mass
    aSound, # isothermal sound speed
    alpha :: Float64; # alpha exponent of CAK Theory
    G :: Float64 = GGrav, # gravitational constant
    )
    z = r^2 * v # modified variable
    GMeff = G * Mstar * (1 - Gamma_e)
    (1 - aSound^2/v^2)*z + (GMeff - 2*aSound*r) - C*z^(1/alpha)

end

fz = F(r,v, 30, a_sound, 0.3)

# Function to find the root (z_root) of the equation F(z) = 0
function z_root(r, v, fz)
    # Find bracket (adjust these bounds if needed)
    z_min = 0.0
    z_max = 100.0

    try
        return brent(fz, z_min, z_max)
    catch e
        # Try with different bounds if first attempt fails
        z_max = 1000.0
        return brent(fz, z_min, z_max)
    end
end


# ODE function
function ode_input_function(r, v, F)
    z = z_root(r, v, F)
    return z / (r^2 * v)
end

# Main solver function
solve_ode(ode_input_function, x0, y0, x_end, h; chosen_method=rk4_step)


#################

NN_max = 10

while NN < NN_max
    
    

end