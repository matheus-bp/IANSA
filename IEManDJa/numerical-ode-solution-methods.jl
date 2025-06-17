function fw_euler_step(
    x::Vector{Float64},
    y::Vector{Float64},
    h::Float64,
    dydx::Float64   
    )
"""
    Makes one step solving the ordinary differential equation (ODE) "dydx = f(y,x)" using the Forward Euler (FE) algorithm.
"""

    return ynew
end


# RK4 implementation
function rk4_step(
    dxdy, 
    x, 
    y, 
    h
    )
    """
    Makes one step solving the ordinary differential equation (ODE) "dy/dx = f(x, y)" using the Runge-Kutta 4th order method.
    """
    k1 = h * dxdy(x, y)
    k2 = h * dxdy(x + h/2, y + k1/2)
    k3 = h * dxdy(x + h/2, y + k2/2)
    k4 = h * dxdy(x + h, y + k3)
    
    return y + (k1 + 2k2 + 2k3 + k4)/6
end


###########################################################

# The solution of the isothermal wind equation is essentially the solution of an ODE. However, the dv/dr is not easily isolatable.
# Therefore One needs to first modify the equation so that the ODE can be solved numerically (e.g. using Euler, RK4, etc).
# 1. Modify variables: z = r^2 * v
# 2. Find the root of the EOM (e.g., using Newton-Raphson, Brent, Secant, etc)
# 3. Recover the dv/dr from the root (z_root) of the EOM
# 4. Solve the ODE using a numerical method (e.g., RK4, Euler, etc)

function F(
    z, #z = r^2 v dvdr
    r, # radial distance
    v, # velocity
    Mstar, # stellar mass
    aSound, # isothermal sound speed
    alpha :: Float64; # stellar mass
    G :: Float64 = physconst.GGrav, # gravitational constant
    )
    GMeff = G * Mstar * (1 - Gamma_e)
    (1 - aSound^2/v^2)*z + (GMeff - 2*aSound*r) - C*z^(1/alpha)

end

# Function to find the root (z_root) of the equation F(z) = 0
function z_root(r, v, fz)
    # Find bracket (adjust these bounds if needed)
    z_min = 0.0
    z_max = 100.0
    try
        return brent(f, z_min, z_max)
    catch e
        # Try with different bounds if first attempt fails
        z_max = 1000.0
        return brent(f, z_min, z_max)
    end
end


# ODE function
function ode_func(x, y)
    z = solve_z(x, y)
    return z / (x^2 * y)
end

# Main solver function
function solve_ode(x0, y0, x_end, h)
    x = x0
    y = y0
    xs = [x]
    ys = [y]
    
    while x < x_end
        y = rk4_step(ode_func, x, y, h)
        x += h
        push!(xs, x)
        push!(ys, y)
    end
    
    return xs, ys
end