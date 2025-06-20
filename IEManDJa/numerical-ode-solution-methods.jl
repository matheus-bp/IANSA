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

