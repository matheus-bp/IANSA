function fw_euler_step(
    y::Vector{Float64},
    h::Float64,
    dydx::Float64
)
"""
    Solves the ordinary differential equation (ODE) "dydx = f(y,x)" using the Forward Euler (FE) algorithm.
"""

    return ynew
end


# RK4 implementation
function rk4_step(f, x, y, h)
    k1 = h * f(x, y)
    k2 = h * f(x + h/2, y + k1/2)
    k3 = h * f(x + h/2, y + k2/2)
    k4 = h * f(x + h, y + k3)
    return y + (k1 + 2k2 + 2k3 + k4)/6
end