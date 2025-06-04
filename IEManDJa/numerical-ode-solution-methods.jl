function solveODE_FE(
    y::Vector{Float64},
    h::Float64,
    dydx::Float64
)
"""
    Solves the ordinary differential equation (ODE) "dydx = f(y,x)" using the Forward Euler (FE) algorithm.
"""

    return ynew
end