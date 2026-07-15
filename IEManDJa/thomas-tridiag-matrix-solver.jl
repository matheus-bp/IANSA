#######################
if abspath(PROGRAM_FILE) == @__FILE__

using Pkg
Pkg.add("CairoMakie")
using CairoMakie
const mak = CairoMakie  

########################
else
########################
"""
    thomas(A, B, C, S)

Solve

    A[i]u[i-1] + B[i]u[i] + C[i]u[i+1] = S[i]

using the Thomas algorithm.

A[1] and C[end] are ignored.
"""
function thomas(A, B, C, S)
    n = length(B)

    cp = similar(C)
    dp = similar(S)

    # Forward sweep
    cp[1] = C[1] / B[1]
    dp[1] = S[1] / B[1]

    for i in 2:n
        denom = B[i] - A[i] * cp[i-1]
        cp[i] = i < n ? C[i] / denom : zero(eltype(C))
        dp[i] = (S[i] - A[i] * dp[i-1]) / denom
    end

    # Back substitution
    u = similar(S)
    u[n] = dp[n]

    for i in (n-1):-1:1
        u[i] = dp[i] - cp[i] * u[i+1]
    end

    return u
end



########################
end
########################