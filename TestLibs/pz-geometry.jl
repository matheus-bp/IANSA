# create a grid with the pz geometry scheme

ND :: Int64 = 5 # number of depthpoints, i.e., the number of atmospheric impact parameters
NC :: Int64 = 5 # number of core rays

r_min :: Float64 = 1.0 # in stellar radii
r_max :: Float64 = 100.0 # also in stellar radii

radius = zeros(Float64, ND) # initialize r_grid
imppar_atm = zeros(Float64, ND-1) # initialize atmospheric p grid
imppar_cor = zeros(Float64, NC) # initialize core p grid

# Grid spacing

radius[1] = r_min

# linear spacing
for (i, r) in enumerate(radius[2:end])
    radius[i+1] += radius[1] + i*(r_max - r_min)/ND
end

imppar_atm = radius[1:end-1]

zaxes_atm = [Vector{Float64}() for _ in 1:ND-1]
zaxes_cor = [Vector{Float64}() for _ in 1:NC]

for (j,p) in enumerate(imppar_atm)
    for i in collect(Int64, j:-1:1)
        z = (radius[j]^2 - imppar_atm[i]^2)^0.5
        push!(zaxes_atm[j], z)
    end
end

########################

if abspath(PROGRAM_FILE) == @__FILE__

try
    using CairoMakie
catch
    using Pkg
    Pkg.add("CairoMakie")
    using CairoMakie
end
const mak = CairoMakie  # Now mak is defined in global scope
    
fig = mak.Figure(size = (650,650))
ax = mak.Axis(fig[1,1])

for p in collect(1:ND-1)
    imppar = fill(imppar_atm[p], length(zaxes_atm[p]))
    zaxis = zaxes_atm[p]
    println(imppar)
    mak.scatter!(ax, zaxis, imppar[end:-1:1])
end

for r in radius
    mak.arc!(ax, 0.0, r, 0, π/2, color=:blue, linewidth=0.5)
end

display(fig)
println("Press Enter to close...")
readline()  # Wait for user to press Enter

########################
else
    println("imported pz-geometry.jl")
end