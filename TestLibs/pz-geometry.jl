# create a grid with the pz geometry scheme

ND :: Int64 = 20 # number of depthpoints, i.e., the number of atmospheric impact parameters
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
    for i in collect(Int64, 1:j)
        z = (radius[j]^2 - imppar_atm[i]^2)^0.5
        push!(zaxes_atm[j], z)
    end
end

println(zaxes_atm)

