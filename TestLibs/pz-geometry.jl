#######################
if abspath(PROGRAM_FILE) == @__FILE__

using Pkg
Pkg.add("CairoMakie")
using CairoMakie
const mak = CairoMakie  

# create a grid with the pz geometry scheme

ND :: Int64 = 20 # number of depthpoints, i.e., the number of atmospheric impact parameters
NC :: Int64 = 5 # number of core rays

r_min :: Float64 = 1.0 # in stellar radii
core_p_max_factor = 0.8  # sets the maximum impact parameter for core rays
r_max :: Float64 = 10.0 # also in stellar radii

radius = zeros(Float64, ND) # initialize r_grid
imppar_atm = zeros(Float64, ND-0) # initialize atmospheric p grid
imppar_cor = zeros(Float64, NC-1) # initialize core p grid

# radius grid spacing
radius = collect(range(r_max, r_min, length=ND))

# impact parameters 
imppar_atm = radius[1:end-0]
imppar_cor = collect(range(core_p_max_factor * r_min, 0, length=NC-1+1)) # assumes a linear spacing of the core rays

zrays_atm = [Vector{Float64}() for _ in 1:ND-0]
zrays_aba =  [Vector{Float64}() for _ in 1:ND-0] # aba := Atmosphere BAck
zrays_cor = [Vector{Float64}() for _ in 1:NC-1+1]

for (j,p) in enumerate(imppar_atm)
    for i in collect(Int64, 1:j)
        z = (radius[i]^2 - p^2)^0.5
        push!(zrays_atm[j], z)
    end
end

for (j,p) in enumerate(imppar_atm)
    for i in collect(Int64, 1:j)
        z = -(radius[i]^2 - p^2)^0.5
        push!(zrays_aba[j], z)
    end
end

for z in zrays_aba
    pop!(z)
end

for (j,p) in enumerate(imppar_cor)
    for (i,r) in enumerate(radius[end:-1:1])
        z = (r^2 - p^2)^0.5
        push!(zrays_cor[j], z)
    end
end

    
fig = mak.Figure(size = (800,400))
ax = mak.Axis(fig[1,1], xgridvisible=false, ygridvisible=false,)

colors_atm = mak.cgrad(:Blues_3, ND-0) 
colors_aba = mak.cgrad(:Greens_3, ND-0) 
colors_cor = mak.cgrad(:BuPu_3, NC-1+1) 

for r in radius
    mak.arc!(ax, 0.0, r, 0, π, color=:gray, linewidth=0.5,alpha=0.5)
end

for p in imppar_atm
    mak.lines!(ax,[-r_max,r_max],[p,p], linewidth=2, alpha=0.2, color=:black)
end

for p in imppar_cor
    z = (r_min^2 - p^2)^0.5
    mak.lines!(ax,[z,r_max],[p,p], linewidth=2, alpha=0.2, color=:black)
end


for p in collect(1:ND-0)
    imppar = fill(imppar_atm[p], length(zrays_atm[p]))
    zaxis = zrays_atm[p]
    mak.scatter!(ax, zaxis, imppar[end:-1:1], color=colors_atm[p])
end

for p in collect(1:ND-0)
    imppar = fill(imppar_atm[p], length(zrays_aba[p]))
    zaxis = zrays_aba[p]
    mak.scatter!(ax, zaxis, imppar[end:-1:1], color=colors_aba[p])
end

for p in collect(1:NC-1+1)
    imppar = fill(imppar_cor[p], length(zrays_cor[p]))
    zaxis = zrays_cor[p]
    mak.scatter!(ax, zaxis, imppar, color = colors_cor[p])
end


display(fig)
println("Press Enter to close...")
readline()  # Wait for user to press Enter

########################
else
########################

"""
    pz_geometry_rays(radius, NC)

TBW
"""
function pz_geometry_rays(radii_atm, radii_cor)

# impact parameters 
imppar_atm = radii_atm[1:end]
imppar_cor = radii_cor[1:end]

ND = length(imppar_atm)
NC = length(imppar_cor)

zrays_atm = [Vector{Float64}() for _ in 1:ND]
zrays_aba =  [Vector{Float64}() for _ in 1:ND] # aba := Atmosphere BAck
zrays_cor = [Vector{Float64}() for _ in 1:NC]

for (j,p) in enumerate(imppar_atm)
    for i in collect(Int64, 1:j)
        z = (radii_atm[i]^2 - p^2)^0.5
        push!(zrays_atm[j], z)
    end
end

for (j,p) in enumerate(imppar_atm)
    for i in collect(Int64, 1:j)
        z = -(radii_atm[i]^2 - p^2)^0.5
        push!(zrays_aba[j], z)
    end
end
for z in zrays_aba
    pop!(z)
end

for (j,p) in enumerate(imppar_cor)
    for (i,r) in enumerate(radii_atm[end:-1:1])
        z = (r^2 - p^2)^0.5
        push!(zrays_cor[j], z)
    end
end

return imppar_atm, imppar_cor, zrays_atm, zrays_aba, zrays_cor

end

########################
end
########################