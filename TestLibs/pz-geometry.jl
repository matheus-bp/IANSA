#######################
if abspath(PROGRAM_FILE) == @__FILE__

Pkg.add("CairoMakie")
using Pkg
using CairoMakie
const mak = CairoMakie  

# create a grid with the pz geometry scheme

ND :: Int64 = 50 # number of depthpoints, i.e., the number of atmospheric impact parameters
NC :: Int64 = 7 # number of core rays

r_min :: Float64 = 1.0 # in stellar radii
r_max :: Float64 = 10.0 # also in stellar radii

radius = zeros(Float64, ND) # initialize r_grid
imppar_atm = zeros(Float64, ND-1) # initialize atmospheric p grid
imppar_cor = zeros(Float64, NC) # initialize core p grid

# radius grid spacing
radius = collect(range(r_max, r_min, length=ND))

# impact parameters 
imppar_atm = radius[1:end-1]
imppar_cor = collect(range(r_min, 0, length=NC)) # assumes a linear spacing of the core rays

zrays_atm = [Vector{Float64}() for _ in 1:ND-1]
zrays_aba =  [Vector{Float64}() for _ in 1:ND-1] # aba := Atmosphere BAck
zrays_cor = [Vector{Float64}() for _ in 1:NC]

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

colors_atm = mak.cgrad(:Blues_3, ND-1) 
colors_aba = mak.cgrad(:Greens_3, ND-2) 
colors_cor = mak.cgrad(:BuPu_3, NC) 

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


for p in collect(1:ND-1)
    imppar = fill(imppar_atm[p], length(zrays_atm[p]))
    zaxis = zrays_atm[p]
    mak.scatter!(ax, zaxis, imppar[end:-1:1], color=colors_atm[p])
end

for p in collect(1:ND-1)
    imppar = fill(imppar_atm[p], length(zrays_aba[p]))
    zaxis = zrays_aba[p]
    mak.scatter!(ax, zaxis, imppar[end:-1:1], color=colors_aba[p])
end

for p in collect(1:NC)
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
function pz_geometry_rays(radius, NC)

# impact parameters 
imppar_atm = radius[1:end-1]
imppar_cor = collect(range(r_min, 0, length=NC)) # assumes a linear spacing of the core rays

zrays_atm = [Vector{Float64}() for _ in 1:ND-1]
zrays_aba =  [Vector{Float64}() for _ in 1:ND-1] # aba := Atmosphere BAck
zrays_cor = [Vector{Float64}() for _ in 1:NC]

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

return imppar_atm, imppar_cor, zrays_atm, zrays_aba, zrays_cor

end











end