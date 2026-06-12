# three_level_radiative_transfer.jl
# Toy 1D plane-parallel radiative transfer + rate equations for a 3-level hydrogen-like atom
# Produces an emergent spectral line (flux vs wavelength)
# Requirements: Julia 1.6+ and Plots.jl
# Usage: julia three_level_radiative_transfer.jl

#import Pkg;
#Pkg.add("Plots")
#Pkg.add("LinearAlgebra")

using LinearAlgebra
using Plots

# ------------------------- Physical constants -------------------------
const c = 2.99792458e10      # cm/s
const h = 6.62607015e-27     # erg*s
const kB = 1.380649e-16      # erg/K
const eV = 1.602176634e-12   # erg

# ------------------------- Model parameters -------------------------
T = 8000.0                   # K (gas temperature)
ne = 1e12                   # electron density cm^-3 (sets collisional rates)
n_total = 1.0                # total number density of the toy atom (arbitrary units)

# Geometry
Nz = 60                      # number of depth points (plane-parallel slab)
depth = range(0.0, stop=1.0, length=Nz)   # geometrical coordinate (arbitrary units)

# Frequency / wavelength grid around a single line (1 <-> 2 transition)
lambda0 = 1215.67e-8         # cm; Lyman-alpha-ish (toy) -> 1215.67 Angstrom
nλ = 201
Δλ = 2e-9                    # cm ~ small window
lam = range(lambda0 - (nλ-1)/2*Δλ, stop=lambda0 + (nλ-1)/2*Δλ, length=nλ)
ν = c ./ lam

# Doppler width (thermal)
mass_H = 1.6735575e-24      # g
vth = sqrt(2*kB*T/mass_H)
Δν_D = (vth/c) * ν[div(nλ,2)+1]

# Line profile: Gaussian (complete redistribution assumed)
function phi_nu(nu, nu0, ΔnuD)
    x = (nu - nu0)/ΔnuD
    return exp(-x^2)/ (sqrt(pi)*ΔnuD)
end

# ------------------------- Atomic data (toy) -------------------------
# Levels: 1 (ground), 2 (excited), 3 (continuum/ionized)
E1 = 0.0
E2 = 10.2 * eV             # ~n=2 energy above ground (erg)
E3 = 13.6 * eV             # ionization energy (erg)
ν12 = (E2 - E1)/h
λ12 = c/ν12
A21 = 6.25e2               # s^-1 (toy; large to make line strong)
# Einstein B coefficients from A (using g-factors simplified)
g1 = 2.0; g2 = 8.0
B21 = A21 * (λ12^3) / (2*h*c)   # using relation A = (2hν^3/c^2) * B21 -> B21 = A c^2/(2hν^3)
B12 = (g2/g1) * B21

# Collisional rates (very simplified approximations)
function C12(T, ne)
    # collisional excitation (s^-1) per particle in level 1
    # simplified Arrhenius-type: C12 = ne * q12; use q12 = 1e-8 * exp(-ΔE/kT) (cm^3 s^-1)
    ΔE = E2 - E1
    q12 = 1e-8 * exp(-ΔE/(kB*T))
    return ne * q12
end
function C21(T, ne)
    # detailed balance: C21 = C12 * (g1/g2) * exp(ΔE/kT)
    ΔE = E2 - E1
    return C12(T, ne) * (g1/g2) * exp(ΔE/(kB*T))
end

# Collisional ionization/recombination toy rates between level 2 and continuum (3)
function C23(T, ne)
    # collisional ionization from 2 -> 3
    return ne * 1e-9   # (very rough)
end
function C32(T, ne)
    # three-body recombination / collisional recombination estimate
    return ne * 1e-12
end

# ------------------------- Angle quadrature (4-point Gauss-Legendre)
mus = [ -0.8611363116, -0.3399810436, 0.3399810436, 0.8611363116 ]
wmu = [ 0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451 ]

# ------------------------- Discretization for Feautrier scheme
# Optical depth grid: we will build optical depth from opacity later; for Feautrier we use tau grid
Nz_tau = Nz
# For now create a geometrical grid and map to optical depth via cumulative κ later

# Helper: tridiagonal solver (Thomas algorithm)
function tridiag_solver(a, b, c, r)
    n = length(b)
    # a: sub-diagonal (n-1), b: diagonal (n), c: super-diagonal (n-1)
    cp = similar(c)
    dp = similar(r)
    cp[1] = c[1] / b[1]
    dp[1] = r[1] / b[1]
    for i in 2:n-1
        denom = b[i] - a[i-1]*cp[i-1]
        cp[i] = c[i] / denom
        dp[i] = (r[i] - a[i-1]*dp[i-1]) / denom
    end
    dp[n] = (r[n] - a[n-1]*dp[n-1]) / (b[n] - a[n-1]*cp[n-1])
    x = similar(r)
    x[n] = dp[n]
    for i in n-1:-1:1
        x[i] = dp[i] - cp[i]*x[i+1]
    end
    return x
end

# ------------------------- Formal solver (Feautrier) for given source function S(τ,ν)
# We solve d^2u/dτ^2 = u - S with boundary: at top (τ=0): du/dτ = u (no incident) ???
# We'll use approximate boundary conditions for zero incident intensity: du/dτ = u at top, du/dτ = -u at bottom (symmetric)

function feautrier(S, τgrid, μsign=1)
    # S: vector length Nz (source function at each depth), τgrid in same length
    N = length(τgrid)
    # set up second-order finite difference for u at grid points
    # coefficients for tridiagonal system: a_i u_{i-1} + b_i u_i + c_i u_{i+1} = r_i
    a = zeros(N-1); b = zeros(N); c = zeros(N-1); r = zeros(N)
    # build using central difference with non-uniform Δτ
    dτ = diff(τgrid)
    # interior points
    for i in 2:N-1
        dτim = dτ[i-1]; dτip = dτ[i]
        a[i-1] = 2/(dτim*(dτim+dτip))
        c[i]   = 2/(dτip*(dτim+dτip))
        b[i]   = -(a[i-1] + c[i]) - 1.0   # because equation: d^2u/dτ^2 - u = -S -> discrete gives ... - u = -S
        r[i]   = -S[i]
    end
    # top boundary i=1 (τ=0) -> approximate: du/dτ = u  => (u2-u1)/dτ1 = u1  => u2 - u1 = dτ1*u1
    # rearrange to get equation for u1
    dτ1 = dτ[1]
    b[1] = -1.0 - 1.0/dτ1
    c[1] = 1.0/dτ1
    r[1] = -S[1]
    # bottom boundary i=N: du/dτ = -u  => (uN - uN-1)/dτ_{N-1} = -uN  => uN - uN-1 = -dτ_{N-1} uN
    dτm = dτ[end]
    a[end] = 1.0/dτm
    b[end] = -1.0 - 1.0/dτm
    r[end] = -S[end]
    # Solve tridiagonal; but note a length conventions. We prepared a[1..N-1], c[1..N-1]
    return tridiag_solver(a, b, c, r)
end

# From u(τ) we can reconstruct intensities at surface using I(0, μ) ≈ u(0) + μ du/dτ|0, but we'll approximate emergent intensity
# Simpler: we will compute specific intensity by short-characteristics formal integrator using S and optical depths for each μ

function formal_solver_short_characteristics(S, κν, τgrid, μ)
    # Integrate along ray from bottom to top for μ>0 and top to bottom for μ<0
    N = length(τgrid)
    I = zeros(N)
    if μ > 0
        # boundary at bottom (no incident): I_bottom = 0
        I[N] = 0.0
        for i in N-1:-1:1
            dτ = (τgrid[i+1] - τgrid[i])/abs(μ)
            # formal solution: I_i = I_{i+1} e^{-dτ} + ∫_0^{dτ} S(t) e^{-(dτ - t)} dt
            # approximate S linear between points
            S_av = 0.5*(S[i+1] + S[i])
            I[i] = I[i+1]*exp(-dτ) + S_av*(1.0 - exp(-dτ))
        end
    else
        # μ < 0, integrate from top to bottom; top boundary I_top = 0
        I[1] = 0.0
        for i in 1:N-1
            dτ = (τgrid[i+1] - τgrid[i])/abs(μ)
            S_av = 0.5*(S[i+1] + S[i])
            I[i+1] = I[i]*exp(-dτ) + S_av*(1.0 - exp(-dτ))
        end
    end
    return I
end

# simple trapezoid integration function
function trapz(x, y)
    s = 0.0
    for i in 1:length(x)-1
        s += 0.5*(y[i] + y[i+1])*(x[i+1] - x[i])
    end
    return s
end


# ------------------------- Main solver loop -------------------------
# Initialize populations with Boltzmann for 3-levels (toy): n_i = g_i * exp(-E_i/kT) normalized
pop0 = [g1*exp(-E1/(kB*T)), g2*exp(-E2/(kB*T)), 1.0] # continuum approx
pop0 = pop0 ./ sum(pop0) * n_total

# Initialize populations at each depth (start LTE-ish)
n1 = fill(pop0[1], Nz)
n2 = fill(pop0[2], Nz)
n3 = fill(pop0[3], Nz)

# Continuum (small) opacity per atom (toy)
κ_cont_per_atom = 1e-20

# Frequency loop will call radiative transfer for each ν. But to speed up we will vectorize over ν in the outer loop
maxiter = 120
tol = 1e-4

τ = zeros(Nz)

println("Starting NLTE iteration...")
for iter = 1:maxiter
    # Build opacity and emissivity at each depth and frequency given current populations
    # We'll build optical depth scale at line center for mapping geometry->tau: κ(τ) = κ_cont + κ_line
    global τ
    global κ_total

    κ_line_center = zeros(Nz)
    for i in 1:Nz
        # absorption coefficient per atom for line at center: α0 = (h ν / 4π) * (n1 B12 - n2 B21) ???
        # We'll use simpler classical: χ = n_lower * σ0 - n_upper * σ0 * (g_lower/g_upper)
        σ0 = (π * eV) * 1e-17  # toy cross-section scale (nonsense number but works to set opacity)
        κ_line_center[i] = max(0.0, n1[i]*σ0 - n2[i]*σ0*(g1/g2))
    end
    # For a monotonic τ scale take τ = cumulative integral of κ_line_center + continuum
    κ_total = κ_line_center .+ κ_cont_per_atom .* (n1 .+ n2 .+ n3)
    # map to optical depth by integrating from top to bottom
#    τ = zeros(Nz)
    for i in 2:Nz
        # simple trapezoid over geometrical coordinate (depth spacing assumed uniform)
        Δz = abs(depth[i] - depth[i-1]) + 1e-12
       τ[i] = τ[i-1] + 0.5*(κ_total[i] + κ_total[i-1]) * Δz
    end

    # For each frequency compute J(ν) by solving RTE
    Jν = zeros(nλ, Nz)
   
     for (inu, νval) in enumerate(ν)  # use 'inu' for the frequency index
        # Build monochromatic opacity including profile
        φs = [phi_nu(νval, ν12, Δν_D) for _ in 1:Nz]  # same φ at each depth
        κν = κ_total .* (1.0 .+ 0.0) .+ (n1 .* B12 .- n2 .* B21) .* (h*νval/(4*pi)) .* φs
        κν = abs.(κν) .+ 1e-30
        emissivity = (h*νval/(4*pi)) .* n2 .* A21 .* φs  # renamed to avoid collision
        S = (emissivity ./ κν) .+ 1e-6

        # Solve intensities for each angle and accumulate J
        Iμ = zeros(length(mus), Nz)
        for (ia, μ) in enumerate(mus)
            Iμ[ia, :] = formal_solver_short_characteristics(S, κν, τ, μ)
        end

        # mean intensity J(ν) at each depth
        for iz in 1:Nz
            Jν[inu, iz] = 0.5 * sum(wmu .* Iμ[:, iz])
        end
    end


    # Compute radiative rates integrated over line profile
    # For bound-bound 1-2: R12 = B12 * J_bar, R21 = A21 + B21 * J_bar, where J_bar = ∫ φ(ν) Jν dν
    Jbar = zeros(Nz)
    for iz in 1:Nz
        # integrate over ν with simple trapezoid
        integrand = [phi_nu(νval, ν12, Δν_D) * Jν[k, iz] for (k, νval) in enumerate(ν)]
        Jbar[iz] = trapz(collect(ν), integrand)
    end

    R12 = B12 .* Jbar
    R21 = A21 .+ B21 .* Jbar

    # For bound-free (2-3) assume a weak radiative ionization rate R23 (toy constant)
    R23 = fill(1e-3, Nz)
    R32 = fill(1e-6, Nz)

    # Now solve statistical equilibrium for each depth: steady-state -> matrix eqn A * n = 0 with normalization
    n1_new = similar(n1); n2_new = similar(n2); n3_new = similar(n3)
    for iz in 1:Nz
        C12v = C12(T, ne); C21v = C21(T, ne);
        C23v = C23(T, ne); C32v = C32(T, ne);
        # Rate matrix P_ij: transitions i->j
        P12 = R12[iz] + C12v
        P21 = R21[iz] + C21v
        P23 = R23[iz] + C23v
        P32 = R32[iz] + C32v
        # steady state: for levels 1,2,3
        # -n1*(P12 + P13) + n2*P21 + n3*P31 = 0
        # -n2*(P21 + P23) + n1*P12 + n3*P32 = 0
        # normalization: n1 + n2 + n3 = n_total
        # For simplicity assume P13=P31=0 (no direct 1<->3 radiative), can be extended
        A = [-(P12)  P21   0.0;
              P12   -(P21+P23) P32;
              1.0    1.0    1.0]
        b = [0.0, 0.0, n_total]
        # solve linear system (3x3)
        nvec = A \ b
        n1_new[iz] = max(nvec[1], 1e-20)
        n2_new[iz] = max(nvec[2], 1e-20)
        n3_new[iz] = max(nvec[3], 1e-20)
    end

    # check convergence
    err = maximum(abs.((n2_new .- n2) ./ (n2 .+ 1e-30)))
    println("iter $iter, max rel change in n2 = $(round(err, sigdigits=4))")
    n1 .= n1_new; n2 .= n2_new; n3 .= n3_new
    if err < tol
        println("Converged after $iter iterations.")
        break
    end
    if iter == maxiter
        println("Reached maximum iterations ($maxiter) without full convergence.")
    end
end

# ------------------------- Compute emergent flux (integrate over μ and wavelength)
# For final populations compute emergent intensity at surface for each ν and angles

Fnu = zeros(nλ)
for (inu, νval) in enumerate(ν)  # inu = integer index
    φs = [phi_nu(νval, ν12, Δν_D) for _ in 1:Nz]
    κν = κ_total .* (1.0 .+ 0.0) .+ (n1 .* B12 .- n2 .* B21) .* (h*νval/(4*pi)) .* φs
    κν = abs.(κν) .+ 1e-30
    emissivity = (h*νval/(4*pi)) .* n2 .* A21 .* φs  # renamed
    S = (emissivity ./ κν)

    Iμ = zeros(length(mus), Nz)
    for (ia, μ) in enumerate(mus)
        Iμ[ia, :] = formal_solver_short_characteristics(S, κν, τ, μ)
    end

    I_at_surface = Iμ[:, 1]
    Fnu[inu] = 2*pi * sum(wmu .* max.(0.0, mus) .* I_at_surface)
end

print(Fnu)

# Convert frequency to wavelength and plot flux vs wavelength
lam_ang = lam .* 1e8   # Angstroms
plot(lam_ang, Fnu, xlabel="Wavelength (Å)", ylabel="Flux (arb)", title="Emergent spectral line (toy 3-level atom)")
savefig("emergent_line.png")
println("Saved emergent_line.png in current directory.")
println("Done.")

