# IAU Resolution B3 2015: https://arxiv.org/abs/1510.07674
# NIST constants: https://physics.nist.gov/cuu/Constants/

# Fundamental constants
cLight :: Float64 = 2.99792458e10      # cm/s                    | Speed of light in vacuum
eCharge :: Float64 = 4.8032068e-10     # esu (statcoulombs)      | Electron charge
NAvogadro :: Float64 = 6.02214076e23     # mol⁻¹                 | Avogadro's number
GGrav :: Float64 = 6.67259e-8  # cm³/(g·s²)                      | Gravitational constant

# Electromagnetic constants
mElectron :: Float64 = 9.1093897e-28       # g                   | Electron mass
mProton :: Float64 = 1.6726231e-24         # g                   | Proton mass
mNeutron :: Float64 = 1.6749286e-24        # g                   | Neutron mass  
atomMassUnit :: Float64 = 1.66053906660e-24  # g                 | Atomic mass unit (amu)

# Quantum constants
hPlanck :: Float64 = 6.6260755e-27     # erg·s                   | Planck constant
hcutPlanck :: Float64 = 1.05457266e-27  # erg·s                  | Reduced Planck constant (hPlanck/2π)
RydConst :: Float64 = 2.1798723611035e-11  #erg                  | Rydberg constant

# Thermodynamic constants
gasConst :: Float64 = 8.314462618e7        # erg/(mol·K)         | Ideal gas constant
kBoltzmann :: Float64 = 1.380649e-16    # erg/K                  | Boltzmann constant
sigmaSB :: Float64 = 5.670374419e-5  # erg/(cm²·s·K⁴)            | Stefan-Boltzmann constant

# Astronomical constants
## Solar values
MSun :: Float64 = 1.98847e33             # g                     | Solar mass
Rsun :: Float64 = 6.957e10          # cm                         | Solar radius
Lsun :: Float64 = 3.828e33      # erg/s                          | Solar luminosity
Ssun :: Float64 = 1.361e6         # erg/(cm²·s)                  | Earth's solar irradiance
TeffSun :: Float64 = 5772.0  # K                                 | Effective blackbody temperature
## Other astronomical constants
astrounit :: Float64 = 1.4959787066e13  # cm                     | Astronomical unit (AU)
parsec :: Float64 = 3.085677581e18             # cm              | Parsec
lightyear :: Float64 = 9.4607304725808e17     # cm               | Light year