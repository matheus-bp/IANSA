# IAU Resolution B3 2015: https://arxiv.org/abs/1510.07674
# NIST constants: https://physics.nist.gov/cuu/Constants/

# Fundamental constants
cLight :: Real = 2.99792458e10      # cm/s                    | Speed of light in vacuum
eCharge :: Real = 4.8032068e-10     # esu (statcoulombs)      | Electron charge
NAvogadro :: Real = 6.02214076e23     # mol⁻¹                 | Avogadro's number
GGrav :: Real = 6.67259e-8  # cm³/(g·s²)                      | Gravitational constant

# Electromagnetic constants
mElectron :: Real = 9.1093897e-28       # g                   | Electron mass
mProton :: Real = 1.6726231e-24         # g                   | Proton mass
mNeutron :: Real = 1.6749286e-24        # g                   | Neutron mass  
atomMassUnit :: Real = 1.66053906660e-24  # g                 | Atomic mass unit (amu)

# Quantum constants
hPlanck :: Real = 6.6260755e-27     # erg·s                   | Planck constant
hcutPlanck :: Real = 1.05457266e-27  # erg·s                  | Reduced Planck constant (hPlanck/2π)
RydConst :: Real = 2.1798723611035e-11  #erg                  | Rydberg constant

# Thermodynamic constants
gasConst :: Real = 8.314462618e7        # erg/(mol·K)         | Ideal gas constant
kBoltzmann :: Real = 1.380649e-16    # erg/K                  | Boltzmann constant
sigmaSB :: Real = 5.670374419e-5  # erg/(cm²·s·K⁴)            | Stefan-Boltzmann constant

# Astronomical constants
## Solar values
MSun :: Real = 1.98847e33             # g                     | Solar mass
Rsun :: Real = 6.957e10          # cm                         | Solar radius
Lsun :: Real = 3.828e33      # erg/s                          | Solar luminosity
Ssun :: Real = 1.361e6         # erg/(cm²·s)                  | Earth's solar irradiance
TeffSun :: Real = 5772.0  # K                                 | Effective blackbody temperature
## Other astronomical constants
astrounit :: Real = 1.4959787066e13  # cm                     | Astronomical unit (AU)
parsec :: Real = 3.085677581e18             # cm              | Parsec
lightyear :: Real = 9.4607304725808e17     # cm               | Light year