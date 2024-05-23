#!/usr/local/bin/python
# -*- coding: utf-8 -*-
"""
File : Constants.py 
For tests before implementing in Prop_Tau
Define cosmology & injection constants
"""

# Cosmological parameters
om_m = 0.3      # Omega matter
om_l = 0.7      # Omega lambda
H0 = 7.17e-11   # Hubble constant yr^{-1} (H0=70 Km/Mpc/s)

# Particle masses (in eV)
mN = 938.27e6   # Proton mass
mpi = 134.977e6 # Pion mass
mmu = 105.658e6 # Muon mass
me = 0.511e6    # Electron mass

# Pion production threshold (in MeV)
epsthr = (mpi + (mpi * mpi) / (2 * mN)) / 1e6

# Conversion factors
Mpc2cm = 3.086e24       # Mpc to cm
SpeedOfLight = 3e8      # Speed of light (m/s)
m2cm = 1e2              # meters to cm
SecInY = 3.1557e7       # year to seconds

# Mathematical constant
PI = 3.14159265358979323846

# Other constants
kREarth = 6371.0                # Earth radius (in km)
kNlat, kNlong = 91, 181         # Number of latitude and longitude points
kNTauolaFraction = 100000       # Tauola fraction constant
kRsun = 8.33                    # Solar radius (in kpc)
kpctokm = 3.086e16              # kpc to km conversion
kRh = 20.0                      # Hubble radius (in kpc)
kNAvogadro = 6.02e23            # Avogadro's number
alpha_min = 0 * PI / 180.0
alpha_max = 3 * PI / 180.0  # Minimum and maximum alpha values (in radians)
