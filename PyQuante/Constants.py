"""\
 Constants.py: Constants for PyQuante

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

# Misc units
Clight=2.99792458e8     # speed of light in m/s
Kboltz=3.166830e-6      # Boltzmann constant
e2 = 14.399             # Coulomb's law coeff if R in \AA and resulting E in eV
planck=6.6260755e-34    # Planck's constant, in Js

# Distance units
bohr2ang = 0.529177249  # Conversion of length from bohr to angstrom
ang2bohr = 1/bohr2ang

# Energy units
hartree2kcal = 627.5095 # Hartree to kcal/mol conversion
kcal2hartree = 1/hartree2kcal

ev2kcal = 23.061        # Conversion of energy in eV to energy in kcal/mol
kcal2ev = 1/ev2kcal

hartree2joule = 4.3597482e-18   # Hatree to Joule conversion factor
joule2hartree = 1/hartree2joule

ev2hartree = hartree2kcal/ev2kcal
hartree2ev = 1/ev2hartree

# Mass units
amu2me = 1822.882       # Conversion from mass in amu to mass in au (m_e)
me2amu = 1/amu2me       # Conversion from mass in au (m_e) to mass in amu 


# Time units
tau2ps = 41341.447      # Conversion from time in au to time in ps
ps2tau = 1/tau2ps       # inverse

# Derived quantities
Rgas = Kboltz*hartree2kcal*1000.0 # gas constant R = 1.98722 cal/mole/K

