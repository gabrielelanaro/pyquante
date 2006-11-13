#!/usr/bin/env python
"Helium using Gaussians"

# Jaguar results for 6-31G** basis
# Energy   -2.855 160
# Eone     -3.882 068
# Etwo      1.026 907

# GAMESS-UK result
# Energy -2.85516043

energy = -2.855223
energy = -2.85516040352   # After DIIS
name = "He"

from PyQuante.hartree_fock import rhf
from PyQuante.Ints import getbasis,getints
from PyQuante.Molecule import Molecule

def main():
    he = Molecule('He',atomlist = [(2,(0,0,0))])
    en,orbe,orbs = rhf(he)
    return en

if __name__ == '__main__': main() 
