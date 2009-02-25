#!/usr/bin/env python
"""
Reproduce the data on the excellent Density Functional Repository web
site http://www.cse.clrc.ac.uk/qcg/dft/.

"""
from PyQuante.dft import dft
from PyQuante.Molecule import Molecule
from PyQuante.basis_dzvp import basis_data

def main():
    he = Molecule('He',atomlist=[(2,(0,0,0))])
    ne = Molecule('Ne',atomlist=[(10,(0,0,0))])
    he3 = Molecule('He3',atomlist=[(2,(0,0,0))],multiplicity=3)

    # S0 functional
    en,orbe,orbs = dft(he,functional='S0',basis=basis_data)
    print "He S0 sb %10.5f is %10.5f" % (-2.7229973821,en)
    en,orbe,orbs = dft(ne,functional='S0',basis=basis_data)
    print "Ne S0 sb %10.5f is %10.5f" % (-127.4597878033,en)
    en,orbe,orbs = dft(he3,functional='S0',basis=basis_data)
    print "He3 S0 sb %10.5f is %10.5f" % (-1.7819689849,en)

    # SVWN functional
    en,orbe,orbs = dft(he,functional='SVWN',basis=basis_data)
    print "He SVWN sb %10.5f is %10.5f" % (-2.834247,en)
    en,orbe,orbs = dft(ne,functional='SVWN',basis=basis_data)
    print "Ne SVWN sb %10.5f is %10.5f" % (-127.203239,en)
    en,orbe,orbs = dft(he3,functional='SVWN',basis=basis_data)
    print "He3 SVWN sb %10.5f is %10.5f" % (-1.833287,en)

    # XPBE functional
    return

if __name__ == '__main__': main()
