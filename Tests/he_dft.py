#!/usr/bin/python
"He using Gaussians to test the DFT module"

from PyQuante.dft import *
from PyQuante.Molecule import Molecule

energy = -2.8266976
name = "He"


def main():
    en,orbe,orbs = dft(Molecule('He',atomlist=[(2,(0,0,0))]))
    return en

if __name__ == '__main__': main()
