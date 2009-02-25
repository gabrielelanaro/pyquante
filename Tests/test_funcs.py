#!/usr/bin/env python
"""\
Quick check for whether the functionals work. Not meant to be exhaustive.
"""
from PyQuante.dft import dft
from PyQuante.Molecule import Molecule

results = {
    ('He','S0')   : -2.7146414,
    ('H2','S0')   : -1.0411670,
    ('Be','S0')   : -14.216448,
    ('Ne','S0')   : -127.395408,
    ('He','SVWN') : -2.8266697,
    ('H2','SVWN') : -1.135311,
    ('Be','SVWN') : -14.442029,
    ('Ne','SVWN') : -128.142338,
    # Note: The BLYP functions only have the first-order potential terms
    ('He','BLYP') : -2.894904,
    ('H2','BLYP') : -1.166493,
    ('Be','BLYP') : -14.653761,
    ('Ne','BLYP') : -128.878195,
    # My first-order PBE results
    ('He','PBE')  : -2.883036,
    ('H2','PBE')  : -1.163945,
    ('Be','PBE')  : -14.623466,
    ('Ne','PBE')  : -128.778327,
    #PBE results from Jaguar
    ('He','PBEJ')  : -2.884456,  
    ('H2','PBEJ')  : -1.164108,
    ('Be','PBEJ')  : -14.625458,
    ('Ne','PBEJ')  : -128.777705,
    }

def is_close(e,e0,tol=1e-5): return abs(e-e0) < tol

def main():
    he = Molecule('He',atomlist=[(2,(0,0,0))])
    be = Molecule('Be',atomlist=[(4,(0,0,0))])
    ne = Molecule('Ne',atomlist=[(10,(0,0,0))])
    h2 = Molecule('H2',atomlist = [(1,(0,0,0.70)),
                                   (1,(0,0,-0.70))])

    #mols = [he,h2,be]
    mols = [he,h2,be,ne]
    #mols = [ne]
    for mol in mols:
        #for functional in ['SVWN','BLYP','PBE']:
        #for functional in ['S0','PBE']:
        for functional in ['PBE']:
            en,orbe,orbs = dft(mol,functional = functional)
            if (mol.name,functional) in results:
                worked = is_close(en,results[mol.name,functional])
            else:
                worked = None
            print mol.name,functional,en,worked
    return
if __name__ == '__main__': main()
