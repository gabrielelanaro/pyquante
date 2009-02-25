#!/usr/bin/env python
"""\
 Test of the EN2 module.
"""

from PyQuante.EN2 import EN2
from PyQuante.Molecule import Molecule, toBohr3
from PyQuante.Element import symbol

def main():
    multn = [None,2,1,2,1,2,3,4,3,2,1]    # multiplicity of neutral
    multc = [None,None,2,1,2,1,2,3,4,3,2] # multiplicity of cation
    for atno in range(1,11):
        print symbol[atno]
        neutral = Molecule(symbol[atno],
                           [(atno,(0,0,0))],
                           multiplicity=multn[atno])
        en1 = EN2(neutral)
        if atno > 1:
            cation = Molecule(symbol[atno]+"+",
                              [(atno,(0,0,0))],
                              charge=1,
                              multiplicity=multc[atno])
            en2 = EN2(cation)
        else:
            en2 = 0
        ie = (en2-en1)*27.2114
        print en1,en2,ie

if __name__ == '__main__': main()
