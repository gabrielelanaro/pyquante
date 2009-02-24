"""\
 NEB.py - Nudged Elastic Band solvers in Python

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from Util import frange

def fhf(geo):
    return Molecule('fhf',
                    [(9,geo[0,0],geo[0,1],geo[0,2]),
                     (1,geo[1,0],geo[1,1],geo[1,2]),
                     (9,geo[2,0],geo[2,1],geo[2,2])],
                    units='Angstrom')

def eb(reactant,product,nimages):
    nbeads = nimages+2
    spacing = 1/float(nimages+1)
    coefs = frange(spacing,1,spacing)

    # Set the initial geometry for each bead
    geos = [None]*nbeads
    geos[0] = fhf(reactant)
    geos[nbeads] = fhf(product)
    for i in range(nimages):
        coef = coefs[i]
        geos[i+1] = fhf((1-coef)*reactant+coef*product)

    # Set the initial E,F for each bead
    Es = [0]*nbeads
    Fs = [0]*nbeads
    for i in range(nbeads):
        Es[i],Fs[i] = get_energy_forces(geos[i],-1)
    for iter in range(20):
        # compute forces on each bead
        for i in range(nimages):
            E[i+1],F[i+1] = get_energy_foces(geos[i+1],-1)

        # add bead-bead forces
        for i in range(nimages):
            geo = geos[i+1]
            geo_plus = geos[i+2]
            geo_minus = geos[i]
            for iatom in range(len(geo)):
                dx2p,dy2p,dz2p = geo[iatom].dist2(geo_plus[iatom])
                dx2m,dy2m,dz2m = geo[iatom].dist2(geo_minus[iatom])
                F[i+1][iatom,0] += 0.5*kbead*(dx2p+dx2m)
                F[i+1][iatom,1] += 0.5*kbead*(dy2p+dy2m)
                F[i+1][iatom,2] += 0.5*kbead*(dz2p+dz2m)
        
                

        
        
    return

def main():
    reactant = array([[-2.0, 0.0, 0.0],
                      [-1.0, 0.0, 0.0],
                      [ 2.0, 0.0, 0.0]])
    product = array([[-2.0, 0.0, 0.0],
                      [ 1.0, 0.0, 0.0],
                      [ 2.0, 0.0, 0.0]])
    eb(reactant,product,3,)

if __name__ == '__main__': main()
