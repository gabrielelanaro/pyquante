"""\
 Minimizers.py: Geometry Minimizers

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

# Todo list:
# * conjugate gradient geometry minimizer
# * numerical hessian/frequencies

from math import sqrt
from IO import write_xyz

def shift_geo(atomlist,atom,dir,amount):
    atno,xyz = atomlist[atom]
    lxyz = list(xyz)
    lxyz[dir] += amount
    atomlist[atom] = atno,tuple(lxyz)
    return atomlist

def NumericForces(atomlist,EnergyFunction):
    "Return the forces on each atom in atomlist via finite differences"
    Forces = []
    nat = len(atomlist)
    for i in xrange(nat):
        plus_x_geo = shift_geo(copy(atomlist),i,0,dx)
        minus_x_geo = shift_geo(copy(atomlist),i,0,-dx)
        plus_y_geo = shift_geo(copy(atomlist),i,1,dx)
        minus_y_geo = shift_geo(copy(atomlist),i,1,-dx)
        plus_z_geo = shift_geo(copy(atomlist),i,2,dx)
        minus_z_geo = shift_geo(copy(atomlist),i,2,-dx)
        fx = (EnergyFunction(plus_x_geo)-EnergyFunction(minus_x_geo))/(2*dx)
        fy = (EnergyFunction(plus_y_geo)-EnergyFunction(minus_y_geo))/(2*dx)
        fz = (EnergyFunction(plus_z_geo)-EnergyFunction(minus_z_geo))/(2*dx)
        Forces.append((fx,fy,fz))
    return Forces

def Frms(F):
    "Compute the RMS value of a vector of forces"
    sqsum = 0
    for fx,fy,fz in F: sqsum += fx*fx+fy*fy+fz*fz
    return sqrt(sqsum/len(F))

def SteepestDescent(atomlist,EnergyForces):
    "Called with a pointer to an energy/force evaluator"
    file = open('out.xyz','w')
    step = 0.0001 # this seems awfully small
    Eold = None
    write_xyz(file,atomlist)
    for i in xrange(50):
        E,F = EnergyForces(atomlist)
        print i,E,Frms(F)
        for i in xrange(len(atomlist)):
            atno,(x,y,z) = atomlist[i]
            fx,fy,fz = F[i]
            x -= step*fx
            y -= step*fy
            z -= step*fz
            atomlist[i] = atno,(x,y,z)
        if Eold:
            dE = E-Eold
            if abs(dE) < 0.1: break
            if dE > 0:
                step *= 0.5
            else:
                step *= 1.2
        write_xyz(file,atomlist)
    file.close()
    return

