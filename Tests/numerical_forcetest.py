#!/usr/bin/env python
from PyQuante.Molecule import Molecule
from PyQuante.Ints import *

bl=1.4
delta=0.000001

mol1 = Molecule('H2',[(1,(0.,0.,0.)),(1,(0.,0.,bl))])	
bfs = getbasis(mol1.atoms,'sto-3g')

Sone,hone = get1ints(bfs,mol1.atoms)
tone = getT(bfs)
vone = getV(bfs,mol1.atoms)	


mol2 = Molecule('H2',[(1,(0.,0.,0.)),(1,(0.,0.,bl+delta))])
bfs = getbasis(mol2.atoms,'sto-3g')

Stwo,htwo = get1ints(bfs,mol2.atoms)
ttwo = getT(bfs)
vtwo = getV(bfs,mol2.atoms)

print "\n***Numerical derivatives***\n"
print "H2 Bondlength: %f" %bl
print "Delta: %f" %delta

print "dHcore/dZa\n",(htwo-hone)/delta

print "dT/dZ\n",(ttwo-tone)/delta

print "dV/dZ\n",(vtwo-vone)/delta

print "dS/dZ\n",(Stwo-Sone)/delta
