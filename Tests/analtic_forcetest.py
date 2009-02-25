#!/usr/bin/env python
from PyQuante.hartree_fock import rhf
from PyQuante.Molecule import Molecule
from PyQuante.force import *


bl=1.4

mol = Molecule('H2',[(1,(0.,0.,0.)),(1,(0.0,0.0,bl))])
en,orb_en,coefs = rhf(mol,MaxIter=20,basis_data="sto-3g")

bfs = getbasis(mol,'sto-3g')

dHcore_dXa,dHcore_dYa,dHcore_dZa = der_Hcore_matrix(1,bfs,mol.atoms)

dS_dXa,dS_dYa,dS_dZa  = der_overlap_matrix(1,bfs)

print "\n*** Analytic gradients ***\n"

print "\ndHcore/dZa\n",dHcore_dZa

print "\ndS/dZa\n",dS_dZa
