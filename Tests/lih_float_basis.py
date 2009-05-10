#!/usr/bin/env

from copy import copy
from math import sqrt
from itertools import izip

from PyQuante import Molecule
from PyQuante.Ints import getbasis,getints,get2JmK
from PyQuante.LA2 import simx,SymOrth,identity
from PyQuante.NumWrap import matrixmultiply,eigh,zeros

def maxdiff_numpy(A,B): return max((A-B).ravel())
def maxdiff(A,B): return max(ai-bi for ai,bi in izip(A,B))

def SymOrthCutoff(S,scut=1e-5):
    """Symmetric orthogonalization of the real symmetric matrix S.
    This is given by Ut(1/sqrt(lambda))U, where lambda,U are the
    eigenvalues/vectors.

    Only eigenvectors with eigenvalues greater that a cutoff are kept.
    This approximation is useful in cases where the basis set has
    linear dependencies.
    """
    val,vec = eigh(S)
    n = vec.shape[0]
    shalf = identity(n,'d')
    for i in range(n):
        if val[i] > scut:
            shalf[i,i] /= sqrt(val[i])
        else:
            shalf[i,i] = 0.
    X = simx(shalf,vec,'T')
    return X

def geigh(H,A,**opts):
    """\
    Generalized eigenproblem using a symmetric matrix H.
    """
    X = SymOrthCutoff(A)
    val,vec = eigh(simx(H,X))
    vec = matrixmultiply(X,vec)
    return val,vec

# We could do this all automatically, but let's make a special purpose
# routine to play with:
def simple_hf(atoms,S,h,Ints):
    from PyQuante.LA2 import mkdens
    from PyQuante.hartree_fock import get_energy
    orbe,orbs = geigh(h,S)
    nclosed,nopen = atoms.get_closedopen()
    enuke = atoms.get_enuke()
    nocc = nclosed
    eold = 0
    for i in range(15):
        D = mkdens(orbs,0,nocc)
        G = get2JmK(Ints,D)
        F = h+G
        orbe,orbs = geigh(F,S)
        energy = get_energy(h,F,D,enuke)
        print "%d %f" % (i,energy)
        if abs(energy-eold) < 1e-4: break
        eold = energy
    print "Final HF energy for system %s is %f" % (atoms.name,energy)
    return energy,orbe,orbs

# Data
Li_x1 = -0.200966
H_x1 = 1.399033
Li_x2 = -0.351691
H_x2 = 2.448309

# Construct a molecule:
LiH1 = Molecule('LiH1',[('Li',(Li_x1,0,0)),('H',(H_x1,0,0))],units='Angs')
bfs1 = getbasis(LiH1)
S1,h1,Ints1 = getints(bfs1,LiH1)
simple_hf(LiH1,S1,h1,Ints1)

# Make another molecule
LiH2 = Molecule('LiH2',[('Li',(Li_x2,0,0)),('H',(H_x2,0,0))],units='Angs')
bfs2 = getbasis(LiH2)
S2,h2,Ints2 = getints(bfs2,LiH2)
simple_hf(LiH2,S2,h2,Ints2)


# Make a superset of the two basis sets:
bfs_big = bfs1 + bfs2
# and make a basis set with it:
S1a,h1a,Ints1a = getints(bfs_big,LiH1)
simple_hf(LiH1,S1a,h1a,Ints1a)
# The energy is slightly lower, which shows the additional functions are
# doing something

# Alternatively, we could do things the other way around
S2a,h2a,Ints2a = getints(bfs_big,LiH2)
simple_hf(LiH2,S2a,h2a,Ints2a)
# The energy is again slightly lower, which shows the additional functions are
# doing something

# The ints are not all the same. The overlaps and the 2e ints should be,
# but the 1e ints are different, since the nuclear attraction ints are
# different
print maxdiff_numpy(S1a,S2a)
print maxdiff_numpy(h1a,h2a)
print maxdiff(Ints1a,Ints2a)


