#!/usr/bin/env python
"""\
Solvers.py - explores ways to use different eigensolvers in Python
"""

from PyQuante import HFSolver, Molecule, logging
from PyQuante.NumWrap import eigh,zeros,matrixmultiply,transpose,dot
from math import sqrt

class SubspaceSolver(HFSolver):
    def __init__(self,molecule,solver,**opts):
        HFSolver.__init__(self,molecule,**opts)
        self.solver = solver
        return
    
    def solve_fock(self):
        from PyQuante.NumWrap import matrixmultiply, eigh, transpose
        from PyQuante.LA2 import simx
        # Transform F to the subspace of the previous iteration's orbitals:
        F = simx(self.F,self.orbs)
        # Solve in the subspace:
        self.orbe,orbs = self.solver(F)
        # Update the orbitals
        self.orbs = matrixmultiply(self.orbs,orbs)
        self.update_density()
        return

class DmatSolver(HFSolver):
    def __init__(self,molecule,solver,**opts):
        HFSolver.__init__(self,molecule,**opts)
        self.solver = solver
        return

    def solve_fock(self):
        self.D = self.solver(self.F, self.S, self.nclosed)
        return

### General functions required for the davidson solver
def appendColumn(A,newVec):
    """\
    Append a column vector onto matrix A; this creates a new
    matrix and does the relevant copying.
    """
    n,m = A.shape
    Anew = zeros((n,m+1),'d')
    Anew[:,:m] = A
    Anew[:,m] = newVec
    return Anew

def orthog(q,qs,**kwopts):
    "Orthonormalize vector q to set of vectors qs"
    nbasis,nvec = qs.shape
    north = kwopts.get('north',nvec)
    #if north==-1: north = nvec
    #north = min(north,nvec)
    for i in range(north):
        olap = dot(q,qs[:,i])
        q -= olap*qs[:,i]
    norm = sqrt(dot(q,q))
    return norm

def davidson(A,nroots,**kwargs):
    etol = kwargs.get('etol',1e-6) # tolerance on the eigenvalues before declaring converged
    ntol = kwargs.get('ntol',1e-10) # tolerance on the norms for adding a vector
    n,m = A.shape
    ninit = max(nroots,2)
    B = zeros((n,ninit),'d')
    for i in range(ninit): B[i,i] = 1.
    #from numpy.random import rand
    #B = rand(n,1)

    nc = 0 # number of converged roots
    eigold = 1e10
    for iter in range(n):
        if nc >= nroots: break
        D = matrixmultiply(A,B)
        #olap = matrixmultiply(transpose(B),B)
        S = matrixmultiply(transpose(B),D)
        m = len(S)
        eval,evec = eigh(S)

        bnew = zeros(n,'d')
        for i in range(m):
            bnew += evec[i,nc]*(D[:,i] - eval[nc]*B[:,i])

        for i in range(n):
            denom = max(eval[nc]-A[i,i],1e-8) # Set a maximum amplification factor
            bnew[i] /= denom

        norm = orthog(bnew,B)
        bnew = bnew / norm

        if abs(eval[nc]-eigold) < etol:
            nc += 1
        eigold = eval[nc]
        if norm > ntol: B = appendColumn(B,bnew)

    E = eval[:nroots]
    nv = len(evec)
    V = matrixmultiply(B[:,:nv],evec)
    return E,V

# Easy way to remove arguments from Davidson, allowing one to make a davidson() call
# into a single-argument call, along the same lines as the new partial function application
# stuff in Py2.5
def init_davidson(nroots,**opts):
    def func(A):
        return davidson(A,nroots,**opts)
    return func

def test():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(message)s")
    h2 = Molecule('H2',atomlist=[(1,(0.35,0,0)),(1,(-0.35,0,0))],units='Angs')

    logging.info("\nRegular eigensolver")
    h2_normal = HFSolver(h2)
    h2_normal.iterate()

    logging.info("\nNormal eigensolver, solved in the subspace of existing orbitals")
    h2_sub = SubspaceSolver(h2,eigh)
    h2_sub.iterate()

    logging.info("\nDavidson eigensolver in the subspace of existing orbitals")
    dav = init_davidson(2) # Have to look for more than 1 root
    h2_dav = SubspaceSolver(h2,dav)
    h2_dav.iterate()

    logging.info("\nDensity Matrix Purification")
    from PyQuante.dmm import DMP
    h2solv = DmatSolver(h2,DMP)
    h2solv.iterate()
    
if __name__ == '__main__': test()
