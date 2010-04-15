"""\
Solvers.py - explores ways to use different eigensolvers in Python
"""

from PyQuante import Molecule
from PyQuante.NumWrap import eigh,zeros,matrixmultiply,transpose,dot,\
     identity,diagonal,array
from math import sqrt


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

def orthog(q,qs,**kwargs):
    "Orthonormalize vector q to set of vectors qs"
    nbasis,nvec = qs.shape
    north = kwargs.get('north',nvec)
    #if north==-1: north = nvec
    #north = min(north,nvec)
    for i in xrange(north):
        olap = dot(q,qs[:,i])
        q -= olap*qs[:,i]
    norm = sqrt(dot(q,q))
    return norm

def davidson(A,nroots,**kwargs):
    etol = kwargs.get('etol',1e-6) # tolerance on the eigenval convergence
    ntol = kwargs.get('ntol',1e-10) # tolerance on the vector norms for addn
    n,m = A.shape
    ninit = max(nroots,2)
    B = zeros((n,ninit),'d')
    for i in xrange(ninit): B[i,i] = 1.

    nc = 0 # number of converged roots
    eigold = 1e10
    for iter in xrange(n):
        if nc >= nroots: break
        D = matrixmultiply(A,B)
        S = matrixmultiply(transpose(B),D)
        m = len(S)
        eval,evec = eigh(S)

        bnew = zeros(n,'d')
        for i in xrange(m):
            bnew += evec[i,nc]*(D[:,i] - eval[nc]*B[:,i])

        for i in xrange(n):
            denom = max(eval[nc]-A[i,i],1e-8) # Maximum amplification factor
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


#  General routine and auxilliary functions for Jacobi
def jacobi(A,**kwargs):
    """\
    E,V = jacobi(A,**kwargs) - Solve the eigenvalues/vectors of matrix
                             A using Jacobi's method.
    Options:
    Name        Default  Definition
    tol         1e-10    The tolerance for an element to be declared zero
    max_sweeps  100      Maximum number of sweeps through the matrix
    """
    max_sweeps = kwargs.get('max_sweeps',100)
    tol = kwargs.get('tol',1e-10)
    n = len(A)
    V = identity(n,'d')
    b = diagonal(A)
    d = diagonal(A)
    z = zeros(n,'d')
    nrot = 0
    for irot in xrange(max_sweeps):
        sm = 0
        for ip in xrange(n-1):
            for iq in xrange(ip+1,n):
                sm += abs(A[ip,iq])
        if sm < tol:
            # Normal return
            return evsort(b,V)

        thresh = 0
        if  irot < 3: thresh = 0.2*sm/n/n

        for ip in xrange(n-1):
            for iq in xrange(ip+1,n):
                g = 100*abs(A[ip,iq])
                if irot > 3 and g < tol:
                    A[ip,iq] = 0
                elif abs(A[ip,iq]) > thresh:
                    h = d[iq]-d[ip]
                    if g < tol:
                        t = A[ip,iq]/h  # t = 1/(2\theta)
                    else:
                        theta = 0.5*h/A[ip,iq] #  eq 11.1.10
                        t = 1./(abs(theta)+sqrt(1.+theta*theta))
                        if theta < 0: t = -t
                    c = 1.0/sqrt(1+t*t)
                    s = t*c
                    tau = s/(1.0+c)
                    h = t*A[ip,iq]
                    z[ip] -= h
                    z[iq] += h
                    d[ip] -= h
                    d[iq] += h
                    A[ip,iq] = 0
                    for j in xrange(ip):
                        A[j,ip],A[j,iq] = rotate(A[j,ip],A[j,iq],s,tau)
                    for j in xrange(ip+1,iq):
                        A[ip,j],A[j,iq] = rotate(A[ip,j],A[j,iq],s,tau)
                    for j in xrange(iq+1,n):
                        A[ip,j],A[iq,j] = rotate(A[ip,j],A[iq,j],s,tau)
                    for j in xrange(n):
                        V[j,ip],V[j,iq] = rotate(V[j,ip],V[j,iq],s,tau)
                    nrot += 1
        for ip in xrange(n):
            b[ip] += z[ip]
            d[ip] = b[ip]
            z[ip] = 0
    else:
        print "Too many iterations"
    return None

def evsort(E,V):
    from PyQuante.NumWrap import argsort
    ind = argsort(E)
    E = E[ind]
    V = V[:,ind]
    return E,V

def rotate(g,h,s,tau): return g-s*(h+g*tau),h+s*(g-h*tau)

# Easy way to remove arguments from multi-argument calls allowing one
# to make them as a single-argument call, using closures
def init_davidson(nroots,**kwargs):
    def func(A):
        return davidson(A,nroots,**kwargs)
    return func

def init_jacobi(**kwargs):
    def func(A):
        return jacobi(A,**kwargs)
    return func

def test():
    import logging
    h2 = Molecule('H2',atomlist=[(1,(0.35,0,0)),(1,(-0.35,0,0))],units='Angs')

    #logging.info("\nRegular eigensolver")
    h2_normal = HFSolver(h2)
    h2_normal.iterate()

    logging.info("\nNormal eigensolver, in subspace of existing orbitals")
    h2_sub = SubspaceSolver(h2,eigh)
    h2_sub.iterate()

    logging.info("\nDavidson eigensolve in subspace of existing orbitals")
    dav = init_davidson(2) # Have to look for more than 1 root
    h2_dav = SubspaceSolver(h2,dav)
    h2_dav.iterate()

    logging.info("\nJacobi eigensolve in subspace of existing orbitals")
    jac = init_jacobi()
    h2_jac = SubspaceSolver(h2,jac)
    h2_jac.iterate()

    logging.info("\nDensity Matrix Purification")
    from PyQuante.DMP import TCP, init_dmat_solver
    solver = init_dmat_solver(TCP)
    h2solv = DmatSolver(h2,solver)
    h2solv.iterate()

    logging.info("\nCanonical Purification")
    from PyQuante.DMP import CP, init_dmat_solver
    solver = init_dmat_solver(CP)
    h2solv = DmatSolver(h2,solver)
    h2solv.iterate()
    
    logging.info("\nMcWeeny Purification")
    from PyQuante.DMP import McWeeny, init_dmat_solver
    solver = init_dmat_solver(McWeeny)
    h2solv = DmatSolver(h2,solver)
    h2solv.iterate()
    
    logging.info("\nTrace Resetting Purification")
    from PyQuante.DMP import TRP, init_dmat_solver
    solver = init_dmat_solver(TRP)
    h2solv = DmatSolver(h2,solver)
    h2solv.iterate()
    
if __name__ == '__main__': test()
