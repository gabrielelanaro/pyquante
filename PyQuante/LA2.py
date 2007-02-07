#!/usr/bin/env python
"""\
 LA2.py: Simple additions to numpy.linalg linear algebra library

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

# Todo
# - update SimilarityTransformation to simx


from math import sqrt
from NumWrap import test_numpy
from NumWrap import matrixmultiply,transpose,diagonal,identity,zeros
if test_numpy:
    from NumWrap import eigh
else:
    from NumWrap import Heigenvectors

# Note: to be really smart in a quantum chemistry program, we would
#  want to only symmetrically orthogonalize the S matrix once, since
#  the matrix doesn't change during the SCF procedure. Thus, we would
#  want to call X = SymOrth(S), and then GHeigenvectorsD(H,X), rather
#  than calling GHeigenvectors(H,S) every time, since the latter
#  recomputes the symmetric orthogonalization every SCF cycle.

def norm(vec):
    "val = norm(vec) : Return the 2-norm of a vector"
    return sqrt(dot(vec,vec))
def sym(A):
    "B = sym(A) : Symmetrize a matrix"
    return 0.5*(A+transpose(A))
def simx(A,B,trans='N'):
    """\
    C = simx(A,B,trans)
    Perform the similarity transformation C = B'*A*B (trans='N') or
    C = B*A*B' (trans='T').
    """
    if trans=='T': return matmul(B,matmul(A,transpose(B)))
    return matmul(transpose(B),matmul(A,B))

def outprod(A):
    "D = outprod(A) : Return the outer product A*A'"
    return matmul(A,transpose(A))

def GHeigenvectors(H,A,**opts):
    """\
    Generalized eigenproblem using a symmetric matrix H.

    Options:
    have_xfrm  False   Need to form canonical transformation from S (default)
               True    A is the canonical transformation matrix
    orthog     'Sym'   Use Symmetric Orthogonalization (default)
               'Can'   Use Canonical Orthogonalization
                
    """
    have_xfrm = opts.get('have_xfrm',False)
    orthog = opts.get('orthog','Sym')
    if not have_xfrm:
        if orthog == 'Can':
            X = CanOrth(A)
        else:
            X = SymOrth(A)
        opts['have_xfrm'] = True
        return GHeigenvectors(H,X,**opts)
    if test_numpy:
        val,vec = eigh(SimilarityTransformT(H,A))
        vec = matrixmultiply(A,vec)
    else:
        val,vec = Heigenvectors(SimilarityTransform(H,A))
        vec = matrixmultiply(vec,A)
    return val,vec

def SymOrth(X):
    """Symmetric orthogonalization of the real symmetric matrix X.
    This is given by Ut(1/sqrt(lambda))U, where lambda,U are the
    eigenvalues/vectors."""
    if test_numpy:
        val,vec = eigh(X)
    else:
        val,vec = Heigenvectors(X)
    n = vec.shape[0]
    shalf = identity(n,'d')
    for i in range(n):
        shalf[i,i] /= sqrt(val[i])
    if test_numpy:
        X = SimilarityTransform(shalf,vec)
    else:
        X = SimilarityTransformT(shalf,vec) 
    return X

def CanOrth(X): 
    """Canonical orthogonalization of matrix X. This is given by
    U(1/sqrt(lambda)), where lambda,U are the eigenvalues/vectors."""
    n = vec.shape[0]
    if test_numpy:
        val,vec = eigh(X)
    else:
        val,vec = Heigenvectors(X)

    if test_numpy:
        for i in range(n):
            vec[:,i] = vec[:,i] / sqrt(val[i])
    else:
        for i in range(n):
            for j in range(n):
                vec[i,j] = vec[i,j]/sqrt(val[i])
    return vec

def TraceProperty(H,D):
    "Return the trace(H*D), used in computing QM energies"
    return sum(diagonal(matrixmultiply(H,D)))
#  ??? Can we just do sum(ravel(H)*ravel(D)) here to make O(N2) ???

def SimilarityTransformT(H,X):
    "Return the similarity transformation XtHX of H"
    return matrixmultiply(transpose(X),matrixmultiply(H,X))

def SimilarityTransform(H,X): 
    "Return the transpose similarity transformation XHXt of H"
    return matrixmultiply(X,matrixmultiply(H,transpose(X)))

def mkdens(c,nstart,nstop):
    "Form a density matrix C*Ct given eigenvectors C[nstart:nstop,:]"
    if test_numpy:
        d = c[:,nstart:nstop]
        Dmat = matrixmultiply(d,transpose(d))
    else:
        d = c[nstart:nstop,:]
        Dmat = matrixmultiply(transpose(d),d)
    return Dmat

def mkdens2(c,nstart,nstop):
    "2*normal density matrix, since that's more common"
    return 2*mkdens(c,nstart,nstop)

def mkdens_spinavg(c,nclosed,nopen):
    """Form a spin averaged density matrix with *nclosed* closed
       shell orbitals and *nopen* open shell orbitals"""
    return mkdens(c,0,nclosed) + 0.5*mkdens(c,nclosed,nclosed+nopen)
    
