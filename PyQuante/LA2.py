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

from PyQuante import logging
from math import sqrt
from NumWrap import matrixmultiply,transpose,diagonal,identity,zeros,eigh

# Note: to be really smart in a quantum chemistry program, we would
#  want to only symmetrically orthogonalize the S matrix once, since
#  the matrix doesn't change during the SCF procedure. Thus, we would
#  want to call X = SymOrth(S), and then geighD(H,X), rather
#  than calling geigh(H,S) every time, since the latter
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
    if trans.lower().startswith('t'):
        return matrixmultiply(B,matrixmultiply(A,transpose(B)))
    return matrixmultiply(transpose(B),matrixmultiply(A,B))

def outprod(A):
    "D = outprod(A) : Return the outer product A*A'"
    return matrixmultiply(A,transpose(A))

def geigh(H,A,**opts):
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
        return geigh(H,X,**opts)
    #val,vec = eigh(SimilarityTransformT(H,A))
    val,vec = eigh(simx(H,A))
    vec = matrixmultiply(A,vec)
    return val,vec

def SymOrth(X):
    """Symmetric orthogonalization of the real symmetric matrix X.
    This is given by Ut(1/sqrt(lambda))U, where lambda,U are the
    eigenvalues/vectors."""
    val,vec = eigh(X)
    n = vec.shape[0]
    shalf = identity(n,'d')
    for i in range(n):
        shalf[i,i] /= sqrt(val[i])
    X = simx(shalf,vec,'T')
    #X = SimilarityTransform(shalf,vec)
    return X

def CanOrth(X): 
    """Canonical orthogonalization of matrix X. This is given by
    U(1/sqrt(lambda)), where lambda,U are the eigenvalues/vectors."""
    n = vec.shape[0]
    val,vec = eigh(X)

    for i in range(n):
        vec[:,i] = vec[:,i] / sqrt(val[i])
    return vec

def trace2(H,D):
    "Return the trace(H*D), used in computing QM energies"
    return sum(sum(H*D)) # O(N^2) version 

# SimilarityTransform/T are deprecated in favor of simx
def SimilarityTransformT(H,X):
    "Return the similarity transformation XtHX of H"
    logging.warning("SimilarityTransformT is deprecated: use simx")
    return simx(H,X)
    #return matrixmultiply(transpose(X),matrixmultiply(H,X))

def SimilarityTransform(H,X): 
    "Return the transpose similarity transformation XHXt of H"
    logging.warning("SimilarityTransform is deprecated: use simx")
    return simx(H,X,'T')
    #return matrixmultiply(X,matrixmultiply(H,transpose(X)))

def mkdens(c,nstart,nstop):
    "Form a density matrix C*Ct given eigenvectors C[nstart:nstop,:]"
    d = c[:,nstart:nstop]
    Dmat = matrixmultiply(d,transpose(d))
    return Dmat

def mkdens2(c,nstart,nstop):
    "2*normal density matrix, since that's more common"
    return 2*mkdens(c,nstart,nstop)

def mkdens_spinavg(c,nclosed,nopen):
    """Form a spin averaged density matrix with *nclosed* closed
       shell orbitals and *nopen* open shell orbitals"""
    return mkdens(c,0,nclosed) + 0.5*mkdens(c,nclosed,nclosed+nopen)
    
