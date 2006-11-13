#!/usr/bin/env python
"""\
 LA2.py: Simple additions to the Python LinearAlgebra library

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from math import sqrt
from NumWrap import matrixmultiply,transpose,diagonal,identity,zeros
from NumWrap import Heigenvectors


# Note: to be really smart in a quantum chemistry program, we would
#  want to only symmetrically orthogonalize the S matrix once, since
#  the matrix doesn't change during the SCF procedure. Thus, we would
#  want to call X = SymOrth(S), and then GHeigenvectorsD(H,X), rather
#  than calling GHeigenvectors(H,S) every time, since the latter
#  recomputes the symmetric orthogonalization every SCF cycle.

def Orthogonalize(H,S):
    X = SymOrth(S)
    return SimilarityTransform(H,X)

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
    val,vec = Heigenvectors(SimilarityTransform(H,A))
    return val,matrixmultiply(vec,A)

def SymOrth(X):
    """Symmetric orthogonalization of the real symmetric matrix X.
    This is given by Ut(1/sqrt(lambda))U, where lambda,U are the
    eigenvalues/vectors."""
    val,vec = Heigenvectors(X)
    n = vec.shape[0]
    shalf = identity(n,'d')
    for i in range(n):
        shalf[i,i] /= sqrt(val[i])
    return SimilarityTransformT(shalf,vec)

def CanOrth(X): 
    """Canonical orthogonalization of matrix X. This is given by
    U(1/sqrt(lambda)), where lambda,U are the eigenvalues/vectors."""
    val,vec = Heigenvectors(X)
    n = vec.shape[0]
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
    d = c[nstart:nstop,:]
    return matrixmultiply(transpose(d),d)

def mkdens2(c,nstart,nstop):
    "2*normal density matrix, since that's more common"
    return 2*mkdens(c,nstart,nstop)

def mkdens_spinavg(c,nclosed,nopen):
    """Form a spin averaged density matrix with *nclosed* closed
       shell orbitals and *nopen* open shell orbitals"""
    return mkdens(c,0,nclosed) + 0.5*mkdens(c,nclosed,nclosed+nopen)
    
