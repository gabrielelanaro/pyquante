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

import logging
from math import sqrt
from NumWrap import matrixmultiply,identity,zeros,eigh,cholesky,inv

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
    return 0.5*(A+A.T)
def simx(A,B,trans='N'):
    """\
    C = simx(A,B,trans)
    Perform the similarity transformation C = B'*A*B (trans='N') or
    C = B*A*B' (trans='T').
    """
    if trans.lower().startswith('t'):
        return matrixmultiply(B,matrixmultiply(A,B.T))
    return matrixmultiply(B.T,matrixmultiply(A,B))

def outprod(A):
    "D = outprod(A) : Return the outer product A*A'"
    return matrixmultiply(A,A.T)

def geigh(H,A,**opts):
    """\
    Generalized eigenproblem using a symmetric matrix H.

    Options:
    have_xfrm  False   Need to form canonical transformation from S (default)
               True    A is the canonical transformation matrix
    orthog     'Sym'   Use Symmetric Orthogonalization (default)
               'Can'   Use Canonical Orthogonalization
               'Chol'  Use a Cholesky decomposition
               'Cut'   Use a symmetric orthogonalization with a cutoff
                
    """
    have_xfrm = opts.get('have_xfrm',False)
    orthog = opts.get('orthog','Chol') 
    if not have_xfrm:
        if orthog == 'Can':
            X = CanOrth(A)
        elif orthog == 'Chol':
            X = CholOrth(A)
        elif orthog == 'Cut':
            X = SymOrthCutoff(A)
        else:
            X = SymOrth(A)
        opts['have_xfrm'] = True
        return geigh(H,X,**opts)
    val,vec = eigh(simx(H,A))
    vec = matrixmultiply(A,vec)
    return val,vec

def SymOrth(S):
    """Symmetric orthogonalization of the real symmetric matrix S.
    This is given by Ut(1/sqrt(lambda))U, where lambda,U are the
    eigenvalues/vectors."""
    val,vec = eigh(S)
    n = vec.shape[0]
    shalf = identity(n,'d')
    for i in xrange(n):
        shalf[i,i] /= sqrt(val[i])
    X = simx(shalf,vec,'T')
    return X

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
    for i in xrange(n):
        if val[i] > scut:
            shalf[i,i] /= sqrt(val[i])
        else:
            shalf[i,i] = 0.
    X = simx(shalf,vec,'T')
    return X

def CanOrth(S): 
    """Canonical orthogonalization of matrix S. This is given by
    U(1/sqrt(lambda)), where lambda,U are the eigenvalues/vectors."""
    val,vec = eigh(S)
    n = vec.shape[0]
    for i in xrange(n):
        vec[:,i] = vec[:,i] / sqrt(val[i])
    return vec

def CholOrth(S):
    """Cholesky orthogonalization of matrix X. This is given by
    LL^T=S; X = transpose(inv(L))"""
    return inv(cholesky(S)).T

def trace2(H,D):
    "Return the trace(H*D), used in computing QM energies"
    return sum(sum(H*D)) # O(N^2) version 

def mkdens(c,nstart,nstop):
    "Form a density matrix C*Ct given eigenvectors C[nstart:nstop,:]"
    d = c[:,nstart:nstop]
    Dmat = matrixmultiply(d,d.T)
    return Dmat

def mkdens_spinavg(c,nclosed,nopen):
    """Form a spin averaged density matrix with *nclosed* closed
       shell orbitals and *nopen* open shell orbitals"""
    return mkdens(c,0,nclosed) + 0.5*mkdens(c,nclosed,nclosed+nopen)

#added by Hatem H Helal 18.07.2007
def pad_out(matrix):
    #this will make debugging matrix operations easier by getting rid of 
    #matrix elements which are ridiculously tiny
    
    print array2string(matrix,max_line_width=200,precision=7,suppress_small=True);    print "\n\n"

    return 0

def diagonal_mat(vector):
    #takes a vector with N components and returns an NXN matrix whose diagonal elements
    #are the components of the vector
    len = vector.shape[0]
        
    matrix = zeros((len,len),'d')
    
    i=0
    for i in xrange(len):
        matrix[i][i]=vector[i]

    return matrix    
    

