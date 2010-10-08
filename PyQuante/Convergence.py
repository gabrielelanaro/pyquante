"""\
 Code for convergence acceleration in PyQuante.
 
 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from PyQuante.NumWrap import dot,ravel,matrixmultiply,zeros
from PyQuante.NumWrap import solve
from PyQuante.LA2 import SymOrth
from math import sqrt

VERBOSE=0

# This is simple density matrix averaging.
class SimpleAverager:
    def __init__(self,alpha=0.5):
        self.alpha = alpha
        self.Dold = None
        return

    def getD(self,D):
        if self.Dold is not None:
            Dreturn = self.alpha*D+(1-self.alpha)*self.Dold
        else:
            Dreturn = D
        self.Dold = D
        return Dreturn

# This is the averager that D.D. Johnson attributes to Anderson.
# PRB 38, 12807 (1988).
class AndersonAverager:
    def __init__(self,alpha=0.5):
        self.alpha = alpha
        self.Dold_in = None 
        self.Dold_out = None # D_{i-2}
        return

    def getD(self,Dold_in,Dold_out):
        if self.Dold_in:
            err = Dold_out-Dold_in
            err12 = err-self.err_old
            beta = dot1d(err,err12)/dot1d(err12,err12)
            print "beta = ",beta
            Dbar_in = (1-beta)*Dold_in + beta*self.Dold_in
            Dbar_out = (1-beta)*Dold_out + beta*self.Dold_out
            D = (1-self.alpha)*Dbar_in + self.alpha*Dbar_out
        else:
            D = (1-self.alpha)*Dold_in + self.alpha*Dold_out
            
        self.Dold_in = Dold_in
        self.Dold_out = Dold_out
        self.err_old = Dold_out-Dold_in
        return D

# Here are my notes on how the full Broyden scheme works.
# it input dens    wfn        output dens
# 1. \rho^1_in -> \psi^1 -> \rho^1_out
# 2. \rho^2_in = 0.3*\rho^1_in + 0.7*\rho^1_out -> \psi^2 -> \rho^2_out
# 3. \rho^3_in = 0.3*\rho^2_in + 0.7*\rho^2_out + gamma*\sum_{i=1}^1 w_i*\rho^i_in
#    \gamma = normalization, w_i = 2*\sqrt{frac{0.01}{\sum Err_i^2}}
# 4. \rho^4_in = 0.3*\rho^3_in + 0.7*\rho^3_out + gamma*\sum_{i=1}^2 ...
# etc.
# at each iteration you need to save \rho^i_in, and Err^i = \rho^i_out-\rho^i_in
#  (or, equivalently \rho^i_in and \rho^i_out)

# Pulay's DIIS
class DIIS:
    def __init__(self,S):
        self.Fs = []
        self.Errs = []
        self.Fold = None
        self.S = S
        #self.X = SymOrth(S)
        # Begin DIIS from iteration 0 in all cases
        self.started = True
        #self.started = 0
        self.errcutoff = 0.1
        return

    def error(self): return self.maxerr

    def getF(self,F,D):
        n,m = F.shape
        err = matrixmultiply(F,matrixmultiply(D,self.S)) -\
              matrixmultiply(self.S,matrixmultiply(D,F))
        err = ravel(err)
        maxerr = max(abs(err))
        self.maxerr = maxerr

        if maxerr < self.errcutoff and not self.started:
            if VERBOSE: print "Starting DIIS: Max Err = ",maxerr
            self.started = 1

        if not self.started:
            # Do simple averaging until DIIS starts
            if self.Fold != None:
                Freturn = 0.5*F + 0.5*self.Fold
                self.Fold = F
            else:
                self.Fold = F
                Freturn = F
            return Freturn

        self.Fs.append(F)
        self.Errs.append(err)
        nit = len(self.Errs)
        a = zeros((nit+1,nit+1),'d')
        b = zeros(nit+1,'d')
        for i in xrange(nit):
            for j in xrange(nit):
                a[i,j] = dot(self.Errs[i],self.Errs[j])
        for i in xrange(nit):
            a[nit,i] = a[i,nit] = -1.0
            b[i] = 0
        #mtx2file(a,'A%d.dat' % nit)
        a[nit,nit] = 0
        b[nit] = -1.0

        # The try loop makes this a bit more stable.
        #  Thanks to John Kendrick!
        try:
            c = solve(a,b)
        except:
            self.Fold = F
            return F
        
        F = zeros((n,m),'d')
        for i in xrange(nit):
            F += c[i]*self.Fs[i]
        return F

class DIIS2:
    # Two-point version of DIIS to save memory
    # This seemed like a good idea at the time, but it doesn't work
    #  particularly well.
    def __init__(self,S):
        self.Fold = None
        self.errold = None
        self.S = S
        self.started = 0
        self.nit = 0
        self.errcutoff = 0.5
        return

    def getF(self,F,D):
        n,m = F.shape
        err = matrixmultiply(F,matrixmultiply(D,self.S)) -\
              matrixmultiply(self.S,matrixmultiply(D,F))
        err = ravel(err)
        maxerr = max(abs(err))

        if maxerr < self.errcutoff and not self.started:
            if VERBOSE: print "Starting DIIS: Max Err = ",maxerr
            self.started = 1

        if not self.started:
            # Do simple averaging until DIIS starts
            if self.Fold:
                Freturn = 0.5*F + 0.5*self.Fold
            else:
                Freturn = F
            self.Fold = F
            return Freturn
        elif not self.errold:
            Freturn = 0.5*F + 0.5*self.Fold
            self.errold = err
            return Freturn

        a = zeros((3,3),'d')
        b = zeros(3,'d')
        a[0,0] = dot(self.errold,self.errold)
        a[1,0] = dot(self.errold,err)
        a[0,1] = a[1,0]
        a[1,1] = dot(err,err)
        a[:,2] = -1
        a[2,:] = -1
        a[2,2] = 0
        b[2] = -1
        c = solve(a,b)

        # Handle a few special cases:
        alpha = c[1]
        print alpha,c
        #if alpha < 0: alpha = 0
        #if alpha > 1: alpha = 1

        F = (1-alpha)*self.Fold + alpha*F
        self.errold = err
        self.Fold = F
        return F

def dot1d(a,b):
    return dot(ravel(a),ravel(b))

def test():
    from Ints import getbasis,getints,get2JmK
    from hartree_fock import get_nel, get_enuke,get_energy
    from LA2 import geigh,mkdens
    from IO import mtx2file
    from Molecule import Molecule

    ConvCriteria = 0.00001
    MaxIt = 30
    h2o = Molecule('h2o',[(8,(0,0,0)),(1,(1.,0,0)),(1,(0,1.,0))],
                   units='Angstrom')
    bfs = getbasis(h2o)
    S,h,Ints = getints(bfs,h2o)
    orbe,orbs = geigh(h,S)
    nel = get_nel(h2o)
    nocc = int(nel/2)
    enuke = get_enuke(h2o)
    eold = 0.
    avg = DIIS2(S)
    for i in xrange(30):
        D = mkdens(orbs,0,nocc)
        mtx2file(D)
        G = get2JmK(Ints,D)
        F = h+G
        F = avg.getF(F,D) # do the DIIS extrapolation
        orbe,orbs = geigh(F,S)
        energy = get_energy(h,F,D,enuke)
        print i+1,energy
        if abs(energy-eold) < ConvCriteria: break
        eold = energy

    return

if __name__ == '__main__': test()
