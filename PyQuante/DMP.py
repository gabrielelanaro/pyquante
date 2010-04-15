"""\
 DMP.py - Implementation of density matrix purification methods, including 
  Niklasson, Tymczak and Challacombe's density matrix minimizer
  JCP 118, 8611 (2003)

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

# Status:
#  PM,TCP,TRS work for h2/* and h2o/*. Appears noisy for
#     h2o/631gss. Also, when lots of iterations are done in the
#     dm convergence, the program gives unreliable results. This
#     normally kicks in around 70 or 80 iterations.
#  MCW works, provided we have a guess for efermi

import logging
from math import sqrt
from PyQuante.Ints import getbasis,getints,get2JmK
from PyQuante.Molecule import Molecule
from PyQuante.LA2 import mkdens,SymOrth,simx
from PyQuante.hartree_fock import get_energy
from PyQuante.NumWrap import matrixmultiply,identity,trace,zeros,eigh,solve

logger = logging.getLogger("pyquante")

class AbstractDMP:
    "AbstractDMP - Functions common to all density matrix purifiers"
    method = "Abstract"
    def __init__(self,F,Ne,S=None,**opts):
        self.tol = opts.get('tol',1e-7)
        self.maxit = opts.get('maxit',100)
        self.do_orth = S is not None

        self.N = F.shape[0]
        self.I = identity(self.N,'d')
        self.Ne = Ne

        if self.do_orth:
            self.X = SymOrth(S)
            self.F = simx(F,self.X)

        self.emin, self.emax = gershgorin_minmax(self.F)
        self.initialize()
        self.print_init_info()
        return

    def iterate(self):
        for self.iter in xrange(self.maxit):
            if self.converged(): break
            self.update()
            self.print_iter_info()
        self.print_iter_end_info()
        if self.do_orth:
            self.D = simx(self.D,self.X,'T')
        return

    def reinitialize(self,F):
        "Used for restarting in a later SCF iteration"
        if self.do_orth:
            self.F = SymOrth(F,self.X)
        else:
            self.F = F
        self.emin,self.emax = gershgorin_minmax(self.F)
        self.initialize()
        self.print_init_info()
        return
            
    def print_iter_end_info(self):
        if self.iter == self.maxit-1:
            logger.warning("Too many iterations taken in %s: %d" %
                            (self.method,self.iter))
        else:
            logger.debug("%s converged in %d iterations" % 
                            (self.method,self.iter))
        return
    def print_iter_info(self): return
    def print_init_info(self): return

    # Functions that must be overloaded
    def initialize(self): print "AbstractDMP.initialize()"
    def update(self): print "AbstractDMP.update()"
    def converged(self): print "AbstractDMP.converged()"

class NOTCP:
    "Nonorthogonal version of Niklasson Trace Correcting Purification"
    method = "NOTCP"
    def __init__(self,F,Ne,S,**opts):
        self.tol = opts.get('tol',1e-7)
        self.maxit = opts.get('maxit',50)
        self.S = S
        self.N = F.shape[0]
        self.I = identity(self.N,'d')
        self.Ne = Ne
        self.F = F
        self.emin, self.emax = lanczos_minmax(self.F,self.S)
        self.initialize()
        self.print_init_info()
        return

    def iterate(self):
        for self.iter in xrange(self.maxit):
            if self.converged(): break
            self.update()
            self.print_iter_info()
        self.print_iter_end_info()
        return

    def reinitialize(self,F):
        "Used for restarting in a later SCF iteration"
        self.F = F
        self.emin,self.emax = lanczos_minmax(self.F,self.S)
        self.initialize()
        self.print_init_info()
        return
            
    def print_iter_end_info(self):
        if self.iter == self.maxit-1:
            logging.warning("Too many iterations taken in %s: %d" %
                            (self.method,self.iter))
        else:
            logging.debug("%s converged in %d iterations" % 
                            (self.method,self.iter))
        return

    def print_iter_info(self):
        #print self.iter,self.Ne,self.Ne_curr
        return
    
    def print_init_info(self): return

    def initialize(self):
        from PyQuante.NumWrap import inv
        self.D = inv(self.F-(self.emin-1)*self.S)
        return
    
    def update(self):
        D2 = matrixmultiply(self.DS,self.D)
        if self.Ne_curr < self.Ne:
            self.D = 2*self.D-D2
        else:
            self.D = D2
        return

    def converged(self):
        self.DS = matrixmultiply(self.D,self.S)
        self.Ne_curr = trace(self.DS)
        return abs(self.Ne_curr - self.Ne) < self.tol

class TCP(AbstractDMP):
    "Niklasson Trace Correcting Purification"
    method = "TCP"
    def initialize(self):
        self.D = (self.emax*self.I-self.F)/(self.emax-self.emin)
        return
    
    def update(self):
        Ne_curr = trace(self.D)
        D2 = matrixmultiply(self.D,self.D)
        Ne2 = trace(D2)
        self.Ne_curr = Ne_curr
        # Anders claims this works better; I didn't see a difference
        #if abs(2*Ne_curr-Ne2-self.Ne) < abs(Ne2-self.Ne):
        if Ne_curr < self.Ne:
            self.D = 2*self.D-D2
        else:
            self.D = D2
        return

    def converged(self):
        return abs(trace(self.D) - self.Ne) < self.tol

    def print_iter_info(self):
        #print self.iter,self.Ne,self.Ne_curr
        return
    
class TRP(TCP):
    "Niklasson/Tymczak/Challacombe Trace Resetting purification"
    method = "TRP"
    def update(self):
        D2 = matrixmultiply(self.D,self.D)
        Df = matrixmultiply(D2,4*self.D-3*D2)
        trf = trace(Df)
        Dp = self.I-self.D
        Dp2 = matrixmultiply(Dp,Dp)
        Dg = matrixmultiply(D2,Dp2)
        trg = trace(Dg)
        gamma = (self.Ne-trf)/trg
        if gamma > 2:
            self.D = 2*self.D-D2
        elif gamma < 0:
            self.D = D2
        else:
            self.D = Df-gamma*Dg
        return

class CP(AbstractDMP):
    "Palser/Manolopolous Canonical Purification"
    method = "CP"
    def initialize(self):
        efermi = trace(self.F)/self.N
        beta = self.Ne/float(self.N)
        alpha = min(self.Ne/(self.emax-efermi),
                    (self.N-self.Ne)/(efermi-self.emin))/float(self.N)
        self.D = alpha*(efermi*self.I-self.F) + beta*self.I
        self.Dsumold = 0
        return

    def update(self):
        D2 = matrixmultiply(self.D,self.D)
        D3 = matrixmultiply(self.D,D2)
        cn = trace(D2-D3)/trace(self.D-D2)
        if cn < 0.5:
            self.D = ((1.0-2.0*cn)*self.D+(1.0+cn)*D2-D3)/(1.0-cn)
        else:
            self.D = ((1+cn)*D2-D3)/cn
        return

    def converged(self):
        Dsum = sum(sum(self.D))
        val = abs(Dsum-self.Dsumold) < self.tol
        self.Dsumold = Dsum
        return val

class McWeeny(AbstractDMP):
    method = "MCW"
    def initialize(self):
        "Set efermi and create the initial D matrix"
        beta = 0.5
        elow = self.emin
        ehigh = self.emax+20
        de = ehigh-elow
        alpha = beta/de
        #nelow = trace(alpha*(elow*I-self.F) + beta*I)
        #nehigh = trace(alpha*(ehigh*I-self.F) + beta*I)
        nelow = self.get_nel(elow,alpha,beta)
        nehigh = self.get_nel(ehigh,alpha,beta)

        for i in xrange(100):
            efermi = 0.5*(elow+ehigh)
            #nefermi = trace(alpha*(efermi*I-F)+ beta*I)
            nefermi = self.get_nel(efermi,alpha,beta)
            if abs(self.Ne-nefermi) < self.tol: break
            if nefermi < self.Ne:
                elow = efermi
                nelow = nefermi
            else:
                ehigh = efermi
                nehigh = nefermi
        alpha = min(beta/(self.emax-efermi),(1-beta)/(efermi-self.emin))
        self.D = alpha*(efermi*self.I-self.F)+beta*self.I
        return

    def get_nel(self,efermi,alpha,beta):
        return trace(alpha*(efermi*self.I-self.F)*beta*self.I)

    def update(self):
        D2 = matrixmultiply(self.D,self.D)
        self.D = 3*D2-2*matrixmultiply(self.D,D2)
        return
    def converged(self): return abs(trace(self.D) - self.Ne) < self.tol


def init_dmat_solver(Method,**opts):
    "Wrapper around Dmat classes to make them work like simple solvers"
    def solver(F,S,Ne):
        solve = Method(F,Ne,S)
        solve.iterate()
        return solve.D
    return solver
    
def gershgorin_minmax(A):
    n,m = A.shape
    mins = []
    maxs = []
    for i in xrange(n):
        offsum = sum(abs(A[i,:]))-abs(A[i,i])
        mins.append(A[i,i]-offsum)
        maxs.append(A[i,i]+offsum)
    return min(mins),max(maxs)

def tridiagmat(alpha,beta):
    N = len(alpha)
    A = zeros((N,N),'d')
    for i in xrange(N):
        A[i,i] = alpha[i]
        if i < N-1:
            A[i,i+1] = A[i+1,i] = beta[i]
    return A

def lanczos_minmax(F,S=None,**kwargs):
    "Estimate the min/max evals of F using a few iters of Lanczos"
    doS = S is not None
    niter = kwargs.get('niter',8)
    N = F.shape[0]
    niter = min(N,niter)
    x = zeros(N,'d')
    x[0] = 1
    q = x
    avals = []
    bvals = []
    if doS:
        r = matrixmultiply(S,q)
    else:
        r = q
    beta = sqrt(matrixmultiply(q,r))
    wold = zeros(N,'d')
    for i in xrange(niter):
        w = r/beta
        v = q/beta
        r = matrixmultiply(F,v)
        r = r - wold*beta
        alpha = matrixmultiply(v,r)
        avals.append(alpha)
        r = r-w*alpha
        if doS:
            q = solve(S,r)
        else:
            q = r
        beta = sqrt(matrixmultiply(q,r))
        bvals.append(beta)
        wold = w
    E,V = eigh(tridiagmat(avals,bvals))
    return min(E),max(E)

def test():
    from PyQuante.PyQuante2 import SCF,DmatSolver
    print "Target energy: ",-1.130501
    h2 = Molecule('H2',atomlist=[(1,(0.35,0,0)),(1,(-0.35,0,0))],
                  units='Angs')
    h2_hf = SCF(h2,method='HF',SolverConstructor=DmatSolver)
    h2_hf.iterate()
    print "Energy:        ",h2_hf.energy
    
if __name__ == '__main__': test()

    
