#!/usr/bin/env python
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
from PyQuante.Ints import getbasis,getints,get2JmK
from PyQuante.Molecule import Molecule
from PyQuante.LA2 import mkdens,SymOrth,simx
from PyQuante.hartree_fock import get_energy
from PyQuante.NumWrap import diagonal,matrixmultiply,identity,trace


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
        for self.iter in range(self.maxit):
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
            logging.warning("Too many iterations taken in %s: %d" %
                            (self.method,self.iter))
        else:
            logging.debug("%s converged in %d iterations" % 
                            (self.method,self.iter))
        return
    def print_iter_info(self): return
    def print_init_info(self): return

    # Functions that must be overloaded
    def initialize(self): print "AbstractDMP.initialize()"
    def update(self): print "AbstractDMP.update()"
    def converged(self): print "AbstractDMP.converged()"

class TCP(AbstractDMP):
    "Niklasson Trace Correcting Purification"
    method = "TCP"
    def initialize(self):
        self.D = (self.emax*self.I-self.F)/(self.emax-self.emin)
        return
    
    def update(self):
        D2 = matrixmultiply(self.D,self.D)
        Ne_curr = trace(self.D)
        if Ne_curr < self.Ne:
            self.D = 2.0*self.D-D2
        else:
            self.D = D2
        return

    def converged(self): return abs(trace(self.D) - self.Ne) < self.tol

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

        for i in range(100):
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
    for i in range(n):
        offsum = sum(abs(A[i,:]))-abs(A[i,i])
        mins.append(A[i,i]-offsum)
        maxs.append(A[i,i]+offsum)
    return min(mins),max(maxs)


### The following two routines are deprecated, and may be removed without warning:

def Dinit_mcw(F,Ne,tol=1e-7,maxit=100):
    #Solve for efermi and D0 using bisection:
    logging.warning("2/20/07: Dinit_mcw is now deprecated. Use DMP.McWeeny instead")
    beta = 0.5
    emin,emax = gershgorin_minmax(F)
    I = identity(F.shape[0],'d')
    elow = emin
    ehigh = emax+20
    de = emax-elow
    alpha = beta/de
    nelow = trace(alpha*(elow*I-F) + beta*I)
    nehigh = trace(alpha*(ehigh*I-F) + beta*I)

    for i in range(100):
        efermi = 0.5*(elow+ehigh)
        nefermi = trace(alpha*(efermi*I-F)+ beta*I)
        print elow,ehigh,nelow,Ne,nehigh
        if abs(Ne-nefermi) < tol: break
        if nefermi < Ne:
            elow = efermi
            nelow = nefermi
        elif nefermi > Ne:
            ehigh = efermi
            nehigh = nefermi
    alpha = min(beta/(emax-efermi),(1-beta)/(efermi-emin))
    return alpha*(efermi*I-F)+beta*I



def DMP(F,S,Ne,Method=0,MaxIter=50,ErrorLimit=1e-12):
    # Density Matrix Purification Methods
    # 0 -> Trace correcting purification (default)
    # 1 -> Trace resetting
    # 2 -> McWeeny purification
    # 3 -> Canonical purification
    logging.warning("2/20/07: DMP is now deprecated. Use DMP.TCP instead")
    methods = ['TCP','TRS','MCW','PM']

    # Step 1: Orthogonalize the Fock matrix:
    X = SymOrth(S)
    F = simx(F,X)

    # Step 2: Initialize the density matrix:
    emin,emax = gershgorin_minmax(F)
    N = F.shape[0]
    I = identity(N,'d')
    if Method == 0 or Method == 1:
        D = (emax*I-F)/(emax-emin)
    elif Method == 2:
        D = Dinit_mcw(F,Ne)
    elif Method == 3:
        efermi = trace(F)/N
        beta = Ne/float(N)
        alpha = min(Ne/(emax-efermi),(N-Ne)/(efermi-emin))/float(N)
        D = alpha*(efermi*I-F) + beta*I
    else:
        raise Exception("Unknown method %d" % Method)

    Dsumold = sum(sum(D))
    # Step 3: Iterate on DM updates:
    for iter in range(MaxIter):
        Ne_curr = trace(D)
        D2 = matrixmultiply(D,D)
        if Method == 0:
            if Ne_curr < Ne:
                D = 2.0*D-D2
            else:
                D = D2
            if abs(Ne_curr-Ne) < ErrorLimit: break
        elif Method == 1:
            Df = matrixmultiply(D2,4*D-3*D2)
            trf = trace(Df)
            Dp = I-D
            Dp2 = matrixmultiply(Dp,Dp)
            Dg = matrixmultiply(D2,Dp2)
            trg = trace(Dg)
            gamma = (Ne-trf)/trg
            if gamma > 2:
                D = 2*D-D2
            elif gamma < 0:
                D = D2
            else:
                D = Df-gamma*Dg
            if abs(Ne_curr-Ne) < ErrorLimit: break
        elif Method == 2:
            D = 3*D2-2*matrixmultiply(D,D2)
            if abs(Ne_curr-Ne) < ErrorLimit: break
        elif Method == 3:
            D3 = matrixmultiply(D,D2)
            cn = trace(D2-D3)/trace(D-D2)
            if cn < 0.5:
                D = ((1.0-2.0*cn)*D+(1.0+cn)*D2-D3)/(1.0-cn)
            else:
                D = ((1+cn)*D2-D3)/cn
            Dsum = sum(sum(D))
            if Dsum-Dsumold  < ErrorLimit: break
            Dsumold = Dsum
    else: print "DMP: Warning MaxIters reached"
    #print "%s converged in %d iters" % (methods[Method],iter)
    D = simx(D,X,'T')
    return D
        
    
