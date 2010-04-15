"""\
Implement Hartree-Fock Greens Function for band energy corrections.
From Szabo/Ostland Chapter 7.
Correcting orbital 5, HF Eorb = -0.620213
-> Eorb = -0.527220 -0.538650
Correcting orbital 6, HF Eorb = -0.571950
-> Eorb = -0.603409 -0.600940
Correcting orbital 7, HF Eorb = -0.571950
-> Eorb = -0.603409 -0.600939
Correcting orbital 8, HF Eorb = 0.119814
-> Eorb = 0.093528 0.095544
Correcting orbital 9, HF Eorb = 0.119814
-> Eorb = 0.093528 0.095544
Correcting orbital 10, HF Eorb = 0.539609
-> Eorb = 0.488557 0.490867


 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from hartree_fock import scf
from Ints import getbasis, getints
from CI import TransformInts
from PyQuante.cints import ijkl2intindex
from NumWrap import zeros
from NumWrap import det

class Sigma2:
    def __init__(self,orbe,orbs,ints,norb,nocc):
        self.e0 = orbe
        self.orbs = orbs
        self.norb = norb
        self.nocc = nocc
        print "Transforming 2e ints"
        self.moints = TransformInts(ints,orbs)
        print "...done"
        return

    def eval0(self,i,E):
        # just do the simple approximation, S/A eqs 7.44-7.46
        term = 0.
        occs = range(self.nocc)
        virts = range(self.nocc,self.norb)
        for a in occs:
            for r in virts:
                for s in virts:
                    iras = self.moints[ijkl2intindex(i,r,a,s)]
                    isar = self.moints[ijkl2intindex(i,s,a,r)]
                    term += iras*(2*iras-isar)/(
                        E+self.e0[a]-self.e0[r]-self.e0[s])
        for a in occs:
            for b in occs:
                for r in virts:
                    iabr = self.moints[ijkl2intindex(i,a,b,r)]
                    ibar = self.moints[ijkl2intindex(i,b,a,r)]
                    term += iabr*(2*iabr-ibar)/(
                        E+self.e0[r]-self.e0[a]-self.e0[b])
        return term

    def eval(self,i):
        del0 = self.eval0(i,self.e0[i])
        dele = del0
        delold = dele
        for j in xrange(60):
            dele = self.eval0(i,self.e0[i]+dele)
            if abs(dele-delold)<1e-5: break
            delold = dele
        else:
            print "Warning: maxiter reached in Sigma2::eval for %d" % i
        return del0,dele

    def fulleval(self,E):
        # This routine doesn't work very well.
        g = zeros((self.norb,self.norb),'d')
        orbs = range(self.norb)
        occs = range(self.nocc)
        virts = range(self.nocc,self.norb)
        term = 0.
        for i in orbs:
            for j in orbs:
                for a in occs:
                    for r in virts:
                        for s in virts:
                            iras = self.moints[ijkl2intindex(i,r,a,s)]
                            jras = self.moints[ijkl2intindex(j,r,a,s)]
                            jsar = self.moints[ijkl2intindex(j,s,a,r)]
                            term += iras*(2*jras-jsar)/(
                                E+self.e0[a]-self.e0[r]-self.e0[s])
                for a in occs:
                    for b in occs:
                        for r in virts:
                            iabr = self.moints[ijkl2intindex(i,a,b,r)]
                            jabr = self.moints[ijkl2intindex(j,a,b,r)]
                            jbar = self.moints[ijkl2intindex(j,b,a,r)]
                            term += iabr*(2*jabr-jbar)/(
                                E+self.e0[r]-self.e0[a]-self.e0[b])
                g[i,j] = -term
        for i in orbs:
            g[i,i] += E - self.e0[i]
        return det(g)

def print_orbe(orbe,nclosed,nvirt,tag=None):
    if tag: print tag
    print "   N  Occ    Energy"
    print "---------------------"
    for i in xrange(nclosed):
        print "%4d  1   %10.4f" % (i+1,orbe[i])
    nv = min(nvirt,10) 
    for i in xrange(nclosed, nclosed+nv):
        print "%4d  0   %10.4f" % (i+1,orbe[i])
    return
    

def HFGF(atoms,charge=0):
    nclosed,nopen = atoms.get_closedopen()
    if nopen: raise Exception("HFGF only works for closed shell cases")
    bfs = getbasis(atoms)
    nvirt = len(bfs)-nclosed
    S,h,Ints = getints(bfs,atoms)
    hf_energy,hf_orbe,hf_orbs = scf(atoms,S,h,Ints,charge)
    print_orbe(hf_orbe,nclosed,nvirt)

    sigma = Sigma2(hf_orbe,hf_orbs,Ints,len(bfs),nclosed)
    
    for i in [4,5,6,7,8,9]:
        print "Correcting orbital %d, HF Eorb = %f" % (i+1,hf_orbe[i])
        del0,dele = sigma.eval(i)
        print "-> Eorb = %f %f" % (hf_orbe[i]+del0,hf_orbe[i]+dele)

    return

if __name__ == '__main__':
    from Molecule import Molecule
    x = 1.0783/2/0.52918
    y = 0.75/2/0.52918
    n2 = Molecule('N2',atomlist=[(7,(x,0,0)),(7,(-x,0,0))])
    h2 = Molecule('H2',atomlist=[(1,(y,0,0)),(1,(-1,0,0))])
    HFGF(n2)
    #HFGF(h2)
