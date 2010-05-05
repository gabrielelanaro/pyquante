"""\
 Implementation of Head-Gordon & Pople's scheme for electron repulsion
  integrals (ref), which, in turn, derives from Saika and Obarra's scheme.

 Routines:
 hrr performs the horizontal recursion relationships
 vrr performs the vertical recursion relationship

 The routines in the accompanying chgp module have the same functions, but
 are written in C to be faster.

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""

from math import sqrt,pi,exp
#from PyQuante.cints import Fgamma
from pyints import gaussian_product_center,Fgamma
#from PyQuante.cints import vrr

def contr_hrr((xa,ya,za),norma,(la,ma,na),aexps,acoefs,
              (xb,yb,zb),normb,(lb,mb,nb),bexps,bcoefs,
              (xc,yc,zc),normc,(lc,mc,nc),cexps,ccoefs,
              (xd,yd,zd),normd,(ld,md,nd),dexps,dcoefs):
    if lb > 0:
        return (contr_hrr((xa,ya,za),norma,(la+1,ma,na),aexps,acoefs,
                    (xb,yb,zb),normb,(lb-1,mb,nb),bexps,bcoefs,
                    (xc,yc,zc),normc,(lc,mc,nc),cexps,ccoefs,
                    (xd,yd,zd),normd,(ld,md,nd),dexps,dcoefs)
                + (xa-xb)*contr_hrr((xa,ya,za),norma,(la,ma,na),aexps,acoefs,
                              (xb,yb,zb),normb,(lb-1,mb,nb),bexps,bcoefs,
                              (xc,yc,zc),normc,(lc,mc,nc),cexps,ccoefs,
                              (xd,yd,zd),normd,(ld,md,nd),dexps,dcoefs)
                )
    elif mb > 0:
        return (contr_hrr((xa,ya,za),norma,(la,ma+1,na),aexps,acoefs,
                    (xb,yb,zb),normb,(lb,mb-1,nb),bexps,bcoefs,
                    (xc,yc,zc),normc,(lc,mc,nc),cexps,ccoefs,
                    (xd,yd,zd),normd,(ld,md,nd),dexps,dcoefs)
                + (ya-yb)*contr_hrr((xa,ya,za),norma,(la,ma,na),aexps,acoefs,
                              (xb,yb,zb),normb,(lb,mb-1,nb),bexps,bcoefs,
                              (xc,yc,zc),normc,(lc,mc,nc),cexps,ccoefs,
                              (xd,yd,zd),normd,(ld,md,nd),dexps,dcoefs)
                )
    elif nb > 0:
        return (contr_hrr((xa,ya,za),norma,(la,ma,na+1),aexps,acoefs,
                    (xb,yb,zb),normb,(lb,mb,nb-1),bexps,bcoefs,
                    (xc,yc,zc),normc,(lc,mc,nc),cexps,ccoefs,
                    (xd,yd,zd),normd,(ld,md,nd),dexps,dcoefs)
                + (za-zb)*contr_hrr((xa,ya,za),norma,(la,ma,na),aexps,acoefs,
                              (xb,yb,zb),normb,(lb,mb,nb-1),bexps,bcoefs,
                              (xc,yc,zc),normc,(lc,mc,nc),cexps,ccoefs,
                              (xd,yd,zd),normd,(ld,md,nd),dexps,dcoefs)
                )
    elif ld > 0:
        return (contr_hrr((xa,ya,za),norma,(la,ma,na),aexps,acoefs,
                    (xb,yb,zb),normb,(lb,mb,nb),bexps,bcoefs,
                    (xc,yc,zc),normc,(lc+1,mc,nc),cexps,ccoefs,
                    (xd,yd,zd),normd,(ld-1,md,nd),dexps,dcoefs)
                + (xc-xd)*contr_hrr((xa,ya,za),norma,(la,ma,na),aexps,acoefs,
                              (xb,yb,zb),normb,(lb,mb,nb),bexps,bcoefs,
                              (xc,yc,zc),normc,(lc,mc,nc),cexps,ccoefs,
                              (xd,yd,zd),normd,(ld-1,md,nd),dexps,dcoefs)
                )
    elif md > 0:
        return (contr_hrr((xa,ya,za),norma,(la,ma,na),aexps,acoefs,
                    (xb,yb,zb),normb,(lb,mb,nb),bexps,bcoefs,
                    (xc,yc,zc),normc,(lc,mc+1,nc),cexps,ccoefs,
                    (xd,yd,zd),normd,(ld,md-1,nd),dexps,dcoefs)
                + (yc-yd)*contr_hrr((xa,ya,za),norma,(la,ma,na),aexps,acoefs,
                              (xb,yb,zb),normb,(lb,mb,nb),bexps,bcoefs,
                              (xc,yc,zc),normc,(lc,mc,nc),cexps,ccoefs,
                              (xd,yd,zd),normd,(ld,md-1,nd),dexps,dcoefs)
                )
    elif nd > 0:
        return (contr_hrr((xa,ya,za),norma,(la,ma,na),aexps,acoefs,
                    (xb,yb,zb),normb,(lb,mb,nb),bexps,bcoefs,
                    (xc,yc,zc),normc,(lc,mc,nc+1),cexps,ccoefs,
                    (xd,yd,zd),normd,(ld,md,nd-1),dexps,dcoefs)
                + (zc-zd)*contr_hrr((xa,ya,za),norma,(la,ma,na),aexps,acoefs,
                              (xb,yb,zb),normb,(lb,mb,nb),bexps,bcoefs,
                              (xc,yc,zc),normc,(lc,mc,nc),cexps,ccoefs,
                              (xd,yd,zd),normd,(ld,md,nd-1),dexps,dcoefs)
                )
    return contr_vrr((xa,ya,za),norma,(la,ma,na),aexps,acoefs,
                     (xb,yb,zb),normb,bexps,bcoefs,
                     (xc,yc,zc),normc,(lc,mc,nc),cexps,ccoefs,
                     (xd,yd,zd),normd,dexps,dcoefs)

def contr_hrr_iter((xa,ya,za),norma,(la,ma,na),aexps,acoefs,
                   (xb,yb,zb),normb,(lb,mb,nb),bexps,bcoefs,
                   (xc,yc,zc),normc,(lc,mc,nc),cexps,ccoefs,
                   (xd,yd,zd),normd,(ld,md,nd),dexps,dcoefs):

    # Non-recursive version, which I implemented in the false hope
    #  that it would reduce the number of calls to vrr. Kept as a
    #  valuable lesson.

    # Precompute all of the required (i0|j0) terms:
    # this is a dictionary pretending to be a 12-d array, which is
    # not the prettiest data structure in the world, and will probably
    # be a challenge when I convert to C.
    hrr_terms = {}
    for i in xrange(la+lb+1):
        for j in xrange(ma+mb+1):
            for k in xrange(na+nb+1):
                for l in xrange(lc+ld+1):
                    for m in xrange(mc+md+1):
                        for n in xrange(nc+nd+1):
                            hrr_terms[i,j,k,0,0,0,l,m,n,0,0,0] = (
                                contr_vrr((xa,ya,za),norma,(i,j,k),aexps,acoefs,
                                          (xb,yb,zb),normb,bexps,bcoefs,
                                          (xc,yc,zc),normc,(l,m,n),cexps,ccoefs,
                                          (xd,yd,zd),normd,dexps,dcoefs)
                                )
    # At this point we have all of the integrals (i0|j0).
    # We now need to use the hrrs to build up all (ab|cd).
    for i in xrange(1,lb+1):
        for j in xrange(ma+mb+1):
            for k in xrange(na+nb+1):
                for l in xrange(lc+ld+1):
                    for m in xrange(mc+md+1):
                        for n in xrange(nc+nd+1):
                            hrr_terms[la+lb-i,j,k,i,0,0,l,m,n,0,0,0] = (
                                hrr_terms[la+lb-i+1,j,k,i-1,0,0,l,m,n,0,0,0]
                                + (xa-xb)*
                                hrr_terms[la+lb-i,j,k,i-1,0,0,l,m,n,0,0,0]
                                )
    
    for i in xrange(1,mb+1):
        for j in xrange(na+nb+1):
            for k in xrange(lc+ld+1):
                for l in xrange(mc+md+1):
                    for m in xrange(nc+nd+1):
                        hrr_terms[la,ma+mb-i,j,lb,i,0,k,l,m,0,0,0] = (
                            hrr_terms[la,ma+mb-i+1,j,lb,i-1,0,k,l,m,0,0,0]
                            + (ya-yb)*
                            hrr_terms[la,ma+mb-i,j,lb,i-1,0,k,l,m,0,0,0]
                            )

    for i in xrange(1,nb+1):
        for j in xrange(lc+ld+1):
            for k in xrange(mc+md+1):
                for l in xrange(nc+nd+1):
                    hrr_terms[la,ma,na+nb-i,lb,mb,i,j,k,l,0,0,0] = (
                        hrr_terms[la,ma,na+nb-i+1,lb,mb,i-1,j,k,l,0,0,0]
                        + (za-zb)*
                        hrr_terms[la,ma,na+nb-i,lb,mb,i-1,j,k,l,0,0,0]
                        )

    for i in xrange(1,ld+1):
        for j in xrange(mc+md+1):
            for k in xrange(nc+nd+1):
                hrr_terms[la,ma,na,lb,mb,nb,lc+ld-i,j,k,i,0,0] = (
                    hrr_terms[la,ma,na,lb,mb,nb,lc+ld-i+1,j,k,i-1,0,0]
                    + (xc-xd)*
                    hrr_terms[la,ma,na,lb,mb,nb,lc+ld-i,j,k,i-1,0,0]
                    )
        
    for i in xrange(1,md+1):
        for j in xrange(nc+nd+1):
            hrr_terms[la,ma,na,lb,mb,nb,lc,mc+md-i,j,ld,i,0] = (
                hrr_terms[la,ma,na,lb,mb,nb,lc,mc+md-i+1,j,ld,i-1,0]
                + (yc-yd)*
                hrr_terms[la,ma,na,lb,mb,nb,lc,mc+md-i,j,ld,i-1,0]
                )

    for i in xrange(1,nd+1):
        hrr_terms[la,ma,na,lb,mb,nb,lc,mc,nc+nd-i,ld,md,i] = (
            hrr_terms[la,ma,na,lb,mb,nb,lc,mc,nc+nd-i+1,ld,md,i-1]
            + (zc-zd)*
            hrr_terms[la,ma,na,lb,mb,nb,lc,mc,nc+nd-i,ld,md,i-1]
            )

    # Done; return the relevant value:
    return hrr_terms[la,ma,na,lb,mb,nb,lc,mc,nc,ld,md,nd]


def contr_vrr((xa,ya,za),norma,(la,ma,na),aexps,acoefs,
              (xb,yb,zb),normb,bexps,bcoefs,
              (xc,yc,zc),normc,(lc,mc,nc),cexps,ccoefs,
              (xd,yd,zd),normd,dexps,dcoefs):
    val = 0.
    for i in xrange(len(aexps)):
        for j in xrange(len(bexps)):
            for k in xrange(len(cexps)):
                for l in xrange(len(dexps)):
                    val = val + acoefs[i]*bcoefs[j]*ccoefs[k]*dcoefs[l]\
                          *vrr((xa,ya,za),norma[i],(la,ma,na),aexps[i],
                               (xb,yb,zb),normb[j],bexps[j],
                               (xc,yc,zc),normc[k],(lc,mc,nc),cexps[k],
                               (xd,yd,zd),normd[l],dexps[l],0)
    return val

def hrr((xa,ya,za),norma,(la,ma,na),alphaa,
        (xb,yb,zb),normb,(lb,mb,nb),alphab,
        (xc,yc,zc),normc,(lc,mc,nc),alphac,
        (xd,yd,zd),normd,(ld,md,nd),alphad):

    if lb > 0:
        return (hrr((xa,ya,za),norma,(la+1,ma,na),alphaa,
                    (xb,yb,zb),normb,(lb-1,mb,nb),alphab,
                    (xc,yc,zc),normc,(lc,mc,nc),alphac,
                    (xd,yd,zd),normd,(ld,md,nd),alphad)
                + (xa-xb)*hrr((xa,ya,za),norma,(la,ma,na),alphaa,
                              (xb,yb,zb),normb,(lb-1,mb,nb),alphab,
                              (xc,yc,zc),normc,(lc,mc,nc),alphac,
                              (xd,yd,zd),normd,(ld,md,nd),alphad)
                )
    elif mb > 0:
        return (hrr((xa,ya,za),norma,(la,ma+1,na),alphaa,
                    (xb,yb,zb),normb,(lb,mb-1,nb),alphab,
                    (xc,yc,zc),normc,(lc,mc,nc),alphac,
                    (xd,yd,zd),normd,(ld,md,nd),alphad)
                + (ya-yb)*hrr((xa,ya,za),norma,(la,ma,na),alphaa,
                              (xb,yb,zb),normb,(lb,mb-1,nb),alphab,
                              (xc,yc,zc),normc,(lc,mc,nc),alphac,
                              (xd,yd,zd),normd,(ld,md,nd),alphad)
                )
    elif nb > 0:
        return (hrr((xa,ya,za),norma,(la,ma,na+1),alphaa,
                    (xb,yb,zb),normb,(lb,mb,nb-1),alphab,
                    (xc,yc,zc),normc,(lc,mc,nc),alphac,
                    (xd,yd,zd),normd,(ld,md,nd),alphad)
                + (za-zb)*hrr((xa,ya,za),norma,(la,ma,na),alphaa,
                              (xb,yb,zb),normb,(lb,mb,nb-1),alphab,
                              (xc,yc,zc),normc,(lc,mc,nc),alphac,
                              (xd,yd,zd),normd,(ld,md,nd),alphad)
                )
    elif ld > 0:
        return (hrr((xa,ya,za),norma,(la,ma,na),alphaa,
                    (xb,yb,zb),normb,(lb,mb,nb),alphab,
                    (xc,yc,zc),normc,(lc+1,mc,nc),alphac,
                    (xd,yd,zd),normd,(ld-1,md,nd),alphad)
                + (xc-xd)*hrr((xa,ya,za),norma,(la,ma,na),alphaa,
                              (xb,yb,zb),normb,(lb,mb,nb),alphab,
                              (xc,yc,zc),normc,(lc,mc,nc),alphac,
                              (xd,yd,zd),normd,(ld-1,md,nd),alphad)
                )
    elif md > 0:
        return (hrr((xa,ya,za),norma,(la,ma,na),alphaa,
                    (xb,yb,zb),normb,(lb,mb,nb),alphab,
                    (xc,yc,zc),normc,(lc,mc+1,nc),alphac,
                    (xd,yd,zd),normd,(ld,md-1,nd),alphad)
                + (yc-yd)*hrr((xa,ya,za),norma,(la,ma,na),alphaa,
                              (xb,yb,zb),normb,(lb,mb,nb),alphab,
                              (xc,yc,zc),normc,(lc,mc,nc),alphac,
                              (xd,yd,zd),normd,(ld,md-1,nd),alphad)
                )
    elif nd > 0:
        return (hrr((xa,ya,za),norma,(la,ma,na),alphaa,
                    (xb,yb,zb),normb,(lb,mb,nb),alphab,
                    (xc,yc,zc),normc,(lc,mc,nc+1),alphac,
                    (xd,yd,zd),normd,(ld,md,nd-1),alphad)
                + (zc-zd)*hrr((xa,ya,za),norma,(la,ma,na),alphaa,
                              (xb,yb,zb),normb,(lb,mb,nb),alphab,
                              (xc,yc,zc),normc,(lc,mc,nc),alphac,
                              (xd,yd,zd),normd,(ld,md,nd-1),alphad)
                )
    return vrr((xa,ya,za),norma,(la,ma,na),alphaa,
               (xb,yb,zb),normb,alphab,
               (xc,yc,zc),normc,(lc,mc,nc),alphac,
               (xd,yd,zd),normd,alphad,0)

def vrr_recursive((xa,ya,za),norma,(la,ma,na),alphaa,
        (xb,yb,zb),normb,alphab,
        (xc,yc,zc),normc,(lc,mc,nc),alphac,
        (xd,yd,zd),normd,alphad,m):
    "Old VRR code, which is called recursively"

    px,py,pz = gaussian_product_center(alphaa,(xa,ya,za),alphab,(xb,yb,zb))
    qx,qy,qz = gaussian_product_center(alphac,(xc,yc,zc),alphad,(xd,yd,zd))
    zeta,eta = alphaa+alphab,alphac+alphad
    wx,wy,wz = gaussian_product_center(zeta,(px,py,pz),eta,(qx,qy,qz))

    if nc:
        val = (qz-zc)*vrr((xa,ya,za),norma,(la,ma,na),alphaa,
                          (xb,yb,zb),normb,alphab,
                          (xc,yc,zc),normc,(lc,mc,nc-1),alphac,
                          (xd,yd,zd),normd,alphad,m) \
              + (wz-qz)*vrr((xa,ya,za),norma,(la,ma,na),alphaa,
                            (xb,yb,zb),normb,alphab,
                            (xc,yc,zc),normc,(lc,mc,nc-1),alphac,
                            (xd,yd,zd),normd,alphad,m+1)
        if nc > 1:
            val = val +\
                  0.5*(nc-1)/eta*(vrr((xa,ya,za),norma,(la,ma,na),alphaa,
                                   (xb,yb,zb),normb,alphab,
                                   (xc,yc,zc),normc,(lc,mc,nc-2),alphac,
                                   (xd,yd,zd),normd,alphad,m)
                               -zeta/(zeta+eta)* 
                               vrr((xa,ya,za),norma,(la,ma,na),alphaa,
                                   (xb,yb,zb),normb,alphab,
                                   (xc,yc,zc),normc,(lc,mc,nc-2),alphac,
                                   (xd,yd,zd),normd,alphad,m+1) )
        if na:
            val = val +\
                  0.5*na/(zeta+eta)*vrr((xa,ya,za),norma,(la,ma,na-1),alphaa,
                                   (xb,yb,zb),normb,alphab,
                                   (xc,yc,zc),normc,(lc,mc,nc-1),alphac,
                                   (xd,yd,zd),normd,alphad,m+1)
        return val

    elif mc:
        val = (qy-yc)*vrr((xa,ya,za),norma,(la,ma,na),alphaa,
                          (xb,yb,zb),normb,alphab,
                          (xc,yc,zc),normc,(lc,mc-1,nc),alphac,
                          (xd,yd,zd),normd,alphad,m) \
              + (wy-qy)*vrr((xa,ya,za),norma,(la,ma,na),alphaa,
                            (xb,yb,zb),normb,alphab,
                            (xc,yc,zc),normc,(lc,mc-1,nc),alphac,
                            (xd,yd,zd),normd,alphad,m+1)
        if mc > 1:
            val = val +\
                  0.5*(mc-1)/eta*(vrr((xa,ya,za),norma,(la,ma,na),alphaa,
                                   (xb,yb,zb),normb,alphab,
                                   (xc,yc,zc),normc,(lc,mc-2,nc),alphac,
                                   (xd,yd,zd),normd,alphad,m)
                               -zeta/(zeta+eta)* 
                               vrr((xa,ya,za),norma,(la,ma,na),alphaa,
                                   (xb,yb,zb),normb,alphab,
                                   (xc,yc,zc),normc,(lc,mc-2,nc),alphac,
                                   (xd,yd,zd),normd,alphad,m+1) )
        if ma:
            val = val +\
                  0.5*ma/(zeta+eta)*vrr((xa,ya,za),norma,(la,ma-1,na),alphaa,
                                   (xb,yb,zb),normb,alphab,
                                   (xc,yc,zc),normc,(lc,mc-1,nc),alphac,
                                   (xd,yd,zd),normd,alphad,m+1)
        return val
    elif lc:
        val = (qx-xc)*vrr((xa,ya,za),norma,(la,ma,na),alphaa,
                          (xb,yb,zb),normb,alphab,
                          (xc,yc,zc),normc,(lc-1,mc,nc),alphac,
                          (xd,yd,zd),normd,alphad,m) \
              + (wx-qx)*vrr((xa,ya,za),norma,(la,ma,na),alphaa,
                            (xb,yb,zb),normb,alphab,
                            (xc,yc,zc),normc,(lc-1,mc,nc),alphac,
                            (xd,yd,zd),normd,alphad,m+1)
        if lc > 1:
            val = val +\
                  0.5*(lc-1)/eta*(vrr((xa,ya,za),norma,(la,ma,na),alphaa,
                                   (xb,yb,zb),normb,alphab,
                                   (xc,yc,zc),normc,(lc-2,mc,nc),alphac,
                                   (xd,yd,zd),normd,alphad,m)
                               -zeta/(zeta+eta)* 
                               vrr((xa,ya,za),norma,(la,ma,na),alphaa,
                                   (xb,yb,zb),normb,alphab,
                                   (xc,yc,zc),normc,(lc-2,mc,nc),alphac,
                                   (xd,yd,zd),normd,alphad,m+1) )
        if la:
            val = val +\
                  0.5*la/(zeta+eta)*vrr((xa,ya,za),norma,(la-1,ma,na),alphaa,
                                   (xb,yb,zb),normb,alphab,
                                   (xc,yc,zc),normc,(lc-1,mc,nc),alphac,
                                   (xd,yd,zd),normd,alphad,m+1)
        return val

    elif na:
        val = (pz-za)*vrr((xa,ya,za),norma,(la,ma,na-1),alphaa,
                               (xb,yb,zb),normb,alphab,
                               (xc,yc,zc),normc,(lc,mc,nc),alphac,
                               (xd,yd,zd),normd,alphad,m) \
                  + (wz-pz)*vrr((xa,ya,za),norma,(la,ma,na-1),alphaa,
                                  (xb,yb,zb),normb,alphab,
                                  (xc,yc,zc),normc,(lc,mc,nc),alphac,
                                  (xd,yd,zd),normd,alphad,m+1)
        
        if na > 1:
            val =  val + \
                  0.5*(na-1)/zeta*(vrr((xa,ya,za),norma,(la,ma,na-2),alphaa,
                                     (xb,yb,zb),normb,alphab,
                                     (xc,yc,zc),normc,(lc,mc,nc),alphac,
                                     (xd,yd,zd),normd,alphad,m)
                               -eta/(zeta+eta)* 
                               vrr((xa,ya,za),norma,(la,ma,na-2),alphaa,
                                     (xb,yb,zb),normb,alphab,
                                     (xc,yc,zc),normc,(lc,mc,nc),alphac,
                                     (xd,yd,zd),normd,alphad,m+1) )
        return val
    elif ma:
        val = (py-ya)*vrr((xa,ya,za),norma,(la,ma-1,na),alphaa,
                               (xb,yb,zb),normb,alphab,
                               (xc,yc,zc),normc,(lc,mc,nc),alphac,
                               (xd,yd,zd),normd,alphad,m) \
                  + (wy-py)*vrr((xa,ya,za),norma,(la,ma-1,na),alphaa,
                                  (xb,yb,zb),normb,alphab,
                                  (xc,yc,zc),normc,(lc,mc,nc),alphac,
                                  (xd,yd,zd),normd,alphad,m+1)
        
        if ma > 1:
            val =  val + \
                  0.5*(ma-1)/zeta*(vrr((xa,ya,za),norma,(la,ma-2,na),alphaa,
                                     (xb,yb,zb),normb,alphab,
                                     (xc,yc,zc),normc,(lc,mc,nc),alphac,
                                     (xd,yd,zd),normd,alphad,m)
                               -eta/(zeta+eta)* 
                               vrr((xa,ya,za),norma,(la,ma-2,na),alphaa,
                                     (xb,yb,zb),normb,alphab,
                                     (xc,yc,zc),normc,(lc,mc,nc),alphac,
                                     (xd,yd,zd),normd,alphad,m+1) )
        return val
    elif la:
        val = (px-xa)*vrr((xa,ya,za),norma,(la-1,ma,na),alphaa,
                               (xb,yb,zb),normb,alphab,
                               (xc,yc,zc),normc,(lc,mc,nc),alphac,
                               (xd,yd,zd),normd,alphad,m) \
                  + (wx-px)*vrr((xa,ya,za),norma,(la-1,ma,na),alphaa,
                                  (xb,yb,zb),normb,alphab,
                                  (xc,yc,zc),normc,(lc,mc,nc),alphac,
                                  (xd,yd,zd),normd,alphad,m+1)
        
        if la > 1:
            val =  val + \
                  0.5*(la-1)/zeta*(vrr((xa,ya,za),norma,(la-2,ma,na),alphaa,
                                     (xb,yb,zb),normb,alphab,
                                     (xc,yc,zc),normc,(lc,mc,nc),alphac,
                                     (xd,yd,zd),normd,alphad,m)
                               -eta/(zeta+eta)* 
                               vrr((xa,ya,za),norma,(la-2,ma,na),alphaa,
                                     (xb,yb,zb),normb,alphab,
                                     (xc,yc,zc),normc,(lc,mc,nc),alphac,
                                     (xd,yd,zd),normd,alphad,m+1) )
        return val

    rab2 = pow(xa-xb,2) + pow(ya-yb,2) + pow(za-zb,2)
    Kab = sqrt(2)*pow(pi,1.25)/(alphaa+alphab)\
          *exp(-alphaa*alphab/(alphaa+alphab)*rab2)
    rcd2 = pow(xc-xd,2) + pow(yc-yd,2) + pow(zc-zd,2)
    Kcd = sqrt(2)*pow(pi,1.25)/(alphac+alphad)\
          *exp(-alphac*alphad/(alphac+alphad)*rcd2)
    rpq2 = pow(px-qx,2) + pow(py-qy,2) + pow(pz-qz,2)
    T = zeta*eta/(zeta+eta)*rpq2
    val = norma*normb*normc*normd*Kab*Kcd/sqrt(zeta+eta)*Fgamma(m,T)
    return val

def vrr((xa,ya,za),norma,(la,ma,na),alphaa,
        (xb,yb,zb),normb,alphab,
        (xc,yc,zc),normc,(lc,mc,nc),alphac,
        (xd,yd,zd),normd,alphad,M):

    px,py,pz = gaussian_product_center(alphaa,(xa,ya,za),alphab,(xb,yb,zb))
    qx,qy,qz = gaussian_product_center(alphac,(xc,yc,zc),alphad,(xd,yd,zd))
    zeta,eta = float(alphaa+alphab),float(alphac+alphad)
    wx,wy,wz = gaussian_product_center(zeta,(px,py,pz),eta,(qx,qy,qz))

    rab2 = pow(xa-xb,2) + pow(ya-yb,2) + pow(za-zb,2)
    Kab = sqrt(2)*pow(pi,1.25)/(alphaa+alphab)\
          *exp(-alphaa*alphab/(alphaa+alphab)*rab2)
    rcd2 = pow(xc-xd,2) + pow(yc-yd,2) + pow(zc-zd,2)
    Kcd = sqrt(2)*pow(pi,1.25)/(alphac+alphad)\
          *exp(-alphac*alphad/(alphac+alphad)*rcd2)
    rpq2 = pow(px-qx,2) + pow(py-qy,2) + pow(pz-qz,2)
    T = zeta*eta/(zeta+eta)*rpq2

    mtot = la+ma+na+lc+mc+nc+M

    Fgterms = [0]*(mtot+1)
    Fgterms[mtot] = Fgamma(mtot,T)
    for im in xrange(mtot-1,-1,-1):
        Fgterms[im]=(2.*T*Fgterms[im+1]+exp(-T))/(2.*im+1)

    # Store the vrr values as a 7 dimensional array
    # vrr_terms[la,ma,na,lc,mc,nc,m]
    vrr_terms = {}
    for im in xrange(mtot+1):
        vrr_terms[0,0,0,0,0,0,im] = (
            #norma*normb*normc*normd*Kab*Kcd/sqrt(zeta+eta)*Fgamma(im,T)
            norma*normb*normc*normd*Kab*Kcd/sqrt(zeta+eta)*Fgterms[im]
            )

    for i in xrange(la):
        for im in xrange(mtot-i):
            vrr_terms[i+1,0,0, 0,0,0, im] = (
                (px-xa)*vrr_terms[i,0,0, 0,0,0, im]
                + (wx-px)*vrr_terms[i,0,0, 0,0,0, im+1]
                )
            if i:
                vrr_terms[i+1,0,0, 0,0,0, im] += (
                    i/2./zeta*( vrr_terms[i-1,0,0, 0,0,0, im]
                               - eta/(zeta+eta)*vrr_terms[i-1,0,0, 0,0,0, im+1]
                               ))

    for j in xrange(ma):
        for i in xrange(la+1):
            for im in xrange(mtot-i-j):
                vrr_terms[i,j+1,0, 0,0,0, im] = (
                    (py-ya)*vrr_terms[i,j,0, 0,0,0, im]
                    + (wy-py)*vrr_terms[i,j,0, 0,0,0, im+1]
                    )
                if j:
                    vrr_terms[i,j+1,0, 0,0,0, im] += (
                        j/2./zeta*(vrr_terms[i,j-1,0, 0,0,0, im]
                                  - eta/(zeta+eta)
                                  *vrr_terms[i,j-1,0, 0,0,0, im+1]
                                  ))


    for k in xrange(na):
        for j in xrange(ma+1):
            for i in xrange(la+1):
                for im in xrange(mtot-i-j-k):
                    vrr_terms[i,j,k+1, 0,0,0, im] = (
                        (pz-za)*vrr_terms[i,j,k, 0,0,0, im]
                        + (wz-pz)*vrr_terms[i,j,k, 0,0,0, im+1]
                        )
                    if k:
                        vrr_terms[i,j,k+1, 0,0,0, im] += (
                            k/2./zeta*(vrr_terms[i,j,k-1, 0,0,0, im]
                                      - eta/(zeta+eta)
                                      *vrr_terms[i,j,k-1, 0,0,0, im+1]
                                      ))

    for q in xrange(lc):
        for k in xrange(na+1):
            for j in xrange(ma+1):
                for i in xrange(la+1):
                    for im in xrange(mtot-i-j-k-q):
                        vrr_terms[i,j,k, q+1,0,0, im] = (
                            (qx-xc)*vrr_terms[i,j,k, q,0,0, im]
                            + (wx-qx)*vrr_terms[i,j,k, q,0,0, im+1]
                            )
                        if q:
                            vrr_terms[i,j,k, q+1,0,0, im] += (
                                q/2./eta*(vrr_terms[i,j,k, q-1,0,0, im]
                                         - zeta/(zeta+eta)
                                         *vrr_terms[i,j,k, q-1,0,0, im+1]
                                         ))
                        if i:
                            vrr_terms[i,j,k, q+1,0,0, im] += (
                                i/2./(zeta+eta)*vrr_terms[i-1,j,k, q,0,0, im+1]
                                )

    for r in xrange(mc):
        for q in xrange(lc+1):
            for k in xrange(na+1):
                for j in xrange(ma+1):
                    for i in xrange(la+1):
                        for im in xrange(mtot-i-j-k-q-r):
                            vrr_terms[i,j,k, q,r+1,0, im] = (
                                (qy-yc)*vrr_terms[i,j,k, q,r,0, im]
                                + (wy-qy)*vrr_terms[i,j,k, q,r,0, im+1]
                                )
                            if r:
                                vrr_terms[i,j,k, q,r+1,0, im] += (
                                    r/2./eta*(vrr_terms[i,j,k, q,r-1,0, im]
                                             - zeta/(zeta+eta)
                                             *vrr_terms[i,j,k, q,r-1,0, im+1]
                                             ))
                            if j:
                                vrr_terms[i,j,k, q,r+1,0, im] += (
                                    j/2./(zeta+eta)*vrr_terms[i,j-1,k,q,r,0,im+1]
                                    )

    for s in xrange(nc):
        for r in xrange(mc+1):
            for q in xrange(lc+1):
                for k in xrange(na+1):
                    for j in xrange(ma+1):
                        for i in xrange(la+1):
                            for im in xrange(mtot-i-j-k-q-r-s):
                                vrr_terms[i,j,k,q,r,s+1,im] = (
                                    (qz-zc)*vrr_terms[i,j,k,q,r,s,im]
                                    + (wz-qz)*vrr_terms[i,j,k,q,r,s,im+1]
                                    )
                                if s:
                                    vrr_terms[i,j,k,q,r,s+1,im] += (
                                        s/2./eta*(vrr_terms[i,j,k,q,r,s-1,im]
                                                 - zeta/(zeta+eta)
                                                 *vrr_terms[i,j,k,q,r,s-1,im+1]
                                                 ))
                                if k:
                                    vrr_terms[i,j,k,q,r,s+1,im] += (
                                        k/2./(zeta+eta)*vrr_terms[i,j,k-1,q,r,s,im+1]
                                        )
    return vrr_terms[la,ma,na,lc,mc,nc,M]

# Implement the interface to coulomb_repulsion and contr_coulomb
coulomb_repulsion = hrr

def contr_coulomb(aexps,acoefs,anorms,xyza,powa,
                  bexps,bcoefs,bnorms,xyzb,powb,
                  cexps,ccoefs,cnorms,xyzc,powc,
                  dexps,dcoefs,dnorms,xyzd,powd):
    return contr_hrr(xyza,anorms,powa,aexps,acoefs,
                     xyzb,bnorms,powb,bexps,bcoefs,
                     xyzc,cnorms,powc,cexps,ccoefs,
                     xyzd,dnorms,powd,dexps,dcoefs)

def test_contr():
    from basis_sto3g import basis_data
    from Molecule import Molecule
    from Ints import getbasis
    from time import time
    
    r = 1/0.52918
    atoms=Molecule('h2o',atomlist = [(8,(0,0,0)),(1,(r,0,0)),(1,(0,0,r))])
    bfs = getbasis(atoms,basis_data)
    
    o_1s = bfs[0]
    o_2s = bfs[1]
    o_px = bfs[2]
    o_py = bfs[3]
    o_pz = bfs[4]
    h1_s = bfs[5]
    h2_s = bfs[6]
    

    t0 = time()
    val = \
        contr_coulomb(h2_s.pexps,h2_s.pcoefs,h2_s.pnorms,
                      h2_s.origin,h2_s.powers,
                      h2_s.pexps,h2_s.pcoefs,h2_s.pnorms,
                      h2_s.origin,h2_s.powers,
                      h2_s.pexps,h2_s.pcoefs,h2_s.pnorms,
                      h2_s.origin,h2_s.powers,
                      o_pz.pexps,o_pz.pcoefs,o_pz.pnorms,
                      o_pz.origin,o_pz.powers)
    t1 = time()
    print val,t1-t0

def test_vrr():
    import time
    xa,ya,za = 0.,0.,0.
    xb,yb,zb = 0.,0.,0.
    xc,yc,zc = 0.,0.,0.
    xd,yd,zd = 0.,0.,0.
    norma = normb = normc = normd = 1.
    alphaa = alphab = alphac = alphad = 1.

    la,ma,na = 0,0,0
    lc,mc,nc = 0,0,0

    M = 0
    t0 = time.time()
    val1 = vrr((xa,ya,za),norma,(la,ma,na),alphaa,
               (xb,yb,zb),normb,alphab,
               (xc,yc,zc),normc,(lc,mc,nc),alphac,
               (xd,yd,zd),normd,alphad,M)
    t1 = time.time()
    val2 = vrr((xc,yc,zc),normc,(lc,mc,nc),alphac,
               (xd,yd,zd),normd,alphad,
               (xa,ya,za),norma,(la,ma,na),alphaa,
               (xb,yb,zb),normb,alphab,M)
    t2 = time.time()
    print "Values:  ",val1,val2
    print "Timings: ",t1-t0,t2-t1
    return

test = test_contr

if __name__ == '__main__':
    doprofile = 0
    if doprofile:
        import profile,pstats
        profile.run('test()','hgpprof.dat')
        profdat = pstats.Stats('hgpprof.dat')
        profdat.strip_dirs().sort_stats('time').print_stats(8)
    else:
        test()
    
