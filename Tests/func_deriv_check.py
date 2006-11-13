#!/usr/bin/env python
"""\
Test the functional derivatives 
"""

from PyQuante.DFunctionals import xpbe,xs,xb,cvwn,clyp

rho = 1.0
gam = 0
d = 1e-5

Funcs = dict(S=xs,B=xb,PBE=xpbe,VWN=cvwn,LYP=clyp)
hasg = dict(S=False,B=True,PBE=True,VWN=False,LYP=True)

for name in ['S','B','PBE']:
    f = Funcs[name]
    print name
    if hasg[name]:
        e,v = f(rho,0)
        e2,v2 = f(rho+d,0)
    else:
        e,v = f(rho)
        e2,v2 = f(rho+d)
    print e,v,(e2-e)/d

for name in ['VWN','LYP']:
    f = Funcs[name]
    print name
    if hasg[name]:
        e,va,vb = f(rho,rho,0,0,0)
        e2,va2,vb2 = f(rho+d,rho+d,0,0,0)
    else:
        e,va,vb = f(rho,rho)
        e2,va2,vb2 = f(rho+d,rho+d)
    print e,va+vb,(e2-e)/d



