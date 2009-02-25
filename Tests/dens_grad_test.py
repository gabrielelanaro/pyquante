#!/usr/bin/env python
"Test the routines that compute the density gradient"
from PyQuante.Molecule import Molecule
from PyQuante.dft import dft
from PyQuante.Ints import getbasis,getints
from PyQuante.LA2 import mkdens_spinavg
from math import sqrt

def rand_xyz(**opts):
    "return 3 random floats in [-a,a]"
    from random import random
    a = opts.get('a',2)
    x = 2*a*random()-a
    y = 2*a*random()-a
    z = 2*a*random()-a
    return x,y,z

def get_rho(x,y,z,D,bfs):
    nbf = len(bfs)
    n,m = D.shape
    assert nbf == n
    assert nbf == m
    rho = 0
    for i in range(nbf):
        iamp = bfs[i].amp(x,y,z)
        for j in range(nbf):
            jamp = bfs[j].amp(x,y,z)
            rho += 2*D[i,j]*iamp*jamp
    return rho

def get_grad_rho(x,y,z,D,bfs):
    nbf = len(bfs)
    n,m = D.shape
    assert nbf == n
    assert nbf == m
    gx = gy = gz = 0
    for i in range(nbf):
        iamp = bfs[i].amp(x,y,z)
        igx,igy,igz = bfs[i].grad(x,y,z)
        for j in range(nbf):
            jamp = bfs[j].amp(x,y,z)
            jgx,jgy,jgz = bfs[j].grad(x,y,z)
            gx += 2*D[i,j]*(iamp*jgx+jamp*igx)
            gy += 2*D[i,j]*(iamp*jgy+jamp*igy)
            gz += 2*D[i,j]*(iamp*jgz+jamp*igz)
    return gx,gy,gz

def test_grad(atoms,**opts):
    basis = opts.get('basis',None)
    verbose = opts.get('verbose',False)

    bfs = getbasis(atoms,basis)
    S,h,Ints = getints(bfs,atoms)

    nel = atoms.get_nel()
    enuke = atoms.get_enuke()

    energy, orbe, orbs = dft(atoms,return_flag=1)

    nclosed,nopen = atoms.get_closedopen()
    D = mkdens_spinavg(orbs,nclosed,nopen)

    # Now set up a minigrid on which to evaluate the density and gradients:
    npts = opts.get('npts',1)
    d = opts.get('d',1e-4)
    # generate any number of random points, and compute
    # analytical and numeric gradients:
    print "Computing Grad Rho for Molecule %s" % atoms.name
    maxerr = -1
    for i in range(npts):
        x,y,z = rand_xyz()
        rho = get_rho(x,y,z,D,bfs)
        rho_px = get_rho(x+d,y,z,D,bfs)
        rho_mx = get_rho(x-d,y,z,D,bfs)
        rho_py = get_rho(x,y+d,z,D,bfs)
        rho_my = get_rho(x,y-d,z,D,bfs)
        rho_pz = get_rho(x,y,z+d,D,bfs)
        rho_mz = get_rho(x,y,z-d,D,bfs)
        gx,gy,gz = get_grad_rho(x,y,z,D,bfs)
        gx_num = (rho_px-rho_mx)/2/d
        gy_num = (rho_py-rho_my)/2/d
        gz_num = (rho_pz-rho_mz)/2/d

        dx,dy,dz = gx-gx_num,gy-gy_num,gz-gz_num
        error = sqrt(dx*dx+dy*dy+dz*dz)
        maxerr = max(error,maxerr)
        print " Point  %10.6f %10.6f %10.6f %10.6f" % (x,y,z,error)
        print "  Numerical %10.6f %10.6f %10.6f" % (gx_num,gy_num,gz_num)
        print "  Analytic  %10.6f %10.6f %10.6f" % (gx,gy,gz)
    print "The maximum error in the gradient calculation is ",maxerr
    return
    
def main():
    he = Molecule('He',atomlist=[(2,(0.0,0.5,0.0))])
    test_grad(he,npts=5,d=1e-9)
    r = 0.70
    h2 = Molecule('h2',
                  atomlist = [(1,(1.,r/2.,0)),
                              (1,(1.,-r/2.,0))],
                  units='Angs')
    test_grad(h2,npts=5,d=1e-9)
    return

if __name__ == '__main__': main()

