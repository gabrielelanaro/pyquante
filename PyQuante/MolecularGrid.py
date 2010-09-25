"""\
 MolecularGrid.py Construct molecular grids from patched atomic
 grids. The technique behind this is based upon:
  A.D. Becke, 'A multicenter numerical integration scheme for
   polyatomic molecules.' J. Chem. Phys 88(4) 1988.

 The atomic grids are constructed from atomic grids that use
 Lebedev grids for the angular part, and Legendre grids for
 the radial parts.

 This program is part of the PyQuante quantum chemistry program suite.

 Copyright (c) 2004, Richard P. Muller. All Rights Reserved. 

 PyQuante version 1.2 and later is covered by the modified BSD
 license. Please see the file LICENSE that is part of this
 distribution. 
"""
from math import sqrt
from AtomicGrid import AtomicGrid, Bragg
from NumWrap import array,reshape,zeros,dot
from PyQuante.cints import dist2

class MolecularGrid:
    "Class to hold grid information from patched atomic grids"
    def __init__(self, atoms, nrad=32, fineness=1,**opts):
        self.version = 1
        self.do_grad_dens = opts.get('do_grad_dens',False)
        self.atoms = atoms
        self.nrad = nrad 
        self.fineness = fineness
        self.make_atom_grids(**opts)
        self.patch_atoms(**opts)
        self.calc_length()
        self._points = self.points()
        self._weights = self.weights()
        return

    def __len__(self):
        return self._length

    def floor_density(self,tol=1e-9):
        for agr in self.atomgrids:
            agr.floor_density(tol)
        return

    def calc_length(self):
        self._length = 0
        for agr in self.atomgrids:
            self._length += len(agr)
        return

    def make_atom_grids(self,**opts):
        self.atomgrids = []
        opts['nrad'] = self.nrad
        opts['fineness'] = self.fineness
        for atom in self.atoms:
            atom.grid = AtomicGrid(atom, **opts)
            self.atomgrids.append(atom.grid)
        return

    def patch_atoms_naive(self,**opts):
        """\
        This was the original PyQuante patching scheme. It simply
        cuts off the grid at the voronai polyhedra. That is, if a
        grid point is closer to another nucleus than it is to its
        parent nucleus, its weight is set to zero.
        """
        nat = len(self.atoms)
        for iat in xrange(nat):
            ati = self.atoms[iat]
            npts = len(self.atomgrids[iat])
            for i in xrange(npts):
                point = self.atomgrids[iat].points[i]
                xp,yp,zp,wp = point.xyzw()
                rip2 = dist2(ati.pos(),(xp,yp,zp))
                for jat in xrange(nat):
                    if jat == iat: continue
                    atj = self.atoms[jat]
                    rjp2 = dist2(atj.pos(),(xp,yp,zp))
                    if rjp2 < rip2: point._w = 0
        return
    
    def patch_atoms(self,**opts):
        """\
        This is Becke's patching algorithm. Attempting to implement
        the normalization that is in eq 22 of that reference.
        """
        nat = len(self.atoms)
        for iat in xrange(nat):
            ati = self.atoms[iat]
            npts = len(self.atomgrids[iat])
            for i in xrange(npts):
                point = self.atomgrids[iat].points[i]
                xp,yp,zp,wp = point.xyzw()
                rip2 = dist2(ati.pos(),(xp,yp,zp))
                rip = sqrt(rip2)
                Pnum = 1
                Pdenom = 0
                for jat in xrange(nat):
                    bap = becke_atomic_grid_p(jat,(xp,yp,zp),self.atoms,**opts)
                    Pdenom += bap
                    if iat == jat: P_iat = bap
                Ptot = P_iat/Pdenom
                point._w *= Ptot
        return

    def points(self):
        "Dynamically form an array of all grid points"
        p = []
        for agr in self.atomgrids: p.extend(agr.points)
        return p

    def set_bf_amps(self,bfs,**opts):
        "Set the basis func amplitude at each grid point"
        for agr in self.atomgrids: agr.set_bf_amps(bfs,**opts)
        self.make_bfgrid()
        if self.do_grad_dens:
            self.make_bfgrads()
        return

    def setdens(self,D,**opts):
        "Set the density at each grid point"
        for agr in self.atomgrids: agr.setdens(D,**opts)
        return

    def weights(self):
        "Return a vector of weights of each point in the grid"
        weights = []
        for agr in self.atomgrids:
            weights.extend(agr.weights())
        return array(weights)
    
    def dens(self):
        "Return the density for each point in the grid"
        ds = []
        for agr in self.atomgrids:
            ds.extend(agr.dens())
        return array(ds)

    def get_gamma(self):
        "Return the density gradient gamma for each point in the grid"
        if not self.do_grad_dens: return None
        gs = []
        for agr in self.atomgrids:
            gs.extend(agr.get_gamma())
        return array(gs)

    def grad(self):
        pts = self._points
        npts = len(pts)
        gr = zeros((npts,3),'d')
        for i in xrange(npts):
            gr[i,:] = pts[i].grad()
        return gr        

    def make_bfgrads(self):
        "Compute gradients over all bfs and all points"
        pts = self._points
        npts = len(pts)
        nbf = len(pts[0].bfgrads[:,0])
        self.bfgrads = zeros((npts,nbf,3),'d')
        for i in xrange(npts):
            self.bfgrads[i,:,:] = pts[i].bfgrads[:,:]
        return

    def grad_bf_prod(self,a,b):
        "Form grad(chia,chib)."
        B = zeros((self._length,3),'d')
        for i in xrange(3):
            B[:,i] = self.bfgrid[:,a]*self.bfgrads[:,b,i] \
                     + self.bfgrid[:,b]*self.bfgrads[:,a,i]
        return B

    def bfs(self,i):
        "Return a basis function over the entire grid"
        bfs = []
        for agr in self.atomgrids:
            bfs.extend(agr.bfs(i))
        return array(bfs)

    def get_nbf(self):
        return self.atomgrids[0].nbf()

    def npts(self):
        npts = 0
        for agr in self.atomgrids: npts += agr.npts()
        return npts

    def renormalize(self,nel):
        "Renormalize the density to sum to the proper number of electrons"
        factor = nel/dot(self.weights(),self.dens())
        if abs(factor-1) > 1e-2:
            print "Warning: large renormalization factor found in grid renormalization"
            print factor
        # Don't know whether it makes more sense to rescale the weights or the dens.
        # The density seems to make more sense, I guess
        #self.scale_weights(factor)
        self.scale_density(factor)
        return

    def scale_weights(self,factor):
        for agr in self.atomgrids:
            agr.scale_weights(factor)
        return

    def scale_density(self,factor):
        for agr in self.atomgrids:
            agr.scale_density(factor)
        return

    def make_bfgrid(self):
        "Construct a matrix with bfs in columns over the entire grid, "
        " so that R[0] is the first basis function, R[1] is the second..."
        bfs = []
        for agr in self.atomgrids:
            bfs.extend(agr.make_bfgrid())
        self.bfgrid = array(bfs)
        npts = self.npts()
        nbf,nrem = divmod(len(self.bfgrid),npts)
        if nrem != 0: raise Exception("Remainder in divmod make_bfgrid")
        nbf2 = self.get_nbf()
        if nbf != nbf2: raise Exception("Wrong # bfns %d %d" % (nbf,nbf2))
        self.bfgrid = reshape(bfs,(npts,nbf))
        return

# These are the functions for the becke projection operator
def fbecke(x,n=3):
    for i in xrange(n): x = pbecke(x)
    return x
def pbecke(x): return 1.5*x-0.5*pow(x,3)
def sbecke(x,n=3): return 0.5*(1-fbecke(x,n))

def becke_atomic_grid_p(iat,(xp,yp,zp),atoms,**opts):
    do_becke_hetero = opts.get('do_becke_hetero',True)
    nat = len(atoms)
    sprod = 1
    ati = atoms[iat]
    rip2 = dist2(ati.pos(),(xp,yp,zp))
    rip = sqrt(rip2)
    for jat in xrange(nat):
        if jat == iat: continue
        atj = atoms[jat]
        rjp2 = dist2(atj.pos(),(xp,yp,zp))
        rjp = sqrt(rjp2)
        rij2 = dist2(ati.pos(),atj.pos())
        rij = sqrt(rij2)
        mu = (rip-rjp)/rij
        # Modify mu based on Becke hetero formulas (App A)
        if do_becke_hetero and ati.atno != atj.atno:
            chi = Bragg[ati.atno]/Bragg[atj.atno]
            u = (chi-1.)/(chi+1.)
            a = u/(u*u-1)
            a = min(a,0.5)
            a = max(a,-0.5)
            mu += a*(1-mu*mu)
        sprod *= sbecke(mu)
    return sprod

if __name__ == '__main__':
    # Test the becke projection grids
    from PyQuante.Molecule import Molecule
    h2 = Molecule('h2',
                  atomlist = [(1,(0.,0.,0.7)),(1,(0.,0.,-0.7))],
                  units = 'Bohr')
    grid = MolecularGrid(h2,do_becke=True)
    
