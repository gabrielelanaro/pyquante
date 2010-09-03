from NumWrap import zeros,dot

class MG2:
    """
    MG2 contains an experimental Molecular Grid for PyQuante.

    Data objects:
    =============
    ng:
      Integer, total number of grid points.

    nbf:
      Integer: number of basis functions

    do_grad:
      Boolean, whether or not we require a density gradient

    xyzw:
      ng x 4 ndarray, the x,y,z, and weigt of each grid point

    iatom:
      ng ndarray, the atom index of each grid point

    density:
      ng x 2 ndarray, the spin-up and spin-down density

    grada:
      ng x 3 ndarray, the three components of the grad of dens a

    gradb:
      ng x 3 ndarray, the three componenet of the grad of denb

    gamma:
      ng x 3 ndarray, the gamma matrix, which contains dot(gradi,gradj)
       for i,j in a,b

    bfgrid:
      ng x nbf ndarray, the basis functions evaluated at the grid points

    bfgrads:
      ng x nbf x 3 ndarray, the gradient of each basis function at each
       grid point. This is the largest storage we require, and at some point
       we need to figure out whether we can reduce this by computing, say,
       basis function gradients on the fly.


    Public Functions:
    =================

    """
    def __init__(self,atoms,nrad=32,fineness=1,**kwargs):
        self.do_grad = kwargs.get('do_grad',False)
        self.atoms = atoms
        self.nrad = nrad
        self.fineness = fineness
        self.make_grid(**kwargs)
        self.zero_density()
        return

    def __len__(self): return self.ng
    def __getitem__(self,item): return self.xyzw[item,:]

    def add_basis(self,bfs):
        # Compute the amplitudes of the basis functions over the grid
        self.nbf = len(bfs)
        self.bfgrid = zeros((self.ng,self.nbf),'d')
        for ibf in xrange(self.nbf):
            for ig in xrange(self.ng):
                x,y,z,w = self.xyzw[ig,:]
                self.bfgrid[ig,ibf] = bfs[ibf].amp(x,y,z)
        if self.do_grad:
            self.bfgrads = zeros((self.ng,self.nbf,3),'d')
            for ibf in xrange(self.nbf):
                for ig in xrange(self.ng):
                    x,y,z,w = self.xyzw[ig,:]
                    self.bfgrads[ig,ibf,:] = bfs[ibf].grad(x,y,z)
        return

    def compute_grid_size(self,atomgrids):
        self.ng = 0
        for ag in atomgrids:
            self.ng += len(ag)
        self._length = self.ng # backwards compatibility
        return

    def make_atomgrids(self,**kwargs):
        from PyQuante.AtomicGrid import AtomicGrid
        atomgrids = []
        for atom in self.atoms:
            atomgrids.append(AtomicGrid(atom,**kwargs))
        self.patch_atoms(atomgrids)
        return atomgrids

    def make_grid(self,**kwargs):
        atomgrids = self.make_atomgrids(**kwargs)
        self.compute_grid_size(atomgrids)
        self.patch_grids(atomgrids)
        return

    def patch_atoms(self,atomgrids,**opts):
        """This is Becke's patching algorithm. Attempting to implement
        the normalization that is in eq 22 of that reference."""
        from PyQuante.MolecularGrid import becke_atomic_grid_p
        from PyQuante.cints import dist2
        from math import sqrt
        nat = len(self.atoms)
        for iat in xrange(nat):
            ati = self.atoms[iat]
            npts = len(atomgrids[iat])
            for i in xrange(npts):
                point = atomgrids[iat].points[i]
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

    def patch_grids(self,atomgrids):
        """Create the xyzw and iatom arrays from the atomgrids, which
        can then be discarded."""
        self.xyzw = zeros((self.ng,4),'d')
        self.iatom = zeros(self.ng,'l')
        ig = 0
        nat = len(atomgrids)
        for iat in xrange(nat):
            npts = len(atomgrids[iat])
            for i in xrange(npts):
                point = atomgrids[iat].points[i]
                self.xyzw[ig,:] = point.xyzw()
                self.iatom[ig] = iat
                ig += 1
        assert ig == self.ng
        return

    def set_density(self,D,Db=None):
        """Given either one density matrix, corresponding to
        a spin unpolarized case, or two density matrices, corresponding
        to a spin polarized case, create the density array and,
        if necessary, the gradients"""

        self.density[:,0] = dot(self.bfgrid,dot(D,self.bfgrid))
        if Db is None: # Spin unpolarized case
            self.density[:,1] = self.density[:,0]
        else:
            self.density[:,1] = dot(self.bfgrid,dot(Db,self.bfgrid))

        if self.do_grad:
            self.grada = dot(self.bfgrid.T,dot(D,self.bfgrads)) +\
                         dot(self.bfgrads.T,dot(D,self.bfgrid))
            self.gamma[:,0] = dot(self.grada.T,self.grada)
            if Db is None:
                self.gradb = self.grada
                self.gamma[:,1] = self.gamma[:,0]
                self.gamma[:,2] = self.gamma[:,0]
            else:
                self.gradb = dot(self.bfgrid.T,dot(D,self.bfgrads)) +\
                             dot(self.bfgrads.T,dot(D,self.bfgrid))
                self.gamma[:,1] = dot(self.gradb.T,self.gradb)
                self.gamma[:,2] = dot(self.grada.T,self.gradb)
        return

    def zero_density(self):
        """Initialize the density matrices for later. If nothing else,
        this assures that I don't have to worry about an undefined
        gamma matrix later on."""
        self.density = zeros((self.ng,2),'d')
        self.grada = zeros((self.ng,3),'d')
        self.gradb = zeros((self.ng,3),'d')
        self.gamma = zeros((self.ng,3),'d')

    # These are some convenience functions to allow MG2 grids to be
    # used in the same way as the old MolecularGrid objects
    def set_bf_amps(self,bfs): self.add_basis(bfs)
    def setdens(self,D): self.set_density(D)
    def dens(self): return self.density[:,0]+self.density[:,1]
    def points(self): return self.xyzw
    def weights(self): return self.xyzw[:,3]
    def grad(self): return self.grada + self.gradb
    def get_gamma(self): return 2*(self.gamma[:,0]+self.gamma[:,1])
    # This can only go so far, since the internal structure is now
    # quite different. But this might avoid a few crashes during
    # porting to the new grids.

def new_grid_tester():
    from PyQuante.TestMolecules import he
    from PyQuante.MolecularGrid import MolecularGrid
    from PyQuante.Ints import getbasis
    from PyQuante import SCF
    mol = he
    gr = MolecularGrid(mol)
    gr2 = MG2(mol)
    print "test_length: ",test_length(gr,gr2)
    print "test_distance: ",test_distance(gr,gr2)

    bfs = getbasis(mol)
    gr.set_bf_amps(bfs)
    gr2.add_basis(bfs)
    print "test_bfgrid: ",test_bfgrid(gr,gr2)

    # This is a little weird, but now use the hf density matrix to
    #  test whether the densities are the same
    hf = SCF(mol)
    hf.iterate()
    gr.setdens(hf.dmat)
    gr2.set_density(hf.dmat)
    print "test_density: ",test_density(gr,gr2)
    # also test gamma

def test_density(old,new):
    d = old.dens()-new.density[:,0]-new.density[:,1]
    return sum(sum(d)) < 1e-5

def test_bfgrid(old,new):
    d = old.bfgrid-new.bfgrid
    return sum(sum(d)) < 1e-5

def test_length(old,new):
    return len(old) == new.ng

def test_distance(old,new):
    from PyQuante.cints import dist2
    points = old.points()
    s = 0
    for i in xrange(new.ng):
        x1,y1,z1,w1 = points[i].xyzw()
        x2,y2,z2,w2 = new[i]
        s += dist2((x1,y1,z1),(x2,y2,z2))
    return s<1e-5
    

if __name__ == '__main__':
    new_grid_tester()

