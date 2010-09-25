from NumWrap import zeros,dot,matrixmultiply

class MG2:
    """
    MG2 contains an experimental Molecular Grid for PyQuante.

    Data objects:
    =============
    version:
      Integer, the type of the grid. 1 was the original, these are v=2
      
    ng:
      Integer, total number of grid points.

    nbf:
      Integer: number of basis functions

    do_grad_dens:
      Boolean, whether or not we require a density gradient

    xyzw:
      ng x 4 ndarray, the x,y,z, and weigt of each grid point

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
        self.version = 2
        self.do_grad_dens = kwargs.get('do_grad_dens',False)
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
        if self.do_grad_dens:
            self.bfgrads = zeros((self.ng,self.nbf,3),'d')
            for ibf in xrange(self.nbf):
                for ig in xrange(self.ng):
                    x,y,z,w = self.xyzw[ig,:]
                    self.bfgrads[ig,ibf,:] = bfs[ibf].grad(x,y,z)
        return

    def floor_density(self,tol=1e-9):
        """
        Set density values below tol to zero
        """
        self.density[self.density<tol] = 0
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
        nat = len(self.atoms)
        for iat in xrange(nat):
            ati = self.atoms[iat]
            npts = len(atomgrids[iat])
            for i in xrange(npts):
                point = atomgrids[iat].points[i]
                xp,yp,zp,wp = point.xyzw()
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
        """Create the xyzw array from the atomgrids, which
        can then be discarded."""
        self.xyzw = zeros((self.ng,4),'d')
        #self.iatom = zeros(self.ng,'l') # Can also keep the atom index, if desired
        ig = 0
        nat = len(atomgrids)
        for iat in xrange(nat):
            npts = len(atomgrids[iat])
            for i in xrange(npts):
                point = atomgrids[iat].points[i]
                self.xyzw[ig,:] = point.xyzw()
                #self.iatom[ig] = iat
                ig += 1
        assert ig == self.ng
        return

    def grad_bf_prod(self,a,b):
        "Form grad(chia,chib)."
        gab = zeros((self.ng,3),'d')
        for i in xrange(3):
            gab[:,i] = self.bfgrid[:,a]*self.bfgrads[:,b,i] \
                     + self.bfgrid[:,b]*self.bfgrads[:,a,i]
        return gab

    def set_density(self,D,Db=None):
        """Given either one density matrix, corresponding to
        a spin unpolarized case, or two density matrices, corresponding
        to a spin polarized case, create the density array and,
        if necessary, the gradients"""

        self.density[:,0] = bdb(self.bfgrid,D)
        if Db is None: # Spin unpolarized case
            self.density[:,1] = self.density[:,0]
        else:
            self.density[:,1] = bdb(self.bfgrid,Db)

        if self.do_grad_dens:
            self.grada = 2*bdg(self.bfgrid,D,self.bfgrads)
            # Note: this code was:
            #self.grada = bdg(self.bfgrid,D,self.bfgrads) +\
            #             gdb(self.bfgrads,D,self.bfgrid)
            # but I can't see that the gdb part does anything different
            # than the bdg
            self.gamma[:,0] = abdot(self.grada,self.grada)
            if Db is None:
                self.gradb = self.grada
                self.gamma[:,1] = self.gamma[:,0]
                self.gamma[:,2] = self.gamma[:,0]
            else:
                self.gradb = 2*bdg(self.bfgrid,Db,self.bfgrads)
                self.gamma[:,1] = abdot(self.gradb,self.gradb)
                self.gamma[:,2] = abdot(self.grada,self.gradb)
        return

    def renormalize(self,nel):
        factor = nel/dot(self.xyzw[:,3],self.density.sum(1))
        if abs(factor-1) > 1e-2:
            print "Warning: large renormalization factor found in grid renormalization"
            print factor
        # Don't know whether it makes more sense to rescale the weights or the dens.
        # The density seems to make more sense, I guess
        #self.scale_weights(factor)
        self.scale_density(factor)
        return

    def scale_density(self,factor): self.density *= factor
    def scale_weights(self,factor): self.xyzw[:,3] *= factor

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
    def get_nbf(self): return self.nbf
    # This can only go so far, since the internal structure is now
    # quite different. But this might avoid a few crashes during
    # porting to the new grids.

# Need to find a faster way to do these, perhaps using tensordot?
def bdb(b,d):
    """Basis x Density x Basis matrix multiply."""
    db = dot(b,d)
    return abdot(db,b)

def bdg(b,d,g):
    """Basis x Density x Gradient matrix multiply."""
    n,m = b.shape
    _bdg = zeros((n,3),'d')
    db = dot(b,d)
    for j in xrange(3):
        _bdg[:,j] = abdot(db,g[:,:,j])
    return _bdg

def abdot(A,B):
    """
    Multiply two n x m matrices together so that the result is a n-length vector
    (i.e. the part over m is accumulated).
    """
    return (A*B).sum(1)

def new_grid_tester():
    from PyQuante.TestMolecules import he,h2
    from PyQuante.MolecularGrid import MolecularGrid
    from PyQuante.Ints import getbasis
    from PyQuante import SCF
    mol = h2
    gr = MolecularGrid(mol,do_grad_dens=True)
    gr2 = MG2(mol,do_grad_dens=True)
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
    print "test_gamma: ",test_gamma(gr,gr2)

def test_density(old,new):
    d = old.dens()-new.density[:,0]-new.density[:,1]
    return sum(d) < 1e-5

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

def test_gamma(old,new):
    d = old.get_gamma()-2*(new.gamma[:,0]+new.gamma[:,1])
    return sum(d) < 1e-5
    

if __name__ == '__main__':
    new_grid_tester()

