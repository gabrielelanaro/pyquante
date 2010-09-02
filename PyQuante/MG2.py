class MG2:
    """
    MG2 contains an experimental Molecular Grid for PyQuante.

    Data objects:
    =============
    ng:
      Integer, total number of grid points.

    nbf:
      Integer: number of basis functions

    do_spin_polarized:
      Boolean, whether or not the density is spin polarized.

    do_grad_dens:
      Boolean, whether or not we require a density gradient

    xyzw:
      ng x 4 ndarray, the x,y,z, and weigt of each grid point

    iatom:
      ng ndarray, the atom index of each grid point

    density:
      ng x 2 ndarray, the spin-up and spin-down density

    gamma:
      ng x 3 ndarray, the gaa, gbb, and gab gradient terms

    grada:
      ng x 3 ndarray, the three components of the grad of dens a

    gradb:
      ng x 3 ndarray, the three componenet of the grad of denb

    bfgrid:
      ng x nbf ndarray, the basis functions evaluated at the grid points


    Public Functions:
    =================

    """
    def __init__(self,atoms,nrad=32,fineness=1,**kwargs):
        self.do_grad_dens = kwargs.get('do_grad_dens',False)
        self.do_spin_polarized = kwargs.get('do_spin_polarized',False)
        self.atoms = atoms
        self.nrad = nrad
        self.fineness = fineness
        self.make_grid(**kwargs)
        return

    def make_grid(self,**kwargs):
        from AtomicGrid import AtomicGrid
        atomgrids = self.make_atomgrids(**kwargs)
        self.compute_grid_size()
        self.patch_grids()

    def make_atomgrids(self,**kwargs):
        atomgrids = []
        for atom in self.atoms:
            atomgrids.append(AtomicGrid(atom,**kwargs))
        return atomgrids

    def compute_grid_size(self):
        self.ng = 0
        for i in xrange(nat):
            self.ng += len(atomgrids[i])
        return


    
