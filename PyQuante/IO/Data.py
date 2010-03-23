class Data(object):
    '''
    Handles the information contained in various types of chemical files,
    it can be produced with a FormatHandler
    
    Attributes:
    - molecule: Molecule instance
    - orbitals: 
    - molecules: [Molecule, Molecule, ... ]
    
    '''
    def __init__(self, molecule = None, molecules = None, orbitals = None):
        self.molecule = molecule
        self.molecules = molecules
        self.orbitals = orbitals
    def build_molecule(self, *a, **kw):
        from PyQuante import Molecule
        self.molecule = Molecule(*a, **kw)
    def has(self, name):
        '''
        Check if the Data has this kind of information saved.
        '''
        if hasattr(self,name):
            attr = getattr(self,name)
            return attr
        else:
            return False
