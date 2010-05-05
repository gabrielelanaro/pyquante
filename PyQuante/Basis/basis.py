from PyQuante.shell import Shell

from PyQuante.Ints import sym2powerlist
from PyQuante.Basis.Tools import get_basis_data

from PyQuante.CGBF import CGBF

class BasisSet(object):
    """
    """
    def __init__(self, atoms, basis_data = None, **opts):
        """
        
        """
        # Option to omit f basis functions from imported basis sets
        omit_f = opts.get('omit_f',False)
        if not basis_data:
            from PyQuante.Basis.p631ss import basis_data
        elif type(basis_data) == type(''):
            # Assume this is a name of a basis set, e.g. '6-31g**'
            #  and import dynamically
            basis_data = get_basis_data(basis_data)
        
        # Creating the function lists, by shell and by arbitrary order
        bfs = []   # Basis list
        shells = []# Shell list
        for atom in atoms:
            bs = basis_data[atom.atno]
            for sym,prims in bs: # Shell Symbol S,P,D,F
                if omit_f and sym == "F": continue
                shell = Shell(sym)
                for power in sym2powerlist[sym]:
                    cgbf = CGBF(atom.pos(), power, atom.atid)
                    
                    #exps,coefs = zip(*prims)
                    #primlist = [PrimitiveGTO(alpha,atom.pos(),power) for alpha in exps]
                    
                    #bf = ContractedGTO(primlist,coefs)
                    #bf.normalize()
                    [cgbf.add_primitive(alpha,coef) for alpha,coef in prims]
                    bfs.append(cgbf) # Normal ordering
                    
                    # Shell ordering
                    bfs_index = len(bfs)-1 # last added
                    shell.append(cgbf, bfs_index) 
                shells.append(shell)
        
        self.bfs = bfs
        self.shells = shells
        #for index,func in enumerate(self.__iter__()):
        #    func.index = index
    def __len__(self):
        return len(self.bfs)
    def __iter__(self):
        return self.bfs
    def __getitem__(self, item):
        return self.bfs[item]
