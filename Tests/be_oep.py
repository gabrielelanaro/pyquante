from PyQuante.Molecule import Molecule
from PyQuante.Ints import getbasis,getints
from PyQuante.hartree_fock import rhf
from PyQuante.OEP import oep_hf_an,oep_hf

name = "Be/OEP"
energy = -7.98104444095

def main(**opts):
    do_oep_an = opts.get('do_oep_an',True)
    mol = Molecule('LiH',[(1,(0,0,1.5)),(3,(0,0,-1.5))],units = 'Bohr')

    bfs = getbasis(mol)
    S,h,Ints = getints(bfs,mol)
    E_hf,orbe_hf,orbs_hf = rhf(mol,bfs=bfs,integrals=(S,h,Ints),
                               DoAveraging=True)
    if do_oep_an:
        E_exx,orbe_exx,orbs_exx = oep_hf_an(mol,orbs_hf,bfs=bfs,
                                            integrals=(S,h,Ints))
    else:
        E_exx,orbe_exx,orbs_exx = oep_hf(mol,orbs_hf,bfs=bfs,
                                         integrals=(S,h,Ints))
    return E_exx

def test():
    from PyQuante import logging
    logging.basicConfig(#filename='testsuite.log',
                        level=logging.INFO,
                        format="%(message)s",
                        filemode='w')
    logging.info("Running %-15s:" % name)
    E = main()
    error = abs(E-energy)
    if error < 1.e-4:
        logging.info("--- E=%12.6f Worked ---" % E)
    else:
        logging.info("*** Warning: E=%12.6f should be %12.6f ***" %
                     (E,energy))
    

if __name__ == '__main__': test()
