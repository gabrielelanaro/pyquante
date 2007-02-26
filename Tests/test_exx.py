#!/usr/bin/env python
"Test the EXX routines"

from PyQuante.OEP import exx
from PyQuante.Molecule import Molecule
from PyQuante.Ints import getbasis,getints
from PyQuante.hartree_fock import rhf
from PyQuante.dft import dft
import logging
import time

def test_exx():
    logging.basicConfig(filename='test_exx.log',
                        level=logging.INFO,
                        format="%(message)s",
                        filemode='w')
    logging.info("Testing EXX functions")
    logging.info(time.asctime())

    h2o = Molecule('H2O',
                   atomlist = [(8,(0,0,0)),    # the geo corresponds to the 
                               (1,(0.959,0,0)),# one that Yang used
                               (1,(-.230,0.930,0))],
                   units = 'Angstrom')

    h2 = Molecule('H2',
                   atomlist = [(1,(0.,0.,0.7)),(1,(0.,0.,-0.7))],
                   units = 'Bohr')
    he = Molecule('He',atomlist = [(2,(0,0,0))])
    ne = Molecule('Ne',atomlist = [(10,(0,0,0))])
    ohm = Molecule('OH-',atomlist = [(8,(0.0,0.0,0.0)),
                                     (1,(0.971,0.0,0.0))],
                   units = 'Angstrom',
                   charge = -1)
    be = Molecule('Be',atomlist = [(4,(0.0,0.0,0.0))],
                  units = 'Angstrom')
    lih = Molecule('LiH',
                   atomlist = [(1,(0,0,1.5)),(3,(0,0,-1.5))],
                   units = 'Bohr')
    #for atoms in [h2,he,be,lih,ohm,ne,h2o]:
    for atoms in [h2,lih]:
        logging.info("--%s--" % atoms.name)

        # Caution: the rydberg molecules (He, Ne) will need a
        #  much much larger basis set than the default returned
        # by getbasis (which is 6-31G**)
        bfs = getbasis(atoms)
        S,h,Ints = getints(bfs,atoms)

        enhf,orbehf,orbshf = rhf(atoms,integrals=(S,h,Ints))
        logging.info("HF  total energy = %f"% enhf)

        enlda,orbelda,orbslda = dft(atoms,
                                   integrals=(S,h,Ints),
                                   bfs = bfs)
        logging.info("LDA total energy = %f" % enlda)

        energy,orbe_exx,orbs_exx = exx(atoms,orbslda,
                                       integrals=(S,h,Ints),
                                       bfs = bfs,
                                       verbose=True,
                                       opt_method="BFGS")
        logging.info("EXX total energy = %f" % energy)
    return


if __name__ == '__main__': test_exx()
