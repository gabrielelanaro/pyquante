from PyQuante import Molecule
from PyQuante.Ints import getbasis

h2 = Molecule('h2',[('H',(0,0,-0.5)),('H',(0,0,0.5))],units='Angs')
basis = getbasis(h2,basis_data='sto-3g')
nbf = len(basis)

for i in range(nbf):
    for j in range(nbf):
        print i,j,basis[i].multipole(basis[j],0,0,1)


