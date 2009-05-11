from PyQuante import Molecule
from PyQuante.Ints import getbasis

h2 = Molecule('h2',[('H',(0,0,-0.5)),('H',(0,0,0.5))],units='Angs')
basis = getbasis(h2,basis_data='sto-3g')
nbf = len(basis)

print "STO-3G z"
for i in range(nbf):
    for j in range(nbf):
        print i,j,basis[i].multipole(basis[j],0,0,1)

# Sheesh! I don't have 6-31g in PyQuante

print '6-31G CH4 Z'
ch4 = Molecule('ch4',
               [('C',(0.0, 0.0, 0.0)),
                ('H',(0.0, 0.0, 1.083658)),
                ('H',(1.021683, 0.0, -0.361219)),
                ('H',(-0.510841, 0.884804, -0.361219)),
                ('H',(-0.510841,-0.884804, -0.361219)),])

basis = getbasis(ch4,basis_data='6-31g')
nbf = len(basis)

for i in range(nbf):
    for j in range(nbf):
        print i,j,basis[i].multipole(basis[j],0,0,1)


