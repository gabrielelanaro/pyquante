from PyQuante.Molecule import Molecule
# Construct some molecules (at HF/6-31G** minimum energy geometries)
#  for use elsewhere
h = Molecule('H',[(1,  (0.00000000,     0.00000000,     0.00000000))],
             multiplicity=2)

h2 = Molecule('H2',
              [(1,  (0.00000000,     0.00000000,     0.36628549)),
               (1,  (0.00000000,     0.00000000,    -0.36628549))],
              units='Angstrom')

h2o = Molecule('H2O',
               [(8,  ( 0.00000000,     0.00000000,     0.04851804)),
                (1,  ( 0.75300223,     0.00000000,    -0.51923377)),
                (1,  (-0.75300223,     0.00000000,    -0.51923377))],
               units='Angstrom')

oh = Molecule('OH',
              [(8,  (0.00000000,     0.00000000,    -0.08687037)),
               (1,  (0.00000000,     0.00000000,     0.86464814))],
               units='Angstrom',
               multiplicity=2)
he = Molecule('He',atomlist = [(2,(0,0,0))])
li = Molecule('Li',atomlist = [(3,(0,0,0))], multiplicity=2)
li_p = Molecule('Li+',atomlist = [(3,(0,0,0))],charge=1)
li_m = Molecule('Li-',atomlist = [(3,(0,0,0))],charge=-1)

lih = Molecule('LiH',
               [(3,    (0.00000000,     0.00000000,    -0.53999756)),
                (1,    (0.00000000,     0.00000000,     1.08999756))],
               units='Angstrom')
               
co = Molecule('CO',
              [(6,  (0.00000000,     0.00000000,    -0.63546711)),
               (8,  (0.00000000,     0.00000000,     0.47832425))],
              units='Angstrom')

ch4 = Molecule('CH4',
               [(6,  ( 0.00000000,     0.00000000,     0.00000000)),
                (1,  ( 0.62558332,    -0.62558332,     0.62558332)),
                (1,  (-0.62558332,     0.62558332,     0.62558332)),
                (1,  ( 0.62558332,     0.62558332,    -0.62558332)),
                (1,  (-0.62558332,    -0.62558332,    -0.62558332))],
               units='Angstrom')

c6h6 = Molecule('C6H6',
                [ (6 ( 0.98735329,     0.98735329,     0.00000000)),
                  (6 ( 1.34874967,    -0.36139639,     0.00000000)),
                  (6 ( 0.36139639,    -1.34874967,     0.00000000)),
                  (6 (-0.98735329,    -0.98735329,     0.00000000)),
                  (6 (-1.34874967,     0.36139639,     0.00000000)),
                  (6 (-0.36139639,     1.34874967,     0.00000000)),
                  (1 ( 1.75551741,     1.75551741,     0.00000000)),
                  (1 ( 2.39808138,    -0.64256397,     0.00000000)),
                  (1 ( 0.64256397,    -2.39808138,     0.00000000)),
                  (1 (-1.75551741,    -1.75551741,     0.00000000)),
                  (1 (-2.39808138,     0.64256397,     0.00000000)),
                  (1 (-0.64256397,     2.39808138,     0.00000000))]
               units='Angstrom')
