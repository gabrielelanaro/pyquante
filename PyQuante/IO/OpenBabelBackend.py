'''
Formats that can handle openbabel:

    Key       Description
 ---------    -----------
 CONTCAR:     VASP format,
 POSCAR:      VASP format,
 acr:         ACR format,
 adfout:      ADF output format,
 alc:         Alchemy format,
 arc:         Accelrys/MSI Biosym/Insight II CAR format,
 bgf:         MSI BGF format,
 box:         Dock 3.5 Box format,
 bs:          Ball and Stick format,
 c3d1:        Chem3D Cartesian 1 format,
 c3d2:        Chem3D Cartesian 2 format,
 caccrt:      Cacao Cartesian format,
 can:         Canonical SMILES format.,
 car:         Accelrys/MSI Biosym/Insight II CAR format,
 ccc:         CCC format,
 cdx:         ChemDraw binary format,
 cdxml:       ChemDraw CDXML format,
 cif:         Crystallographic Information File,
 ck:          ChemKin format,
 cml:         Chemical Markup Language,
 cmlr:        CML Reaction format,
 crk2d:       Chemical Resource Kit diagram(2D),
 crk3d:       Chemical Resource Kit 3D format,
 ct:          ChemDraw Connection Table format,
 cub:         Gaussian cube format,
 cube:        Gaussian cube format,
 dat:         Generic Output file format,
 dmol:        DMol3 coordinates format,
 dx:          OpenDX cube format for APBS,
 ent:         Protein Data Bank format,
 fch:         Gaussian formatted checkpoint file format,
 fchk:        Gaussian formatted checkpoint file format,
 fck:         Gaussian formatted checkpoint file format,
 feat:        Feature format,
 fract:       Free Form Fractional format,
 fs:          FastSearching,
 g03:         Gaussian Output,
 g09:         Gaussian Output,
 g92:         Gaussian Output,
 g94:         Gaussian Output,
 g98:         Gaussian Output,
 gal:         Gaussian Output,
 gam:         GAMESS Output,
 gamess:      GAMESS Output,
 gamin:       GAMESS Input,
 gamout:      GAMESS Output,
 gpr:         Ghemical format,
 gukin:       GAMESS-UK Input,
 gukout:      GAMESS-UK Output,
 gzmat:       Gaussian Z-Matrix Input,
 hin:         HyperChem HIN format,
 inchi:       InChI format,
 inp:         GAMESS Input,
 ins:         ShelX format,
 jout:        Jaguar output format,
 log:         Generic Output file format,
 mcdl:        MCDL format,
 mcif:        Macromolecular Crystallographic Information,
 mdl:         MDL MOL/SDF format,
 ml2:         Sybyl Mol2 format,
 mmcif:       Macromolecular Crystallographic Information,
 mmd:         MacroModel format,
 mmod:        MacroModel format,
 mol:         MDL MOL/SDF format,
 mol2:        Sybyl Mol2 format,
 mold:        Molden format,
 molden:      Molden format,
 moo:         MOPAC Output format,
 mop:         MOPAC Cartesian format,
 mopcrt:      MOPAC Cartesian format,
 mopin:       MOPAC Internal,
 mopout:      MOPAC Output format,
 mpc:         MOPAC Cartesian format,
 mpo:         Molpro output format,
 mpqc:        MPQC output format,
 msi:         Accelrys/MSI Cerius II MSI format,
 nwo:         NWChem output format,
 out:         Generic Output file format,
 outmol:      DMol3 coordinates format,
 output:      Generic Output file format,
 pc:          PubChem format,
 pcm:         PCModel Format,
 pdb:         Protein Data Bank format,
 png:         PNG files with embedded data,
 pqr:         PQR format,
 pqs:         Parallel Quantum Solutions format,
 prep:        Amber Prep format,
 qcout:       Q-Chem output format,
 res:         ShelX format,
 rsmi:        Reaction SMILES format,
 rxn:         MDL RXN format,
 sd:          MDL MOL/SDF format,
 sdf:         MDL MOL/SDF format,
 smi:         SMILES format,
 smiles:      SMILES format,
 sy2:         Sybyl Mol2 format,
 t41:         ADF TAPE41 format,
 tdd:         Thermo format,
 therm:       Thermo format,
 tmol:        TurboMole Coordinate format,
 txt:         Title format,
 unixyz:      UniChem XYZ format,
 vmol:        ViewMol format,
 xml:         General XML format,
 xtc:         XTC format,
 xyz:         XYZ cartesian coordinates format,
 yob:         YASARA.org YOB format 
'''
from Handlers import StringHandler,FileHandler,FormatUnsupported
from Data import Data
import openbabel,pybel

class BabelStringHandler(StringHandler):
    def read(self, string, format):
        try:
            mol = pybel.readstring(format, string)
        except ValueError:
            raise  FormatUnsupported("Can't recognize the format '%s' supplied"%format)
        atoms = [(atom.atomicnum,atom.coords) for atom in mol.atoms]
        data = Data()
        data.build_molecule(atomlist = atoms,
                            multiplicity = mol.spin,
                            charge = mol.charge )
        return data

    def write(self, data, format):
        mol = openbabel.OBMol()
        molecule = data.molecule
        for atom in molecule.atoms:
            a = mol.NewAtom()
            a.SetAtomicNum(atom.atno)   # carbon atom
            a.SetVector(*atom.r) # coordinates
            # Bonds, do we need them?
        mol.SetTotalCharge (molecule.charge)
        mol.SetTotalSpinMultiplicity (molecule.multiplicity)
        pybelmol = pybel.Molecule(mol)
        try:
             return pybelmol.write(format)
        except ValueError:
            raise FormatUnsupported("Can't recognize the format '%s' supplied"%format)


class BabelFileHandler(FileHandler):
    def __init__(self):
        self.string_handler=BabelStringHandler()
    def guess_format(self,filename):
        conv = openbabel.OBConversion()
        format = conv.FormatFromExt(filename)
        if format==None:
            raise FormatUnsupported("Can't recognize the format of this file")
        else:
            return format
