import unittest
import os,sys
import difflib

#DIR = os.path.abspath(os.path.join(os.path.dirname(__file__),"../../"))
#sys.path.insert(0,DIR)
from PyQuante import Molecule
SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR,"test_io_data")

CML = os.path.join(DATA_DIR,"h2.cml")
GAMIN = os.path.join(DATA_DIR,"h2.imp")

INEXISTENT =os.path.join(DATA_DIR,"inexistent_file.cml")
UNSUPPORTED = os.path.join(DATA_DIR,"unsupported.unsupported")

PYQ = Molecule("h2",
               [(1,(0,0,0)),
                (1,(1,0,0))]
               )

PYQ_STR='2\n\nH          0.00000        0.00000        0.00000\nH          1.00000        0.00000        0.00000\n'

PYQ_FILE=os.path.join(DATA_DIR,"h2_pyq.xyz")

class TestReadOpenBabel(unittest.TestCase):
    def testLoadFile(self, ):
        """
        """
        mol = Molecule.from_file(CML,format="cml")
        self.assertEqual(mol[0].atno,1)
        self.assertRaises(IOError,Molecule.from_file,INEXISTENT,format="cml")
    def testGuess(self, ):
        """
        """
        mol = Molecule.from_file(CML)
        self.assertRaises(IOError,Molecule.from_file,INEXISTENT)
    def testFormatUnsupported(self):
        from PyQuante.IO.Handlers import FormatUnsupported
        self.assertRaises(FormatUnsupported,Molecule.from_file,UNSUPPORTED)
        self.assertRaises(FormatUnsupported,Molecule.from_file,UNSUPPORTED,format="unsupported")
    def testMultiOptions(self):
        mol = Molecule.from_file(CML,format="cml")
        self.assertEqual(mol.charge,0)
        self.assertEqual(mol.multiplicity,1)
    def testLoadString(self):
        mol = Molecule.from_string(open(CML).read(),"cml")
        self.assertEqual(mol.charge,0)
        self.assertEqual(mol.multiplicity,1)

class TestWriteOpenBabel(unittest.TestCase):
    def setUp(self):
        if os.path.exists(PYQ_FILE):
            os.remove(PYQ_FILE)
    def testWriteString(self):
        mol = PYQ
        res = mol.as_string() # Default x,y,z?
        ratio = difflib.SequenceMatcher(None, res, PYQ_STR).ratio()
        self.assert_(ratio > 0.6)
    def testsDump(self):
        mol = PYQ
        mol.dump(PYQ_FILE) # Default molecule.name.xyz ?
        mol.dump(GAMIN,"gamin")
        self.assertTrue(os.path.exists(PYQ_FILE))
        self.assertEqual(open(PYQ_FILE).read(),PYQ_STR)
    def testFormatUnsupported(self):
        from PyQuante.IO.Handlers import FormatUnsupported
        mol = PYQ
        self.assertRaises(FormatUnsupported,mol.dump,UNSUPPORTED)
        self.assertRaises(FormatUnsupported,mol.dump,UNSUPPORTED,format="unsupported")

from  PyQuante.IO.PyQuanteBackend import PyQuanteStringHandler as SH

class TestReadPyQuante(unittest.TestCase):
    def testReadString(self):
        sh = SH()
        data = sh.read(PYQ_STR, format = "xyz")
        self.assertEqual(len(data.molecule.atoms),2 )
from PyQuante.IO.Data import Data
class TestWritePyQuante(unittest.TestCase):
    def testWriteString(self):
        sh = SH()
        data = Data(molecule=PYQ)
        sh.write(data, format="cml")
        

if __name__ == '__main__':
    unittest.main()
