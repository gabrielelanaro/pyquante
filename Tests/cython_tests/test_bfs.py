from PyQuante.PGBF import PGBF
from PyQuante.CGBF import CGBF
import unittest

class TestCGBF(unittest.TestCase):
    """
    """
    def setUp(self):
        self.gto = CGBF((0,0,0),(1,0,0))
    def testAddGTO(self):
        pgto = PGBF(0.4,(0,1,1),(0,0,0))
        self.gto.add_primitive(0.4,0.5)

class TestPGBF(unittest.TestCase):
    def setUp(self):
        self.gto = PGBF(0.4,(0,1,1),(0,0,0))

if __name__ == '__main__':
    unittest.main()
        
