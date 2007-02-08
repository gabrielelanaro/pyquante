import sys
import unittest

testmodules = ['h2','he','he_dft','h2o','h2o_mindo','oh_mindo','h2o_dft',
               'ne','no_uhf','h2_cis','h2_mp2','h2_dft','lih_dft','h_dft',
               'li_dft','no_dft','h2_ft_dft','li_ft_dft','be_oep',
               'ricks_unit']

def importname(modulename, name):
    """Import from a module whose name is determined at runtime.

    (Python Cookbook 2nd ed.)
    """
    module = __import__(modulename, globals(), locals(), [name])
    if not module:
        raise ImportError
    return getattr(module, name)

def suite():
    fullsuite = unittest.TestSuite()
    for testmodule in testmodules:
        try:
            fullsuite.addTest(importname(testmodule, "suite")())
        except ImportError:
            print "%s failed!" % testmodule
    return fullsuite

if __name__=="__main__":
    try:
        import psyco
        psyco.full()
        print "Using Psyco!"
    except ImportError:
        pass

    unittest.TextTestRunner(verbosity=2).run(suite())