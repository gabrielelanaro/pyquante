import sys
import unittest

testmodules = ['be_oep', 'h2', 'he', 'h2o', 'h2o_mindo', 
               'oh_mindo']
dfttests = ['h2_ft_dft', 'he_dft', 'h2o_dft']

def importname(modulename, name):
    """Import from a module whose name is determined at runtime.

    (Python Cookbook 2nd ed.)
    """
    try:
        module = __import__(modulename, globals(), locals(), [name])
    except ImportError:
        return None
    return getattr(module, name)

if __name__=="__main__":
    fullsuite = unittest.TestSuite()
    for testmodule in testmodules:
        fullsuite.addTest(importname(testmodule, "suite")())
    unittest.TextTestRunner(verbosity=2).run(fullsuite)
        
        
