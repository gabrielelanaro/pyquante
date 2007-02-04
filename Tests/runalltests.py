import sys
import unittest

testmodules = ['h2', 'he', 'be_oep', 'oh_mindo']

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
        
        
