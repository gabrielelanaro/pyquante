'''
This module contains the settings related to PyQuante, useful to
selecting things as backends. In this way the settings can also be
modifed at runtime.
'''
import sys

try:
    import openbabel
    openbabel_enabled = True
except ImportError:
    openbabel_enabled = False
    print >> sys.stderr, "openbabel not found in path, switching to PyQuante backend"

try:
    import clibint
    libint_enabled = True
except ImportError:
    libint_enabled = False
    print >> sys.stderr, "libint extension not found, switching to normal ERI computation"
