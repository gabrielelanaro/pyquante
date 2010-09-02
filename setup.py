"""PyQuante: Quantum Chemistry in Python

PyQuante is an open-source suite of programs for developing quantum
chemistry methods. The program is written in the Python programming
language, and has many 'rate-determining' modules also written in C
for speed. The resulting code is not nearly as fast as Jaguar,
Gaussian, or GAMESS, but the resulting code is much easier to
understand and modify. The goal of this software is not necessarily to
provide a working quantum chemistry program (although it will
hopefully do that), but rather to provide a well-engineered set of
tools so that scientists can construct their own quantum chemistry
programs without going through the tedium of having to write every
low-level routine. More information, including links to the download
page, is available at http://pyquante.sourceforge.net.
"""

classifiers = """\
Development Status :: 3 - Alpha
Intended Audience :: Developers
Intended Audience :: Science/Research
License :: OSI Approved :: BSD License
Programming Language :: Python
Topic :: Scientific/Engineering
Topic :: Software Development :: Libraries :: Python Modules
Operating System :: Microsoft :: Windows
Operating System :: Unix
Operating System :: MacOS
"""
from distutils.core import setup, Extension
from setuptools import setup, Extension
import distutils.sysconfig
sysconfig = distutils.sysconfig.get_config_vars()
import os, sys, glob, stat

libs = []
if sys.platform != 'win32':
    return_code = os.system(('%s config/libm_test.c -o' % sysconfig['CC']) +
                            ' config/libm_test >& /dev/null')
    if return_code != 0:
        libs.append('m')

for fname in glob.glob('Tests/*.py'):
    os.chmod(fname,stat.S_IRWXU|stat.S_IRGRP|stat.S_IXGRP|stat.S_IROTH|stat.S_IXOTH)

doclines = __doc__.split("\n")

# Custom library for ease of C coding
lib_includes = ["Src/lib"]

lib_utils = ["Src/lib/utils/math.c","Src/lib/utils/swap.c"] # Generic stuff
lib_pgto = ["Src/lib/primitive-gto.c"] + lib_utils  # PrimitiveGTO code
lib_cgto = ["Src/lib/contracted-gto.c"] + lib_pgto  # ContractedGTO code
lib_shell =["Src/lib/shell.c"] + lib_cgto           # Shell code

# Python extension plus their dependencies
cgto_ext = ["Src/PyQuante/contracted_gto.c"] + lib_cgto
pgto_ext = ["Src/PyQuante/primitive_gto.c"] + lib_pgto
shell_ext = ["Src/PyQuante/shell.c"] + lib_shell

ext_modules=[Extension("PyQuante.cints",["Src/cints.c"],libraries=libs),
             Extension("PyQuante.chgp",["Src/chgp.c"],libraries=libs),
             Extension("PyQuante.crys",["Src/crys.c"],libraries=libs),
             Extension("PyQuante.contracted_gto",
                       cgto_ext,
                       include_dirs=lib_includes,
                       ),
             Extension("PyQuante.primitive_gto",
                       pgto_ext,
                       include_dirs=lib_includes,
                       ),
             Extension("PyQuante.shell",
                       shell_ext,
                       include_dirs = lib_includes,)
             ]

# Fetching command line option
if "--enable-libint" in sys.argv:
    enable_libint = True
    sys.argv.remove("--enable-libint")
else:
    enable_libint = False

if enable_libint:
    # Libint Extension compilation stuff
    libint_static = [] # Libint object files
    for part in ["libint","libderiv","libr12"]:
        libint_static += glob.glob("libint-1.1.4/src/lib/%s/tmp/%s/.libs/*.o"%(part,part))
        
    # Error checking
    if libint_static == []:
        raise Exception("Object files not found, have you compiled libint with configure --enable-shared ?")
    
    # Preparing the extension
    libint_includes = ["libint-1.1.4/include/libint"] # Libint include dir
    
    lib_clibint =["Src/lib/clibint.c"] + lib_shell      # Thin libint wrapper (actually one function, shell_compute)
    clibint_ext_src = ["Src/PyQuante/clibint.c"] + lib_clibint # Need to add libint_static when linking
    
    clibint_ext = Extension("PyQuante.clibint",
                            clibint_ext_src,
                            include_dirs = lib_includes + libint_includes ,
                            libraries = ["stdc++"],
                            extra_objects = libint_static)
    ext_modules.append(clibint_ext)



setup(name="PyQuante",
      version="1.6.3",
      description = doclines[0],
      author = "Rick Muller",
      author_email = "rmuller@sandia.gov",
      long_description = "\n".join(doclines[2:]),
      url = "http://pyquante.sourceforge.net",
      licence = "BSD",
      platforms = ["any"],
      classifiers = filter(None,classifiers.split("\n")),
      packages = ['PyQuante','PyQuante.Basis','PyQuante.IO','PyQuante.IO.FormatHandlers'],
      ext_modules = ext_modules
      )
