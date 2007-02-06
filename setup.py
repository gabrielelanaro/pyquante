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
import distutils.sysconfig
sysconfig = distutils.sysconfig.get_config_vars()
import os, sys

# Thanks to Konrad Hinsen for showing me the light about these files!

libs = []
if sys.platform != 'win32':
    return_code = os.system(('%s config/libm_test.c -o' % sysconfig['CC']) +
                            ' config/libm_test >& /dev/null')
    if return_code != 0:
        libs.append('m')

doclines = __doc__.split("\n")

setup(name="PyQuante",
      version="1.6.0",
      description = doclines[0],
      author = "Rick Muller",
      author_email = "rmuller@sandia.gov",
      long_description = "\n".join(doclines[2:]),
      url = "http://pyquante.sourceforge.net",
      licence = "BSD",
      platforms = ["any"],
      classifiers = filter(None,classifiers.split("\n")),
      packages = ['PyQuante'],
      ext_modules=[Extension("PyQuante.cints",["Src/cints.c"],libraries=libs),
                   Extension("PyQuante.chgp",["Src/chgp.c"],libraries=libs),
                   Extension("PyQuante.crys",["Src/crys.c"],libraries=libs)])
