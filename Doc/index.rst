PyQuante: Python Quantum Chemistry
==================================

.. toctree::
   :maxdepth: 2

Overview
========
PyQuante_ (`Sourceforge Project Page`_) is an open-source suite of
programs for developing quantum chemistry methods. The program is
written in the Python programming language, but has many
"rate-determining" modules also written in C for speed. The resulting
code, though not as fast as Jaguar_, NWChem_, Gaussian_, or MPQC_, is
much easier to understand and modify. The goal of this software is not
necessarily to provide a working quantum chemistry program (although
it will hopefully do that), but rather to provide a well-engineered
set of tools so that scientists can construct their own quantum
chemistry programs without going through the tedium of having to write
every low-level routine.

.. _PyQuante: http://pyquante.sourceforge.net
.. _Jaguar: http://www.schrodinger.org
.. _NWChem: http://www.nwchem.org
.. _Gaussian: http://www.gaussian.com
.. _MPQC: http://mpqc.net
.. _`Sourceforge Project Page`: http://sourceforge.net/projects/pyquante

Current features
----------------
* Hartree-Fock: Restriced closed-shell HF and unrestricted open-shell
  HF; 
* DFT: LDA (SVWN, Xalpha) and GGA (BLYP) functionals;
* Optimized-effective potential DFT;
* Two electron integrals computed using Huzinaga, Rys, or
  Head-Gordon/Pople techniques; C and Python interfaces to all of
  these programs; 
* MINDO/3 semiempirical energies and forces;
* CI-Singles excited states;
* DIIS convergence acceleration;
* Second-order Moller-Plesset (MP2) perturbation theory.

Getting Started
---------------

See the cookbook_ for short snippets to get started, and also see the
tests subdirectory of the code distribution. Subscription to the
pyquante-users_ mailing list is highly recommended for further
nsupport. Additionally, people interested in day-to-day development
issues in PyQuante_ are urged to subscribe to the pyquante-devel_
mailing list.

.. _pyquante-users: https://lists.sourceforge.net/lists/listinfo/pyquante-users
.. _pyquante-devel: https://lists.sourceforge.net/lists/listinfo/pyquante-devel

License
-------
The software is released under the modified BSD license, which means
that everyone is free to download, use, and modify the code without
charge.

Download and Installation
==========================
The program is available in tarball form from the PyQuante
`Sourceforge Project Page`_.
The SVN archive for the program is also at Sourceforge, and is
recommended for anyone wanting to stay current with the codebase;
information on how to access the SVN archive is available here.

Prerequisites
-------------
You will need to have a recent (2.3 or later) version of Python
installed on your computer; this can be installed by following
instructions on the Python Website, and binary packages exist for most
platforms. Additionally, you need the multidimensional arrays and
linear algebra libraries provided by either Numeric or numpy; both of
these can be downloaded from the Numpy Website.

Given the choice between Numeric and Numpy, please use Numpy. The
library is much newer, and future developments will center on this
library. Although we will make every effort to continue supporting
Numeric, Numpy is what we will use as the default, and should be much
more stable with PyQuante.

Note: If you're installing on a Macintosh, I recommend the Framework
binaries of Python and Numpy from the Pythonmac Packages Directory.

On Ubuntu and other Linux distros, several additional packages may be
required, including python-dev, python-numpy, and possibly
python-matplotlib.

Building the Code
-----------------
Much of the code is written in Python, and thus is platform
independent. The rest of the code now uses the python distutils
procedures for building the C modules::

    % (sudo) python setup.py install

and the code should build and install properly. I've tested this on
Linux, Windows/Cygwin, and Macintosh OS X. Installing the module like
tsetxkbmap -option ctrl:nocaps       # Make Caps Lock a Control key
his will insure that the modules are installed where Python can find
them.

However, the above assumes that you have write priviledges in the
Python site-packages directory, possibly via the _sudo_ command. If
you do not have access to these directories, create an install
directory where you do have access and pass this directory to the
install process::

    % mkdir ~/pyquante_lib
    % python setup.py install --install-lib=~/pyquante_lib
    % export PYTHONPATH=~/pyquante_lib:$PYTHONPATH

Email me (rmuller@sandia.gov) if you need additional help.

Using the Code
==============

Cookbook
--------

Specifying a molecule
.....................

You can specify a molecule using the following code::

    >>> from PyQuante.Molecule import Molecule
    >>> h2 = Molecule('h2',[(1,(0,0,0)),(1,(1.4,0,0))])


Simple HF calculation
.....................

Using the above definition, you can run a simple HF calculation via::

    >>> from PyQuante.hartree_fock import rhf
    >>> en,orbe,orbs = rhf(h2)
    >>> print "HF Energy = ",en

We're also working on a more object-oriented interface (called
*PyQuante 2*) that will hopefully decrease the amount of duplicated
code in the project. You can run the same calculation as above using::

    >>> from PyQuante import SCF
    >>> solver = SCF(h2,method="HF")
    >>> solver.iterate()
    >>> print "HF Result = ",solver.energy

With the new solvers, you can run with an alternate basis set (say,
STO-3G), via::

    >>> solver = SCF(h2,method="HF",basis="sto-3g")
    >>> solver.iterate()
    >>> print "HF Result = ",solver.energy

At this geometry (R=1.4 bohr) this should produce an energy of -1.1313
hartrees. Note that you may also import Molecule and SCF directly
from the PyQuante namespace now.

Simple DFT calculation
......................

We can run a DFT calculation on the same molecule by running
the commands::

    >>> from PyQuante.dft import dft
    >>> en,orbe,orbs = dft(h2)

This will produce an energy of -1.1353 hartrees. Again, the 6-31G**
basis set is used by default. In DFT calculations, the functional
defaults to SVWN (LDA). To use a different functional, you can type::

    >>> en,orbe,orbs = dft(h2,functional='BLYP')

which will produce an energy of -1.1665 hartrees. 

With the new solvers this calculation could be run via::

    >>> from PyQuante import SCF
    >>> lda = SCF(h2,method="DFT")
    >>> lda.iterate()
    >>> blyp = SCF(h2,method="DFT",functional="BLYP")
    >>> blyp.iterate()
    >>> print "DFT Results: LDA = %f  BLYP = %f" % (lda.energy,blyp.energy)

Users Guide
-----------

Here we provide a bit more detail about the PyQuante functions.

Molecules
.........
PyQuante programs use the Molecule object to contain the information
about the molecule - the atoms, the charge, and the multiplicity. The
syntax for Molecule is::

  >>> Molecule(name,atomlist,**opts)

The atomlist is a list of atomic numbers and a tuple with the x,y,z
coordinates. Here's an example for constructing a molecule object for
water::

  >>> h2o=Molecule('h2o',[(8,(0,0,0)),(1,(1.0,0,0)),(1,(0,1.0,0))],units = 'Angstrom')

(of course the bond-angle is 90 degrees here, and is thus completely
wrong, but this is only an example). Here's an example for the
hydroxide ion that shows the use of the charge field:: 

  >>> oh = Molecule('OH-',[(8,(0,0,0)),(1,(0.96,0,0))], units = 'Angstrom', charge=-1)

Here's an example for the NO molecule that shows the use of the multiplicity field::

  >>> no = Molecule('NO', [(7,(0,0,0)),(8,(2.12955,0,0))],multiplicity=2)

As of version 1.5.1, you may construct molecules using the atomic symbol instead of the atomic number, e.g.::

  >>> h2o=Molecule('h2o',[('O',(0,0,0)), ('H',(1.0,0,0)),('H',(0,1.0,0))],units = 'Angstrom')

Currently, the semiempirical code uses an extended verion of the
Molecule object that adds a variety of additional features. Upcoming
releases will hopefully unify the use of the Molecule between the HF,
DFT, and semiempirical codes.

As of version 1.6, you can directly import Molecule from the PyQuante
module::

  >>> from PyQuante import Molecule


Basis Sets
..........

Basis functions are constructed using the CGBF (contracted Gaussian
basis function) object, which, in turn, uses the PGBF (primitive
Gaussian basis function) object. Basis sets are simply lists of
CGBF's. In the Ints module there is a convenience function getbasis
that constructs basis sets for different molecules. The syntax for the
getbasis function is:: 

  >>> bfs = getbasis(molecule,basis_data=None)

The basis data can be input from a number of data files in the
PyQuante suite. Here are some of the more commonly used basis sets: 

* basis_631ss The Pople 6-31G** basis set
* basis_sto3g The Pople STO-3G basis set
* basis_321 The Pople 3-21G basis set
* basis_ccpvtz The Dunning cc-pVTZ basis set
* basis_ccpvtzmf The Dunning cc-pVTZ(-f) basis set (cc-pVTZ without
  f-functions) 

For, example, to construct a basis set for the h2o object for water
created above, we would call:: 

  >>> from basis_631ss import basis_data
  >>> bfs = getbasis(h2o,basis_data)

If the basis_data argument is omitted, the program will default to
6-31G**. As of version 1.6, you may now specify a string for the
basis_data argument, e.g. "6-31G**" so you may do things like::

  >>> bfs = getbasis(h2o,'6-31G**')
  >>> bfs2 = getbasis(h2o,'sto-3g')

and so on. Use PyQuante.Basis.Tools.basis_map.keys() for a list of the
supported basis strings.

Output
..................

If you want to check the output of PyQuante and to have a feedback of
what's happening, you can use the `configure_output` routine before
launching the calculation::

    >>> from PyQuante import configure_output
    >>> configure_output()

This will display to the standar output the log of PyQuante.  You can
also specify other ways to handle the output::

    >>> configure_output("calc.log") # save also in a log file and display on stdout 
    >>> configure_output(stream=sys.stderr) # Redirect the output on a generic stream
    >>> configure_output(filename="h2.log", stream=None) # Suppress stream output

Internally, the output is handled using the `logging` module
(available in the standard library), you can access the logger used by
PyQuante in this way::

    >>> import logging
    >>> logger = logging.getLogger("pyquante")

Further information about using Loggers:
http://docs.python.org/library/logging.html

Integrals
.........
One-electron integrals
++++++++++++++++++++++
The one-electron integrals consist of the overlap matrix S, the
kinetic energy matrix T, and the nuclear attraction matrix V. The
latter two are often combined into the one-electron Hamiltonian h. 

There are a number of helper functions in the Ints module:

* getT: Form the kinetic energy matrix T
* getS: Form the overlap matrix S
* getV: Form the nuclear attraction matrix V
* get1ints: Form and return S,h
* getints: Form and return S,h,Ints, where Ints are the two-electron
  integrals (see below). 

These functions actually call instance functions of the CGBF objects,
which can themselves be used individually. Some of the instance
functions in the CGBF module use functions in the pyints and cints
modules. 

Two-electron integrals
++++++++++++++++++++++

The two-electron integrals consist of the electron-electron Coulomb
repulsion interactions. The easiest way to construct these is to use
the functions in the Ints module 

* get2ints: Form and return Ints, where Ints are the two-electron
  integrals (see below). 
* getints: Form and return S,h,Ints, where Ints are the two-electron
  integrals. 

The Ints object (not to be confused with the Ints module, which is
just a collection of helper functions) consists of a list of the
two-electron integrals. 

There are actually three different methods to computing the
two-electron integrals, and the helper functions in the Ints module
default to one of these functions. 

* Huzinaga's original method for computing integrals over Gaussians.
* Rys quadrature.
* Head-Gordon and Pople's recurrance relations.

There are python versions of these methods implemented in the modules
pyints, rys, and hgp, respectively. For speed, there are also
C-versions of these modules in cints, crys, and chgp. The program
defaults to the Coulomb repulsion routines in crys, since these are
the fastest (although recent improvements to chgp make it not much
slower). 

The coulomb_repulsion function is the same in all six modules::

  >>> integral = coulomb_repulsion((xa,ya,za),norma,(la,ma,na),alphaa, 
			(xb,yb,zb),normb,(lb,mb,nb),alphab, 
			(xc,yc,zc),normc,(lc,mc,nc),alphac, 
			(xd,yd,zd),normd,(ld,md,nd),alphad)

This routine computes the Coulomb repulsion integral between four
basis functions. The terms xi, yi, zi are the origins of the different
basis functions. The terms normi are the normalization constants for
the basis function. The terms li, mi, ni are the exponents of the
Cartesian powers for the basis function. And the terms alphai are the
Gaussian exponents.


Hartree-Fock Calculations
.........................
This section is under construction. In the meantime, see the PyQuante
Cookbook recipe `Simple HF calculation`_ for an example of a HF
calculation using PyQuante.

Density Functional Theory Calculations
......................................
This section is under construction. In the meantime, see the PyQuante
Cookbook recipe `Simple DFT calculation`_ for an example of a DFT
calculation using PyQuante.


Semiempirical Calculations
..........................
This section is under construction.




The material on this page is copyright (c) 2009, Richard P. Muller.
Reuse of the material on this page is permitted under either the Gnu
Free Documentation License or the CC-BY-SA License.

