=====================================
Moving from basis functions to shells
=====================================

PyQuante_ currently generates a list of integrals from a list of basis functions, with no knowledge of where those basis functions came from. Some improvements in efficiency can be made by calling an integral evaluation routine with an entire shell rather than with individual basis functions. For example, the libint_ program computes integrals in this manner. This document contains some notes on how PyQuante can make use of integral generation routines of this type.

How PyQuante currently generates basis sets and integrals
=========================================================
A basis set in PyQuante is currently stored as a Python dictionary. For example, here's the 6-31G** basis set for C::

    [('S',
      [(3047.5248999999999, 0.0018347000000000001),
       (457.36950999999999, 0.014037300000000001),
       (103.94869, 0.068842600000000004),
       (29.210155, 0.23218440000000001),
       (9.2866630000000008, 0.4679413),
       (3.1639270000000002, 0.36231200000000002)]),
     ('S',
      [(7.8682724000000004, -0.11933240000000001),
       (1.8812884999999999, -0.1608542),
       (0.54424930000000005, 1.1434564)]),
     ('P',
      [(7.8682724000000004, 0.068999099999999994),
       (1.8812884999999999, 0.31642399999999998),
       (0.54424930000000005, 0.74430830000000003)]),
     ('S', [(0.16871439999999999, 1.0)]),
     ('P', [(0.16871439999999999, 1.0)])],

This contains the primitive exponents and contraction coefficients for each basis shell. **Here we will define a shell as a set of exponents/coefficients with the same angular momentum**, though other definitions can be used as well. A basis set is constructed via the Ints.getbasis() function. Each time a C atom is added by that function, it adds the above functions. When it adds the s-functions, nothing additional needs to be done, but when it adds the p-functions, it adds the px, py, and pz functions separately. The Ints.sym2powerlist function is what does this::

    sym2powerlist = {
      'S' : [(0,0,0)],
      'P' : [(1,0,0),(0,1,0),(0,0,1)],
      'D' : [(2,0,0),(0,2,0),(0,0,2),(1,1,0),(0,1,1),(1,0,1)],
      'F' : [(3,0,0),(2,1,0),(2,0,1),(1,2,0),(1,1,1),(1,0,2),
             (0,3,0),(0,2,1),(0,1,2), (0,0,3)]
    }

After the basis functions are expanded in this way, no further information on what shell they came from is preserved. The basis exists as only a list of CGBF objects from this point onward. Given this list of basis functions, PyQuante then constructs one- and two-electron integrals. The one-electron integrals are computed as NxN matrices for the overlap (S), kinetic (T), and nuclear attraction potential (V) integrals. The two-electron integrals are compressed so that only the N^4/8 unique integrals are computed. This is achieved via the ijkl2intindex indexing function::

	def ijkl2intindex(i,j,k,l):
	   "Indexing into the get2ints long array"
	   if i<j: i,j = j,i
	   if k<l: k,l = l,k
	   ij = i*(i+1)/2+j
	   kl = k*(k+1)/2+l
	   if ij < kl: ij,kl = kl,ij
	   return ij*(ij+1)/2+kl
	
Thus, rather than directly access the (i,j,k,l) integral, you access the Ints[ijkl2intindex(i,j,k,l)], which takes into account the fact that (i,j,k,l)=(j,i,k,l)=(j,i,l,k), etc. The routine that computes the two-electron integrals is the Ints.get2ints() function::

	def get2ints(bfs):
	    """Store integrals in a long array in the form (ij|kl) (chemists
	    notation. We only need i>=j, k>=l, and ij <= kl"""
	    from array import array
	    nbf = len(bfs)
	    totlen = nbf*(nbf+1)*(nbf*nbf+nbf+2)/8
	    Ints = array('d',[0]*totlen)
	    for i in xrange(nbf):
	        for j in xrange(i+1):
	            ij = i*(i+1)/2+j
	            for k in xrange(nbf):
	                for l in xrange(k+1):
	                    kl = k*(k+1)/2+l
	                    if ij >= kl:
	                        Ints[intindex(i,j,k,l)] = coulomb(bfs[i],bfs[j],
	                                                          bfs[k],bfs[l])
	    return Ints

Each call to coulomb() returns a single two-electron integral, which is inserted into the long integral list.

How a shell-based integral program works
========================================
A shell-based integral program computes an entire shell of integrals. For example, in the above 6-31G** basis set for C, we have 3 s-shells, and 2 p-shells. Some portion of the math can be reused in computing, say, the (s,px,s,px) integral and the (s,py,s,py) integral, and thus it is generally faster to make a single call to an integral routine for the entire (s,p,s,p) shell. Pseudocode for this might look something like::

    def get2ints_shells(shells):
        nsh = len(shells)
		for i,j,k,l in unique_shells(len(shells)):
			integral_shell = coulomb_shell(shells[i],shells[j],shells[k],shells[l])
			for ijkl,int in enumerate(integral_shell):
				Ints[shell2intindex(ijkl,i,j,k,l)]=int
		return Ints
	def unique_shells(n):
		def pair(i,j): return i*(i+1)//2+j
		for i in xrange(n):
			for j in xrange(i+1):
				ij = pair(i,j)
				for k in xrange(n):
					for l in xrange(k+1)
						kl = pair(k,l)
						if ij > kl: 
							yield i,j,k,l
		return
	def shell2intindex(ijkl,i,j,k,l):
		# Haven't figured this one out yet.

This obviously requires generating a list of shells rather than a list of basis functions. Seems like the thing to do is write a function that given a list of atoms and some basis data, generates the shells. From the shells one can generate the basis set and the one- and two-electron integrals. A shell object could look something like::

	class Shell:
		def __init__(self,location,type,prims):
			self.location = location
			self.type = type
			assert self.type in ['S','P','D','F']
			self.prims = prims
			
		def generate_cgbf(self):
			cgbfs = []
			for power in sym2powers[self.type]:
				cbfg = CGBF(location,powers,#some way to hack atomid)
				cgbfs.append(cgbf)
				for exp,coef in self.prims:
					cgbf.add_prim(exp,coef)
			return cgbfs

.. _PyQuante: http://pyquante.sourceforge.net
.. _libint: http://www.files.chem.vt.edu/chem-dept/valeev/software/libint/libint.html

