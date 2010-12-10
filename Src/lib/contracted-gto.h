/**
 * @file contracted-gto.h
 * @author Gabriele Lanaro <gabriele.lanaro@gmail.com>
 *
 * This source code is part of the PyQuante Quantum Chemistry suite.
 * 
 * Here are defined the contracted gaussian basis function and related
 * functions.
 *
 * @todo the primitive gto is PART of the contracted gto, most
 * functions pretends that the pgtos are yet instantiated, however we
 * can implement the memory management of the pgto directly in the
 * cgto functions.
 * 
 *
 * Written by Gabriele Lanaro, 2009-2010
 * Copyright (c) 2009-2010, Gabriele Lanaro
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 */

#ifndef _CONTRACTED_GTO_H_
#define _CONTRACTED_GTO_H_

#include "primitive-gto.h"

/** 
 * struct representing a Contracted Gaussian Basis Function 
 */
typedef struct {
  PrimitiveGTO	**primitives;	/**< Array of pointers to
				     PrimitiveGTO */
  int		  nprims;	/**< Number of primitives stored in
				     `primitives` */
  int		  l,m,n;	/**< angular momentum numbers, the
				     same of the primitives */
  double	  x0,y0,z0;	/**< center, the same of the
				     primitives */
  int		  atid;		/**< id of the atom in which this
				     function is centered */
  double	  norm;		/**< Normalization constant */
} ContractedGTO; 


/**
 * Allocation function, initialize also void data structures.
 */
ContractedGTO *contracted_gto_new(void);

/**
 * Constructor function
 *
 * @param atid Atom identification number, if
 * you don't need it (for example if you haven't to access the atom on
 * which the gto is centered ), set it to 0.
 */ 
void contracted_gto_init (ContractedGTO *cgto, double x, double y, double z, int l, int m, int n, int atid);

/**
 * Alternate constructor function, internally it copies the primitive gto passed:
 * 
 * @param pgtos Array of pgtos
 * @param coefs Array of coefficients for each gto
 * @param ngtos Length of the supplied arrays
 */
void contracted_gto_from_primitives(ContractedGTO *cgto, PrimitiveGTO ** pgtos, int ngtos);

/**
 * Add a primitive to the cgto to the linear combination with the coef.
 *
 * @note the pgto is copied
 */
void contracted_gto_add_primitive(ContractedGTO *cgto, PrimitiveGTO *pgto, double coef);

/**
 * Deallocate the gto
 */
void contracted_gto_free(ContractedGTO *cgto);

/**
 * Normalize the gto
 */
void contracted_gto_normalize(ContractedGTO *cgto);

/**
 * Compute the overlap between two cgtos
 */
double contracted_gto_overlap(ContractedGTO *cgto1, ContractedGTO *cgto2);

/**
 * Compute the value of the cgto in the point x,y,z
 */
double contracted_gto_amp(ContractedGTO *cgto, double x, double y, double z);

/**
 * Compute the electron repulsion integral (1234).
 */
/*
double contracted_gto_coulomb(ContractedGTO *cgto1, ContractedGTO *cgto2,
			      ContractedGTO *cgto3, ContractedGTO *cgto4);
*/

/**
 * Recenters the CGTO
 */
void contracted_gto_recenter(ContractedGTO *cgto, double x, double y, double z);

/**
 * Utility to set l,m,n at a time
 */
void contracted_gto_set_powers(ContractedGTO *cgto, int l, int m, int n);

/**
 * Renormalization prefactor part required by the libint library. It is
 * defined as:
 *
 * cgto(powers=(a,b,c)).norm / cgto(powers=(a+b+c, 0,0)).norm
 *
 * @todo renaming it to libint_renorm_part?
 */
double contracted_gto_libint_renorm(ContractedGTO *cgto);

/**
 * Complete renormalization prefactor given by the product of the
 * renorm prefactors.
 *
 * @todo renaming it to libint_renorm?
 */ 
double contracted_gto_renorm_prefactor( ContractedGTO *cgto1, ContractedGTO *cgto2,
					ContractedGTO *cgto3, ContractedGTO *cgto4);

/**
 * @return Angular momentum l+m+n
 */
int contracted_gto_angular_momentum(ContractedGTO *cgto);

/**
 * Utility function, copy x0,y0,z0 in the output array.
 */
void contracted_gto_R(ContractedGTO *cgto, double output[3]);

#endif /* _CONTRACTED_GTO_H_ */
