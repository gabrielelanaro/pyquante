/**
 * @file primitive-gto.h
 * @author Gabriele Lanaro (gabriele.lanaro@gmail.com)
 *
 * This source code is part of the PyQuante Quantum Chemistry suite.
 * 
 * This file is very important, defines all the functions related to
 * the primitive gaussian basis functions.
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

#ifndef _PRIMITIVE_GTO_H_
#define _PRIMITIVE_GTO_H_

/**
 * Struct representing the ContractedGTO
 */
typedef struct {
    double	alpha;		/**< exponent */
    double	x0,y0,z0;	/**< center coordinates */
    int		l,m,n;		/**< angular momentum values */
    double	norm;		/**< normalization factor */
    double	coef;		/**< coefficient of the expansion, if
				     in a ContracedGTO */
} PrimitiveGTO;

/**
 * Allocates the GTO and returns the initializated structure
 */
PrimitiveGTO *primitive_gto_new( void );

/**
 * Constructor function:
 *
 * Initializate the pgto with the parameters
 * 
 * @param gto a pointer for the new initializated structure
 */
void primitive_gto_init(PrimitiveGTO *gto, double alpha,
			double x0, double y0, double z0,
			int l, int m, int n,
			double coef);

/**
 * @return a copy of the PrimitiveGTO
 */
PrimitiveGTO *primitive_gto_copy(PrimitiveGTO *pgto);

/**
 * Deallocate the gto
 */
void primitive_gto_free(PrimitiveGTO *pgto);

/**
 * Normalize the pgto
 *
 * At least updates the PrimitiveGTO.norm normalization factor.
 */
void primitive_gto_normalize(PrimitiveGTO *pgto);

/**
 * Changes the center of the gto, utility function.
 */
void primitive_gto_recenter(PrimitiveGTO *pgto, double x,double y,double z);

/**
 * Compute the overlap integral with another gto.
 */
double primitive_gto_overlap(PrimitiveGTO *pgto1, PrimitiveGTO *pgto2);

/**
 * Compute the value of the gto in the point x,y,z
 */
double primitive_gto_amp(PrimitiveGTO *pgto, double x, double y, double z);

/**
 * Compute the electron repulsion integral (12|34) using Rys
 * polynomials method.
 */
/*
double primitive_gto_coulomb(PrimitiveGTO *pgto1, PrimitiveGTO *pgto2, PrimitiveGTO *pgto3, PrimitiveGTO *pgto4);
*/

/**
 * Compute the angular momentum of the pgto defined as pgto->l +
 * pgto->m + pgto->n
 */
int primitive_gto_angular_momentum(PrimitiveGTO *pgto);

/**
 * Utility function, copies pgto->x0,y0,z0 in output array
 */
void primitive_gto_R(PrimitiveGTO *pgto, double output[3]);

/**
 * Utility function to set l,m,n at a time
 */
void primitive_gto_set_powers(PrimitiveGTO * pgto, int l, int m, int n);

/**********************************************************************
 * Internal stuff pending for removal
 */

double overlap_1D(int l1, int l2, double PAx,
		  double PBx, double gamma);

 double product_center_1D(double alphaa, double xa,  
			  double alphab, double xb);

double overlap(double alpha1, int l1, int m1, int n1,
	       double xa, double ya, double za,
	       double alpha2, int l2, int m2, int n2,
	       double xb, double yb, double zb);

#endif /* _PRIMITIVE_GTO_H_ */
