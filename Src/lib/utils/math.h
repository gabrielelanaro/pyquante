/**********************************************************************
 * 
 * utils.h
 * 
 *   This module contains widely used generic math functions and 
 *   constants.                                                  
 *
 * This source code is part of the PyQuante Quantum Chemistry suite.
 *  
 * Written by Gabriele Lanaro, 2009-2010
 * Copyright (c) 2009-2010, Gabriele Lanaro
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * 
 **********************************************************************/

#ifndef MATHUTIL_H
#define MATHUTIL_H
/**********************************************************************
 * Constants
 */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define SMALL 0.00000001




/**********************************************************************
 * lgamma(x)
 *  
 *  portable definition under different compilers
 */


#if defined(_WIN32) && !defined(_MSC_VER) /* Not required for MSVC
					     since the code is
					     included below */
double lgamma(double x);
#endif

/* lgamma not included in ANSI standard and so not available in
   MSVC  */
#if defined(_MSC_VER)
double lgamma(double z) {
    double c[7];
    double x,y ,tmp, ser, v;
    int i;

    if (z<=0) return 0;

    c[0]=2.5066282746310005;
    c[1]=76.18009172947146;
    c[2]=-86.50532032941677;
    c[3]=24.01409824083091;
    c[4]=-1.231739572450155;
    c[5]=0.1208650973866179e-2;
    c[6]=-0.5395239384953e-5;
   
    x   = z;
    y   = x;
    tmp = x+5.5;
    tmp = (x+0.5)*log(tmp)-tmp;
    ser = 1.000000000190015;
    for (i=1; i<7; i++) {
        y   += 1.0;
        ser += c[i]/y;
        }
    v = tmp+log(c[0]*ser/x);
    return v;
    }
#endif

/**********************************************************************
 * MAX
 * 
 *   Macro that computes the max between two numbers
 */
#define MAX(a,b) (a > b ? a : b)


/**********************************************************************
 * Fgamma(m,x)
 * 
 *   .
 */
double Fgamma(double m, double x);

/**********************************************************************
 * boys function
 */
#define IS_SMALL(x) (x-SMALL) < EPS

double Fm(int m, double x);


/**********************************************************************
 * Factorial functions:
 * fact
 * fact2
 */

int fact(int n);
int fact2(int n);

/**********************************************************************
 * dist2
 */

double dist2(double x1, double y1, double z1,
	     double x2, double y2, double z2);

double dist(double x1, double y1, double z1,
	    double x2, double y2, double z2);

double binomial_prefactor(int s, int ia, int ib, double xpa, double xpb);

int binomial(int a, int b);

/**********************************************************************
 * 
 *  Vector utilities
 * 
 **********************************************************************/

/**********************************************************************
 * vec_dist2(veca[3], vecb[3])
 * 
 * Compute the square of the norm of |veca-vecb|^2
 */
double vec_dist2( double a[3], double b[3]);

/**********************************************************************
 * vec_subtract(A,B,output)
 * 
 *  Subtracts the 3dimensional vectors A, B and put the output vector
 *  in output
 */
void vec_subtract(double A[3],double B[3], double output[3]);

#endif /*MATHUTIL_H */


/**********************************************************************
 * 4 stuff processing utilitiels
 */


int max4(int a, int b, int c , int d);
