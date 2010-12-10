/*
 * This source code is part of the PyQuante Quantum Chemistry suite.
 *  
 * Written by Gabriele Lanaro, 2009-2010
 * Copyright (c) 2009-2010, Gabriele Lanaro
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 */
#include <math.h>
#include <stdio.h>
#include <assert.h>

#include "math.h"

/* Constants */


/* Prototypes for internal use */
static double gamm_inc(double a, double x);
static void gser(double *gamser, double a, double x, double *gln);
static void gcf(double *gammcf, double a, double x, double *gln);


double Fgamma(double m, double x){
  double val;
  if (fabs(x) < SMALL) x = SMALL;
  val = gamm_inc(m+0.5,x);
  /* if (val < SMALL) return 0.; */ /* Gives a bug for D orbitals. */
  return 0.5*pow(x,-m-0.5)*val; 
}

static double gamm_inc(double a, double x){ /* Taken from NR routine gammap */
  double gamser,gammcf,gln;
  
  assert (x >= 0.);
  assert (a > 0.);
  if (x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return exp(gln)*gamser;
  } else {
    gcf(&gammcf,a,x,&gln);
    return exp(gln)*(1.0-gammcf);
  }
}
 
static void gser(double *gamser, double a, double x, double *gln){
  int n;
  double sum,del,ap;

  *gln=lgamma(a);
  if (x <= 0.0) {
    assert(x>=0.);
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	*gamser=sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    printf("a too large, ITMAX too small in routine gser");
    return;
  }
}
 
static void gcf(double *gammcf, double a, double x, double *gln){
  int i;
  double an,b,c,d,del,h;
  
  *gln=lgamma(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  assert(i<=ITMAX);
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
}

int max4(int a, int b, int c , int d) {
  int ret;
  ret = 0;
  ret = MAX(ret,a);
  ret = MAX(ret,b);
  ret = MAX(ret,c);
  ret = MAX(ret,d);
  
  return ret;
  }

double vec_dist2(double a[3], double b[3]){
  return pow(a[0]-b[0],2) + pow(a[1]-b[1],2) + pow(a[2]-b[2],2); 
}
void vec_subtract(double A[3],double B[3], double output[3]){
  output[0] = A[0] - B[0];
  output[1] = A[1] - B[1];
  output[2] = A[2] - B[2];
}


int fact(int n){
  if (n <= 1) return 1;
  return n*fact(n-1);
}

int fact2(int n){ /* double factorial function = 1*3*5*...*n */
  if (n <= 1) return 1;
  return n*fact2(n-2);
}


double dist2(double x1, double y1, double z1,
	     double x2, double y2, double z2){
  return (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2);
}

double dist(double x1, double y1, double z1,
	    double x2, double y2, double z2){
  return sqrt(dist2(x1,y1,z1,x2,y2,z2));
}

double binomial_prefactor(int s, int ia, int ib, double xpa, double xpb){
  int t;
  double sum=0.;
  for (t=0; t<s+1; t++)
    if ((s-ia <= t) && (t <= ib)) 
      sum += binomial(ia,s-t)*binomial(ib,t)*pow(xpa,ia-s+t)*pow(xpb,ib-t);
  return sum;
} 

int binomial(int a, int b){ 
  return fact(a)/(fact(b)*fact(a-b)); 
}


double Fm(int m, double x){
  double precedent, current;
  int k;
  
  if ( IS_SMALL(x) ){
    k = 1;

    current = 1/(2*m+1);
    do  {
      precedent = current;
      current += pow(-x,k) / ( fact(k) * (2*m + 2*k +1));
      k++;
    }
    while (fabs(precedent - current) > EPS);
    
    return current;
  }
  else {
    return fact2(2*m -1)/pow(2.0,m+1) * sqrt( M_PI/pow( x, 2*m+1 ));
  }
}
