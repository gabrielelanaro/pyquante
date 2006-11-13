/* This is a template file that is read in by make_hf_driver.py and
 * populated with data. You can write it pretty much like normal C
 * code, except when you see a percent sign. If you're printing out
 * data from the python driver, you need to use the Python "dictionary
 * style" of expanding strings. If you're putting in a print statement
 * to be executed by the C program, you need to double the percent
 * signs, because Python is reading and then rewriting the text.
 *
 * This program is part of PyQuante, and is distributed under a
 * modified BSD license. Copyright (c) 2006, Richard P. Muller.
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "cints.h"
#include "qc.h"
#include "linalg.h"

int main(){
  /* Simple example program for using cints.c routines */

  int i=0;
  int nat = %(nat)d;
  int nbf = %(nbf)d;
  int nocc=%(nclosed)d,info;
  int iter,maxiter=%(maxiter)d;
  double h[%(nbf2)d],S[%(nbf2)d],T[%(nbf2)d],V[%(nbf2)d],Ints[%(nint)d];
  double D[%(nbf2)d],J[%(nbf2)d],K[%(nbf2)d];
  double orbe[%(nbf)d],orbs[%(nbf2)d],temp[%(nbf2)d];
  double temp1[%(nbf)d], temp2[%(nbf)d];
  double H[%(nbf2)d];
  double energy,eone,ej,ek;
  double enuke=%(enuke)f;

  /*  Data over atoms  */
%(atno_array)s
%(x_array)s
%(y_array)s
%(z_array)s

  /*  Data over contracted bfns  */
%(xcenter_array)s  
%(ycenter_array)s  
%(zcenter_array)s  
%(lpower_array)s
%(mpower_array)s
%(npower_array)s
%(normc_array)s
%(istart_array)s
%(nprim_array)s

  /*  Data over primitive bfns  */
%(alpha_array)s
%(coef_array)s
%(normp_array)s

  getS(nbf,nprim,normc,istart,xcenter,ycenter,zcenter, 
       lpower,mpower,npower,normp,coef,alpha,S);
  getT(nbf,nprim,normc,istart,xcenter,ycenter,zcenter, 
       lpower,mpower,npower,normp,coef,alpha,T);
  getV(nbf,nprim,normc,istart,xcenter,ycenter,zcenter, 
       lpower,mpower,npower,normp,coef,alpha,
       nat,atno,x,y,z,V);

  getInts(nbf,nprim,normc,istart,xcenter,ycenter,zcenter, 
	  lpower,mpower,npower,normp,coef,alpha,Ints);

  /* Form h */
  for (i=0; i<nbf*nbf; i++) h[i] = T[i] + V[i];

  /* Form the density matrix from the guess orbs */
  getD(nbf,nocc,orbs,D);

  /* diagonalize h to get guess orbs */
  for (i=0; i<nbf*nbf; i++) {
    orbs[i] = h[i];
    H[i] = h[i];
    temp[i] = S[i];
  }
  info = gjacobi(H,temp,nbf,orbe,orbs,temp1,temp2);

  for (iter=0; iter<maxiter; iter++){
    /* Form the density matrix from the guess orbs */
    getD(nbf,nocc,orbs,D);

    /* Form the J and K matrices */
    getJ(nbf,D,Ints,J);
    getK(nbf,D,Ints,K);

    /* Form the Fock matrix in the orbs space */
    for (i=0; i<nbf*nbf; i++) {
      H[i] = h[i] + 2*J[i] - K[i];
      temp[i] = S[i];
    };
    /* Diagonalize F */
    info = gjacobi(H,temp,nbf,orbe,orbs,temp1,temp2);

    /* Compute the HF energy */
    eone = 0.;
    for (i=0; i<nbf*nbf; i++) eone += h[i]*D[i];
    ej = 0.;
    for (i=0; i<nbf*nbf; i++) ej += D[i]*J[i];
    ek = 0.;
    for (i=0; i<nbf*nbf; i++) ek += D[i]*K[i];
    energy = enuke+2*eone+2*ej-ek;
    printf("%%d %%f %%f %%f %%f %%f\n",iter,energy,enuke,eone,ej,ek);

  }
  return 0;
}
