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
  int nat = 2;
  int nbf = 10;
  int nocc=1,info;
  int iter,maxiter=5;
  double h[100],S[100],T[100],V[100],Ints[1540];
  double D[100],J[100],K[100];
  double orbe[10],orbs[100],temp[100];
  double temp1[10], temp2[10];
  double H[100];
  double energy,eone,ej,ek;
  double enuke=1.428571;

  /*  Data over atoms  */
  int atno[] = {
    1,
    1,
  };
  double x[] = {
    0.0,
    0.0,
  };
  double y[] = {
    0.0,
    0.0,
  };
  double z[] = {
    -0.35,
    0.35,
  };

  /*  Data over contracted bfns  */
  double xcenter[] = {
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
  };  
  double ycenter[] = {
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
  };  
  double zcenter[] = {
    -0.35,
    -0.35,
    -0.35,
    -0.35,
    -0.35,
    0.35,
    0.35,
    0.35,
    0.35,
    0.35,
  };  
  int lpower[] = {
    0,
    0,
    1,
    0,
    0,
    0,
    0,
    1,
    0,
    0,
  };
  int mpower[] = {
    0,
    0,
    0,
    1,
    0,
    0,
    0,
    0,
    1,
    0,
  };
  int npower[] = {
    0,
    0,
    0,
    0,
    1,
    0,
    0,
    0,
    0,
    1,
  };
  double normc[] = {
    1.00000010954,
    1.0,
    1.0,
    1.0,
    1.0,
    1.00000010954,
    1.0,
    1.0,
    1.0,
    1.0,
  };
  int istart[] = {
    0,
    3,
    4,
    5,
    6,
    7,
    10,
    11,
    12,
    13,
  };
  int nprim[] = {
    3,
    1,
    1,
    1,
    1,
    3,
    1,
    1,
    1,
    1,
  };

  /*  Data over primitive bfns  */
  double alpha[] = {
    18.731137,
    2.825394,
    0.640122,
    0.161278,
    1.1,
    1.1,
    1.1,
    18.731137,
    2.825394,
    0.640122,
    0.161278,
    1.1,
    1.1,
    1.1,
  };
  double coef[] = {
    0.033495,
    0.234727,
    0.813757,
    1.0,
    1.0,
    1.0,
    1.0,
    0.033495,
    0.234727,
    0.813757,
    1.0,
    1.0,
    1.0,
    1.0,
  };
  double normp[] = {
    6.4170171128,
    1.55317129726,
    0.510043429599,
    0.181380852627,
    1.60576114265,
    1.60576114265,
    1.60576114265,
    6.4170171128,
    1.55317129726,
    0.510043429599,
    0.181380852627,
    1.60576114265,
    1.60576114265,
    1.60576114265,
  };

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
    printf("%d %f %f %f %f %f\n",iter,energy,enuke,eone,ej,ek);

  }
  return 0;
}
