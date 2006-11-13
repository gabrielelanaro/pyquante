#include <stdio.h>
#include <stdlib.h>
#include "cints.h"

void getS(int nbf, int *nprim, double *normc, int *istart,
	  double *xcenter, double *ycenter, double *zcenter, 
	  int *lpower,int *mpower, int *npower, 
	  double *normp, double *coef, double *alpha, 
	  double *S){
  int i,j,iprim,jprim,jindex,iindex;
  for (j=0; j<nbf; j++){
    for (i=0; i<nbf; i++){
      S[i+j*nbf] = 0.;
      for (jprim=0; jprim<nprim[j]; jprim++){
	jindex = istart[j]+jprim;
	for (iprim=0; iprim<nprim[i]; iprim++){
	  iindex = istart[i]+iprim;
	  S[i+j*nbf] += normc[i]*normc[j]*normp[iindex]
	    *normp[jindex]*coef[iindex]*coef[jindex]
	    *overlap(alpha[iindex],lpower[i],mpower[i],
		     npower[i],xcenter[i],ycenter[i],
		     zcenter[i],
		     alpha[jindex],lpower[j],mpower[j],
		     npower[j],xcenter[j],ycenter[j],
		     zcenter[j]);
	}
      }
    }
  }
}

void getT(int nbf, int *nprim, double *normc, int *istart,
	  double *xcenter, double *ycenter, double *zcenter, 
	  int *lpower,int *mpower, int *npower, 
	  double *normp, double *coef, double *alpha, 
	  double *T){
  int i,j,iprim,jprim,jindex,iindex;
  for (j=0; j<nbf; j++){
    for (i=0; i<nbf; i++){
      T[i+j*nbf] = 0.;
      for (jprim=0; jprim<nprim[j]; jprim++){
	jindex = istart[j]+jprim;
	for (iprim=0; iprim<nprim[i]; iprim++){
	  iindex = istart[i]+iprim;
	  T[i+j*nbf] += normc[i]*normc[j]*normp[iindex]
	    *normp[jindex]*coef[iindex]*coef[jindex]
	    *kinetic(alpha[iindex],lpower[i],mpower[i],
		     npower[i],xcenter[i],ycenter[i],
		     zcenter[i],
		     alpha[jindex],lpower[j],mpower[j],
		     npower[j],xcenter[j],ycenter[j],
		     zcenter[j]);
	}
      }
    }
  }
}


void getV(int nbf, int *nprim, double *normc, int *istart,
	  double *xcenter, double *ycenter, double *zcenter, 
	  int *lpower,int *mpower, int *npower, 
	  double *normp, double *coef, double *alpha, 
	  int nat, int *atno, double *x, double *y, double *z,
	  double *V){
  int i,j,iprim,jprim,jindex,iindex,iat;
  for (j=0; j<nbf; j++){
    for (i=0; i<nbf; i++){
      V[i+j*nbf] = 0.;
      for (jprim=0; jprim<nprim[j]; jprim++){
	jindex = istart[j]+jprim;
	for (iprim=0; iprim<nprim[i]; iprim++){
	  iindex = istart[i]+iprim;
	  for (iat=0; iat<nat; iat++){
	    V[i+j*nbf] += normc[i]*normc[j]*coef[iindex]*coef[jindex]
              *atno[iat]
	      *nuclear_attraction(xcenter[i],ycenter[i],
                                  zcenter[i],normp[iindex],
				  lpower[i],mpower[i],
                                  npower[i],alpha[iindex],
				  xcenter[j],ycenter[j],
                                  zcenter[j],normp[jindex],
				  lpower[j],mpower[j],
                                  npower[j],alpha[jindex],
				  x[iat],y[iat],z[iat]);
	  }
	}
      }
    }
  }
}

void getInts(int nbf, int *nprim, double *normc, int *istart,
	     double *xcenter, double *ycenter, double *zcenter, 
	     int *lpower,int *mpower, int *npower, 
	     double *normp, double *coef, double *alpha, double *Ints){
  int i,j,k,l;
  int i_istart, j_istart, k_istart, l_istart;
  int ij, kl, ijkl;
  for (i=0; i<nbf; i++){
    i_istart = istart[i];
    for (j=0; j<i+1; j++){
      j_istart = istart[j];
      ij = i*(i+1)/2+j;
      for (k=0; k<nbf; k++){
	k_istart = istart[k];
	for (l=0; l<k+1; l++){
	  l_istart = istart[l];
	  kl = k*(k+1)/2+l;
	  if (ij >= kl) {
	    ijkl = ijkl2intindex(i,j,k,l);
	    Ints[ijkl] = normc[i]*normc[j]*normc[k]*normc[l]*
	      contr_coulomb(nprim[i],&alpha[i_istart],&coef[i_istart],
                            &normp[i_istart],
			    xcenter[i],ycenter[i],zcenter[i],
			    lpower[i],mpower[i],npower[i],
			    nprim[j],&alpha[j_istart],&coef[j_istart],
                            &normp[j_istart],
			    xcenter[j],ycenter[j],zcenter[j],
			    lpower[j],mpower[j],npower[j],
			    nprim[k],&alpha[k_istart],&coef[k_istart],
                            &normp[k_istart],
			    xcenter[k],ycenter[k],zcenter[k],
			    lpower[k],mpower[k],npower[k],
			    nprim[l],&alpha[l_istart],&coef[l_istart],
                            &normp[l_istart],
			    xcenter[l],ycenter[l],zcenter[l],
			    lpower[l],mpower[l],npower[l]);
	  }
	}
      }
    }
  }
}

void getD(int nbf,int nocc,double *orbs,double *D){
  int iorb,ibf,jbf;
  double val;
  for (ibf=0; ibf<nbf; ibf++){
    for (jbf=0; jbf<nbf; jbf++){
      val = 0;
      for (iorb=0; iorb<nocc; iorb++)
	val += orbs[ibf+iorb*nbf]*orbs[jbf+iorb*nbf];
      D[ibf+jbf*nbf] = val;
    }
  }
}

void getJ(int nbf,double *D,double *Ints,double *J){
  int i,j,k,l,kl,ijkl;
  double *tmp;
  double val=0.;
  tmp = (double *)malloc(nbf*nbf*sizeof(double));
  for (i=0; i<nbf; i++){
    for (j=0; j<i+1; j++){
      kl = 0;
      for (k=0; k<nbf; k++){
	for (l=0; l<nbf; l++){
	  ijkl = ijkl2intindex(i,j,k,l);
	  tmp[kl] = Ints[ijkl];
	  kl += 1;
	}
      }
      val = 0;
      for (kl=0; kl<nbf*nbf; kl++) val += tmp[kl]*D[kl];
      J[i+j*nbf] = val;
      J[j+i*nbf] = val;
    }
  }
  free(tmp);
}

void getK(int nbf,double *D,double *Ints,double *K){
  int i,j,k,l,kl,ikjl,iljk;
  double *tmp;
  double val = 0.;
  tmp = (double *)malloc(nbf*nbf*sizeof(double));
  for (i=0; i<nbf; i++){
    for (j=0; j<i+1; j++){
      kl = 0;
      for (k=0; k<nbf; k++){
	for (l=0; l<nbf; l++){
	  ikjl = ijkl2intindex(i,k,j,l);
	  iljk = ijkl2intindex(i,l,j,k);
	  tmp[kl] = 0.5*(Ints[ikjl]+Ints[iljk]);
	  kl += 1;
	}
      }
      val = 0;
      for (kl=0; kl<nbf*nbf; kl++) val += tmp[kl]*D[kl];
      K[i+j*nbf] = val;
      K[j+i*nbf] = val;
    }
  }
  free(tmp);
}
/*
void printmat(int nbf, double *A, char *tag, int ilow, int ihi, 
	      int jlow, int jhi){
  int i,j;
  for (i=ilow; i<ihi; i++){
    for (j=jlow; j<jhi; j++){
      printf("%s[%d,%d] = %f\n",tag,i,j,A[i+j*nbf]);
    }
  }
}
*/
