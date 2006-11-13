void getS(int nbf, int *nprim, double *normc, int *istart,
	  double *xcenter, double *ycenter, double *zcenter, 
	  int *lpower,int *mpower, int *npower, 
	  double *normp, double *coef, double *alpha, 
	  double *S);

void getT(int nbf, int *nprim, double *normc, int *istart,
	  double *xcenter, double *ycenter, double *zcenter, 
	  int *lpower,int *mpower, int *npower, 
	  double *normp, double *coef, double *alpha, 
	  double *T);

void getV(int nbf, int *nprim, double *normc, int *istart,
	  double *xcenter, double *ycenter, double *zcenter, 
	  int *lpower,int *mpower, int *npower, 
	  double *normp, double *coef, double *alpha, 
	  int nat, int *atno, double *x, double *y, double *z,
	  double *V);

void getInts(int nbf, int *nprim, double *normc, int *istart,
	     double *xcenter, double *ycenter, double *zcenter, 
	     int *lpower,int *mpower, int *npower, 
	     double *normp, double *coef, double *alpha, double *Ints);

void getD(int nbf,int nocc,double *orbs,double *D);
void getJ(int nbf,double *D,double *Ints,double *J);
void getK(int nbf,double *D,double *Ints,double *K);
//void printmat(int nbf, double *A, char *tag, int ilow, int ihi, 
//	      int jlow, int jhi);
