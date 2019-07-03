#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "api_scilab.h"
#include "Scierror.h"
#include "BOOL.h"
#include "localization.h"


//[xopt,fopt] = qcqp(x,H,f,A,b,Aeq,beq,Q,c,r,lb,ub)

struct prob
{
/* data */
//initial point : x
double * x;
//Objective function : x'Hx+f'x
double ** H;
double * f;
//No of variables
int n;

//Linear inequality constraint Ax <= b
double ** A; // m x n
double * b; // m x 1
//No of Linear inequality Constraints
int m;

//Linear equality constraint Aeqx = b
double ** Aeq;  
double * beq;
//No of Linear equality constraint 
int p;

//Quadratic inequality constraint x'Qx + c'x <= r
double *** Q; 
double ** c;
double *r; 
//No of Quadratic Constraint
int q;

//Lower and Upper bounds
double * lb;
double * ub;

};
struct prob problem;

void c_algencan(void *myevalf, void *myevalg, void *myevalh, void *myevalc,
	void *myevaljac, void *myevalhc, void *myevalfc, void *myevalgjac,
	void *myevalgjacp, void *myevalhl, void *myevalhlp, int jcnnzmax,
       	int hnnzmax,double *epsfeas, double *epsopt, double *efstin,
	double *eostin, double *efacc, double *eoacc, char *outputfnm,
	char *specfnm, int nvparam,char **vparam, int n, double *x,
	double *l, double *u, int m, double *lambda, _Bool *equatn,
	_Bool *linear, _Bool *coded, _Bool checkder, double *f,
	double *cnorm, double *snorm, double *nlpsupn,int *inform);


/**********************************************************************
   sci_gateway
************************************************************************/
  
static const char fname[] = "qcqp";
/* ==================================================================== */




/* ******************************************************************
    Objective Function. Must be Quadratic.
    f(x) = x'.H.x + f'.x 
   ****************************************************************** */

void myevalf(int n, double *x, double *f_, int *flag) {
   int i = 0;
   int j = 0;

   *flag = 0;
   *f_ = 0;

   for( i = 0 ; i < n ; i++ )
   {
      for( j = 0 ; j < n ; j++ ){
         *f_ += x[i]*problem.H[i][j]*x[j];
      }
   }   

   for( i = 0 ; i < n ; i++ ){
      *f_ += x[i]*problem.f[i];
   }

}

/* ******************************************************************
    Gradient of Objective function.
    f'(x)  
   ****************************************************************** */

void myevalg(int n, double *x, double *g, int *flag) {
   int i = 0;
   int j = 0;
   *flag = 0;
   for( i = 0 ; i < n ; i++ ){
      g[i] = 0;
      for( j = 0 ; j < n ; j++ ){
            g[i] += problem.H[i][j]*x[j];
      }
      g[i] += problem.H[i][i]*x[i] + problem.f[i]; 
   }
}

/* ******************************************************************
    Hessian of Objective function.
    f''(x) 
   ****************************************************************** */

void myevalh(int n, double *x, int *hrow, int *hcol, double *hval, int *hnnz,
	     int lim, _Bool *lmem, int *flag) {

   *flag = 0;
   *lmem = 0;

   int i = 0;
   int j = 0;

   *hnnz = (n*(n+1))/2;
   if( *hnnz > lim ) {
     *lmem = 1;
     return;
   }
   int temp = 0;
   for( i =  0 ; i < n ; i++ ){
      for(j = 0 ; j <= i ; j++ ){
         hrow[temp] = i;
         hcol[temp] = j;
         hval[temp] = 2*problem.H[i][j];
         temp++;
      }
   }
}

/* ******************************************************************
    Equality and Inequality Constraints
   ****************************************************************** */

void myevalc(int n, double *x, int ind, double *c, int *flag) {

   *flag = 0;
   int i = 0;
   int j = 0;
   *c = 0;
   if(ind < problem.m){
      for( i = 0 ; i < n ; i++ ){
         *c += problem.A[ind][i] * x[i];
      } *c -= problem.b[ind];
   }else if(ind < problem.m + problem.p){
      for( i = 0 ; i < n ; i++ ){
         *c += problem.Aeq[ind-problem.m][i] * x[i];
      } *c = problem.beq[ind-problem.m];
   }else{
      for( i = 0 ; i < n ; i++ )
      {
         for( j = 0 ; j < n ; j++ ){
            *c += x[i]*problem.Q[ind - problem.m - problem.p][i][j]*x[j];
         }
      }   

      for( i = 0 ; i < n ; i++ ){
         *c += x[i]*problem.f[ind - problem.m - problem.p];
      }
      *c -= problem.r[ind - problem.m - problem.p];
   }
}

/* ******************************************************************
    Gradient of Constraints
   ****************************************************************** */

void myevaljac(int n, double *x, int ind, int *jcvar, double *jcval,
	       int *jcnnz, int lim, _Bool *lmem, int *flag) {
  *flag = 0;
  *lmem = 0;
   int i = 0;
   int j = 0;
   
   *jcnnz = n;
   if( *jcnnz > lim ) {
     *lmem = 1;
     return;
   }

   if(ind < problem.m){
      for( i = 0 ; i < n ; i++ ){
         jcvar[i] = i;
         jcval[i] = problem.A[ind][i];
      }
   }else if(ind < problem.p+problem.m){
      for( i = 0 ; i < n ; i++ ){
         jcvar[i]=i;
         jcval[i] = problem.Aeq[ind-problem.m][i];
      }
   }else{
      for( i = 0 ; i < n ; i++ ){
         jcval[i] = 0;
         jcvar[i]=i;
         for( j = 0 ; j < n ; j++ ){
               jcval[i] += problem.Q[ind-problem.m-problem.p][i][j]*x[j];
         }
         jcval[i] += problem.Q[ind-problem.m-problem.p][i][i]*x[i] + problem.c[ind-problem.m-problem.p][i]; 
      }
   }
}

/* ******************************************************************
    Hessian of Constraints
   ****************************************************************** */

void myevalhc(int n, double *x, int ind, int *hcrow, int *hccol, double *hcval,
	      int *hcnnz, int lim, _Bool *lmem, int *flag) {
  *flag = 0;
  *lmem = 0;
   int i = 0;
   int j = 0;
   int temp = 0;
   
   
   if(ind < problem.m + problem.p){
      *hcnnz = 0;
   }else{
      *hcnnz = (n*(n+1))/2;
      if( *hcnnz > lim ) {
        *lmem = 1;
         return;
      }
      for( int i =  0 ; i < n ; i++ ){
         for( int j = 0 ; j <= i ; j++ ){
            hcrow[temp] = i;
            hccol[temp] = j;
            hcval[temp] = 2*problem.Q[ind-problem.m-problem.p][i][j];
            temp++;
         }
      }  
   }  
}

/* *****************************************************************
   ***************************************************************** */

void myevalfc(int n, double *x, double *f, int m, double *c, int *flag) {

   *flag = -1;
}

/* *****************************************************************
   ***************************************************************** */

void myevalgjac(int n, double *x, double *g, int m, int *jcfun, int *jcvar,
		double *jcval, int *jcnnz, int lim, _Bool *lmem, int *flag) {

   *flag = -1;
}

/* *****************************************************************
   ***************************************************************** */

void myevalgjacp(int n, double *x, double *g, int m, double *p, double *q,
		 char work, _Bool *gotj, int *flag) {

   *flag = -1;
}

/* *****************************************************************
   ***************************************************************** */

void myevalhl(int n, double *x, int m, double *lambda, double scalef,
	      double *scalec, int *hlrow, int *hlcol, double *hlval,
	      int *hlnnz, int lim, _Bool *lmem, int *flag) {

   *flag = -1;
}

/* *****************************************************************
   ***************************************************************** */

void myevalhlp(int n, double *x, int m, double *lambda, double scalef,
	       double *scalec, double *p, double *hp, _Bool *goth, 
	       int *flag) {

   *flag = -1;
}

/***************************************************************************************************
sci_gateway
***************************************************************************************************/

int sci_qcqp(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt* opt, int nout, scilabVar* out)
{
  _Bool  checkder;
  int    hnnzmax,jcnnzmax,inform,nvparam,ncomp,m;
  double cnorm,efacc,efstin,eoacc,eostin,epsfeas,epsopt,f,nlpsupn,snorm;
  int count = 0;
  
  char   *specfnm, *outputfnm, **vparam;
  _Bool  coded[11],*equatn,*linear;
  double *l,*lambda,*u;

   int i = 0;
   int j = 0;
   int k = 0;
   int temp = 0;
   int temp1 = 0;
   double * Qtemp;

   if(nin!=12){
      Scierror(999,"%s : Wrong number of inputs, %d expected",fname,12);
      return 1;
   }
    
    if(nout!=4){
        Scierror(999,"%s : Wrong number of output argument, %d expected",fname,4);
        return 1;
    }

   for(i = 0 ; i < 12 ; i++){
      if( (scilab_isDouble(env, in[i]) == 0 || scilab_isMatrix2d(env, in[i]) == 0) && i != 7  ) {
        Scierror(999,"%s : Wrong type for input argument %d, A double Matrix expected",fname,i+1);
        return 1;
      }
   }
    
    
    scilab_getDim2d(env,in[0],&problem.n,&temp);
    scilab_getDoubleArray(env,in[0],&problem.x);
    
    scilab_getDoubleArray(env,in[1],&Qtemp);

   
    problem.H = (double **)malloc(sizeof(double*)*problem.n);
    for(i = 0 ; i < problem.n ; i++){
       problem.H[i] = (double *)malloc(sizeof(double)*problem.n);
    }
    
    for(i = 0 ; i < problem.n ; i++){
        for(j = 0 ; j < problem.n ; j++){
            problem.H[j][i] = Qtemp[j+problem.n*i];
        }
    }
    
   scilab_getDoubleArray(env,in[2],&problem.f);
    
   scilab_getDim2d(env,in[3],&problem.m,&temp);
   scilab_getDoubleArray(env,in[3],&Qtemp);

   problem.A = (double **)malloc(sizeof(double*)*problem.m);
   for(i = 0 ; i < problem.m ; i++){
      problem.A[i] = (double *)malloc(sizeof(double)*problem.n);
   }
   scilab_getDoubleArray(env,in[3],&Qtemp);
   
   for(i = 0 ; i < problem.n ; i++){
      for(j = 0 ; j < problem.m ; j++){
         problem.A[j][i] = Qtemp[j+problem.m*i];
      }
   }
   scilab_getDoubleArray(env,in[4],&problem.b);
   
   scilab_getDim2d(env,in[5],&problem.p,&temp);
   problem.Aeq = (double **)malloc(sizeof(double*)*problem.p);
   for(i = 0 ; i < problem.p ; i++){
       problem.Aeq[i] = (double *)malloc(sizeof(double)*problem.n);
   }
   scilab_getDoubleArray(env,in[5],&Qtemp);
   for(i = 0 ; i < problem.n ; i++){
      for(j = 0 ; j < problem.p ; j++){
         problem.Aeq[j][i] = Qtemp[j+problem.p*i];
      }
   }
   scilab_getDoubleArray(env,in[6],&problem.beq);

   scilab_getDim2d(env,in[8],&problem.q,&temp);
   scilab_getDoubleArray(env,in[7],&Qtemp);
   
   problem.Q = (double ***)malloc(sizeof(double**)*problem.q);
    for(i = 0 ; i < problem.q ; i++){
      problem.Q[i] = (double **)malloc(sizeof(double*)*problem.n);
    }

    for(i = 0 ; i < problem.q ; i++){
       for(j = 0 ; j < problem.n ; j++){
          problem.Q[i][j] = (double *)malloc(sizeof(double)*problem.n);
      }
    }
   
   for(i = 0 ; i < problem.q ; i++){
      for(j = 0 ; j < problem.n ; j++){
         for(k = 0 ; k < problem.n ; k++){
            problem.Q[i][k][j] = Qtemp[k+problem.n*(j+problem.n*i) ];
         }
      }
    }
    scilab_getDoubleArray(env,in[8],&Qtemp);
   
    problem.c = (double **)malloc(sizeof(double*)*problem.q);
    for(i = 0 ; i < problem.q ; i++){
       problem.c[i] = (double *)malloc(sizeof(double)*problem.n);
    }
    for(i = 0 ; i < problem.n ; i++){
      for(j = 0 ; j < problem.q ; j++){
         problem.c[j][i] = Qtemp[j+problem.q*i];
      }
    }
    
    scilab_getDoubleArray(env,in[9],&problem.r);
    scilab_getDoubleArray(env,in[10],&problem.lb);
    scilab_getDoubleArray(env,in[11],&problem.ub);
   
   m = problem.m + problem.p + problem.q;
  /* Memory allocation */
  lambda = (double *) malloc(m * sizeof(double));
  equatn = (_Bool  *) malloc(m * sizeof(_Bool ));
  linear = (_Bool  *) malloc(m * sizeof(_Bool ));
    
  /* For each constraint i, set equatn[i] = 1. if it is an equality
     constraint of the form c_i(x) = 0, and set equatn[i] = 0 if it is
     an inequality constraint of the form c_i(x) <= 0. */
  for( i = 0 ; i < m ; i++ ){
     equatn[i] = (i - problem.m) < problem.p && i >= problem.m ? 1 : 0;
   }

  /* For each constraint i, set linear[i] = 1 if it is a linear
     constraint, otherwise set linear[i] = 0 */
  for( i = 0 ; i < m ; i++ ){
     linear[i] = i < problem.m + problem.p ? 1 : 0;
   }
  // Lagrange multipliers approximation. 
  for( i = 0; i < m; i++ ) lambda[i] = 0.0;
  
   // In this C interface evalf, evalg, evalh, evalc, evaljac and
   //   evalhc are present. evalfc, evalgjac, evalhl and evalhlp are
   //   not. 
  
  coded[0]  = 1; // fsub     
  coded[1]  = 1; // gsub     
  coded[2]  = 1; // hsub     
  coded[3]  = 1; // csub     
  coded[4]  = 1; // jacsub   
  coded[5]  = 1; // hcsub    
  coded[6]  = 0; // fcsub    
  coded[7]  = 0; // gjacsub  
  coded[8]  = 0; // gjacpsub 
  coded[9]  = 0; // hlsub    
  coded[10] = 0; // hlpsub   
 
  // Upper bounds on the number of sparse-matrices non-null elements 
  jcnnzmax = 2*problem.n*problem.n;
  hnnzmax  = 2*problem.n*problem.n;

  // Check derivatives? 
  checkder = 0;

  // Parameters setting 
  epsfeas = 1.0e-08;
  epsopt  = 1.0e-08;
  efstin  = sqrt( epsfeas );
  eostin  = pow( epsopt, 1.5 );
  efacc   = sqrt( epsfeas );
  eoacc   = sqrt( epsopt );

  outputfnm = "algencan.out";
  specfnm   = "";

  nvparam = 1;
  
  // Allocates VPARAM array 
  vparam = ( char ** ) malloc( nvparam * sizeof( char * ) );

//   Set algencan parameters 
  vparam[0] = "ITERATIONS-OUTPUT-DETAIL 10";
   
   c_algencan(&myevalf,&myevalg,&myevalh,&myevalc,&myevaljac,&myevalhc,&myevalfc,
	     &myevalgjac,&myevalgjacp,&myevalhl,&myevalhlp,jcnnzmax,hnnzmax,
	     &epsfeas,&epsopt,&efstin,&eostin,&efacc,&eoacc,outputfnm,specfnm,
	     nvparam,vparam,problem.n,problem.x,problem.lb,problem.ub,m,lambda,equatn,linear,coded,checkder,
	     &f,&cnorm,&snorm,&nlpsupn,&inform);

   out[0] = scilab_createDoubleMatrix2d(env, problem.n, 1, 0);
   scilab_setDoubleArray(env,out[0],problem.x);
   out[1] = scilab_createDouble(env, f);
   out[2] = scilab_createDoubleMatrix2d(env,m,1,0);
   scilab_setDoubleArray(env,out[2],lambda);
   out[3] = scilab_createDouble(env, inform);

   return 0;
}




// /* ******************************************************************
      // Solving example using c_algencan
      // Reference : http://www.minlplib.org/nvs11.p1.html
//    ****************************************************************** */

// int main() {
// int sci_qcqp(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt* opt, int nout, scilabVar* out){
//    _Bool  checkder;
//   int    hnnzmax,jcnnzmax,inform,nvparam,ncomp,m;
//   double cnorm,efacc,efstin,eoacc,eostin,epsfeas,epsopt,f,nlpsupn,snorm;
  
//   char   *specfnm, *outputfnm, **vparam;
//   _Bool  coded[11],*equatn,*linear;
//   double *l,*lambda,*u;

//    int i;
//    int j;
//    int k;
//    int temp;
//    int temp1;
//    double * Qtemp;
   
//    problem.n = 3;
//    problem.m = 0;
//    problem.p = 0;
//    problem.q = 3;

//     m = problem.m + problem.p + problem.q;
//   /* Memory allocation */
//   problem.x = (double *) malloc(problem.n * sizeof(double)); 
//   lambda = (double *) malloc(m * sizeof(double));
//   equatn = (_Bool  *) malloc(m * sizeof(_Bool ));
//   linear = (_Bool  *) malloc(m * sizeof(_Bool ));
  
//      problem.H = (double * []){(double[]){7.0,0.0,-3.0},(double[]){0.0,6.0,2.0},(double[]){-3.0,2.0,8.0}};
//    problem.f = (double[]){-15.8,-93.2,-63.0};

//    problem.lb = (double[]){0.0,0.0,0.0};
//    problem.ub = (double[]){200.0,200.0,200.0};
//    problem.Q = (double**[]){(double*[]){(double[]){9.0,5.0,3.0},
//                                        (double[]){5.0,8.0,5.0},
//                                        (double[]){3.0,5.0,5.0}},
//                            (double*[]){(double[]){6.0,4.0,1.0},
//                                        (double[]){4.0,6.0,1.0},
//                                        (double[]){1.0,1.0,4.0}},
//                            (double*[]){(double[]){9.0,-1.0,0.0},
//                                        (double[]){-1.0,6.0,-1},
//                                        (double[]){0.0,-1.0,8.0}}};
//    problem.c = (double*[]){(double[]){0.0,0.0,0.0},(double[]){0.0,0.0,0,0},(double[]){0.0,0.0,0.0}};
//    problem.r = (double[]){1000.0,550.0,340.0};

//   /* Initial point */
//   //for(i = 0; i < problem.n; i++) problem.x[i] = 1.0;
//   problem.x[0] = 2.18;
//   problem.x[1] = 6.94;
//   problem.x[2] = 2.99;
    
//   /* For each constraint i, set equatn[i] = 1. if it is an equality
//      constraint of the form c_i(x) = 0, and set equatn[i] = 0 if it is
//      an inequality constraint of the form c_i(x) <= 0. */
//   for( i = 0 ; i < m ; i++ ){
//      equatn[i] = (i - problem.m) < problem.p && i >= problem.m ? 1 : 0;
//    }

//   /* For each constraint i, set linear[i] = 1 if it is a linear
//      constraint, otherwise set linear[i] = 0 */
//   for( i = 0 ; i < m ; i++ ){
//      linear[i] = i - problem.m + problem.p < 1 ? 1 : 0;
//    }
//   /* Lagrange multipliers approximation. */
//   for( i = 0; i < m; i++ ) lambda[i] = 0.0;
  
//   /* In this C interface evalf, evalg, evalh, evalc, evaljac and
//      evalhc are present. evalfc, evalgjac, evalhl and evalhlp are
//      not. */
  
//   coded[0]  = 1; /* fsub     */
//   coded[1]  = 1; /* gsub     */
//   coded[2]  = 1; /* hsub     */
//   coded[3]  = 1; /* csub     */
//   coded[4]  = 1; /* jacsub   */
//   coded[5]  = 1; /* hcsub    */
//   coded[6]  = 0; /* fcsub    */
//   coded[7]  = 0; /* gjacsub  */
//   coded[8]  = 0; /* gjacpsub */
//   coded[9]  = 0; /* hlsub    */
//   coded[10] = 0; /* hlpsub   */
 
//   /* Upper bounds on the number of sparse-matrices non-null
//      elements */
//   jcnnzmax = 2*problem.n*problem.n;
//   hnnzmax  = 2*problem.n*problem.n;


//   /* Check derivatives? */
//   checkder = 0;

//   /* Parameters setting */
//   epsfeas = 1.0e-08; 
//   epsopt  = 1.0e-08;
//   efstin  = sqrt( epsfeas );
//   eostin  = pow( epsopt, 1.5 );
//   efacc   = sqrt( epsfeas );
//   eoacc   = sqrt( epsopt );

//   outputfnm = "algencan.out";
//   specfnm   = "";

//   nvparam = 1;
  
//   /* Allocates VPARAM array */
//   vparam = ( char ** ) malloc( nvparam * sizeof( char * ) );

//   /* Set algencan parameters */
//   vparam[0] = "ITERATIONS-OUTPUT-DETAIL 10";

//    c_algencan(&myevalf,&myevalg,&myevalh,&myevalc,&myevaljac,&myevalhc,&myevalfc,
// 	     &myevalgjac,&myevalgjacp,&myevalhl,&myevalhlp,jcnnzmax,hnnzmax,
// 	     &epsfeas,&epsopt,&efstin,&eostin,&efacc,&eoacc,outputfnm,specfnm,
// 	     nvparam,vparam,problem.n,problem.x,problem.lb,problem.ub,m,lambda,equatn,linear,coded,checkder,
// 	     &f,&cnorm,&snorm,&nlpsupn,&inform);


//   // scilab_setDoubleArray(env,out[0],&problem.x);
//    //scilab_setDoubleArray(env,out[1],&f);

//   /* Memory deallocation */
//   printf("\nOptimal Solution at ");
//   for(i = 0 ; i < problem.n ; i++){
//      printf("%f ",problem.x[i]);
//   }
//   printf("\n");

//   free(lambda);
//   free(equatn);
//   free(linear);
//   return 0;
// }