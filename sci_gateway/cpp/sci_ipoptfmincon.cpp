// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: R.Vidyadhar & Vignesh Kannan
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#include "IpIpoptApplication.hpp"
#include "minconNLP.hpp"
#include <IpSolveStatistics.hpp>

extern "C"
{
#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <localization.h>
#include <sciprint.h>
#include <wchar.h>


using namespace std;

const char fname[] = "solveminconp";
/* ==================================================================== */
int sci_solveminconp(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out) 
{

	using namespace Ipopt;


	//Function pointers, input matrix(Starting point) pointer, flag variable 
	double *x0ptr=NULL, *lbptr=NULL, *ubptr=NULL,*Aptr=NULL, *bptr=NULL, *Aeqptr=NULL, *beqptr=NULL;
	double nonlinCon=0,nonlinIneqCon=0;
        

    	// Input arguments
	static unsigned int nVars = 0,nCons = 0;
	int x0_rows=0, x0_cols=0, lb_rows=0, lb_cols=0, ub_rows=0, ub_cols=0, A_rows=0, A_cols=0, b_rows=0, b_cols=0, Aeq_rows=0, Aeq_cols=0, beq_rows=0, beq_cols=0;
	
	// Output arguments	
	double  ObjVal=0,iteration=0,cpuTime=0,fobj_eval=0;
	double dual_inf, constr_viol, complementarity, kkt_error;
	const double *fX = NULL, *fGrad =  NULL;
	const double *fHess =  NULL;
	const double *fLambda = NULL;
	const double *fZl=NULL;
	const double *fZu=NULL;
	int rstatus = 0;
	int int_fobj_eval, int_constr_eval, int_fobj_grad_eval, int_constr_jac_eval, int_hess_eval;

	if (nin !=16)  //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname, 16);
        	return STATUS_ERROR; 
	}
	
	if (nout !=12) //Checking the output arguments

	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname, 12);
		return 1;
	}

	////////// Manage the input argument //////////

	
	//Getting matrix representing linear inequality constraints 

	if (scilab_isDouble(env, in[1]) == 0 || scilab_isMatrix2d(env, in[1]) == 0)
	{
		Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 2);
		return 1;
	}
	
	scilab_getDoubleArray(env, in[1], &Aptr);
	int size1 = scilab_getDim2d(env, in[1], &A_rows, &A_cols);

	//Getting matrix representing bounds of linear inequality constraints 
	
	if (scilab_isDouble(env, in[2]) == 0 || scilab_isMatrix2d(env, in[2]) == 0)
	{
		Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 3);
		return 1;
	}
	
	scilab_getDoubleArray(env, in[2], &bptr);
	size1 = scilab_getDim2d(env, in[2], &b_rows, &b_cols);
		
	//Getting matrix representing linear equality constraints 
	if (scilab_isDouble(env, in[3]) == 0 || scilab_isMatrix2d(env, in[3]) == 0)
	{
		Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 4);
		return 1;
	}
	
	scilab_getDoubleArray(env, in[3], &Aeqptr);
	size1 = scilab_getDim2d(env, in[3], &Aeq_rows, &Aeq_cols);

	//Getting matrix representing bounds of linear inequality constraints 
	if (scilab_isDouble(env, in[4]) == 0 || scilab_isMatrix2d(env, in[4]) == 0)
	{
		Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 5);
		return 1;
	}
	
	scilab_getDoubleArray(env, in[4], &beqptr);
	size1 = scilab_getDim2d(env, in[4], &beq_rows, &beq_cols);

	//Getting matrix representing linear inequality constraints 
	if (scilab_isDouble(env, in[5]) == 0 || scilab_isMatrix2d(env, in[5]) == 0)
	{
		Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 6);
		return 1;
	}
	
	scilab_getDoubleArray(env, in[5], &lbptr);

	//Getting matrix representing linear inequality constraints 
	if (scilab_isDouble(env, in[6]) == 0 || scilab_isMatrix2d(env, in[6]) == 0)
	{
		Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 7);
		return 1;
	}
	
	scilab_getDoubleArray(env, in[6], &ubptr);

	//Number of non-linear constraints
	if (scilab_isDouble(env, in[7]) == 0 || scilab_isScalar(env, in[7]) == 0)
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: An int expected.\n", fname, 8);
    	return 1;
	}
	
	scilab_getDouble(env, in[7], &nonlinCon);

	//Number of non-linear inequality constraints
	
	if (scilab_isDouble(env, in[8]) == 0 || scilab_isScalar(env, in[8]) == 0)
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: An int expected.\n", fname, 9);
    	return 1;
	}
	
	scilab_getDouble(env, in[8], &nonlinIneqCon);

	//x0(starting point) matrix from scilab

	if (scilab_isDouble(env, in[13]) == 0 || scilab_isMatrix2d(env, in[13]) == 0)
	{
		Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 14);
		return 1;
	}
	
	scilab_getDoubleArray(env, in[13], &x0ptr);
	size1 = scilab_getDim2d(env, in[13], &x0_rows, &x0_cols);

    //Getting options

	if (scilab_isList(env, in[14]) == 0)
    {
        Scierror(999, "%s: Wrong type for input argument #%d: A list expected.\n", fname, 15);
        return 1;
    }
	scilabVar temp1 = scilab_getListItem( env, in[14], 1);
	scilabVar temp2 = scilab_getListItem( env, in[14], 3);

	double nIters = 0,cpu_Time =0;
	scilab_getDouble(env, temp1, &nIters);
	scilab_getDouble(env, temp2, &cpu_Time);

	int maxIters = (int)nIters;
	int cpu_time = (int)cpu_Time;

	


	//Number of variables and constraints
	nVars = x0_cols;
	nCons = A_rows + Aeq_rows + nonlinCon;


        
        // Starting Ipopt

	SmartPtr<minconNLP> Prob = new minconNLP(env, in, nVars, nCons, x0ptr, Aptr, bptr, Aeqptr, beqptr, A_rows, A_cols, b_rows, b_cols, Aeq_rows, Aeq_cols, beq_rows, beq_cols, lbptr, ubptr, nonlinCon, nonlinIneqCon);
	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

	////////// Managing the parameters //////////

	app->Options()->SetNumericValue("tol", 1e-6);
	app->Options()->SetIntegerValue("max_iter", maxIters);
	app->Options()->SetNumericValue("max_cpu_time", cpu_time);
	//app->Options()->SetStringValue("hessian_approximation", "limited-memory");

	///////// Initialize the IpoptApplication and process the options /////////
	ApplicationReturnStatus status;
	status = app->Initialize();
	if (status != Solve_Succeeded) 
	{
	  	sciprint("\n*** Error during initialization!\n");
   		return (int) status;
 	 }

	 // Ask Ipopt to solve the problem
	 status = app->OptimizeTNLP((SmartPtr<TNLP>&)Prob);
	 //Get the solve statistics
	 cpuTime = app->Statistics()->TotalCPUTime();
	 app->Statistics()->NumberOfEvaluations(int_fobj_eval, int_constr_eval, int_fobj_grad_eval, int_constr_jac_eval, int_hess_eval);
	 app->Statistics()->Infeasibilities(dual_inf, constr_viol, complementarity, kkt_error);
	 rstatus = Prob->returnStatus();
	 fobj_eval=(double)int_fobj_eval;

	////////// Manage the output argument //////////

	fX = Prob->getX();
	fGrad = Prob->getGrad();
	fHess = Prob->getHess();
	fLambda = Prob->getLambda();
	fZl = Prob->getZl();
	fZu = Prob->getZu();
	ObjVal = Prob->getObjVal();
	iteration = (double)app->Statistics()->IterationCount();




	out[0] = scilab_createDoubleMatrix2d(env, 1, nVars, 0);
	scilab_setDoubleArray(env, out[0], fX);

	out[1] = scilab_createDouble(env, ObjVal);

	out[2] = scilab_createDouble(env, rstatus);
	out[3] = scilab_createDouble(env, iteration);
	out[4] = scilab_createDouble(env, cpuTime);
	out[5] = scilab_createDouble(env, fobj_eval);
	
	out[6] = scilab_createDouble(env, dual_inf);

	out[7] = scilab_createDoubleMatrix2d(env, 1, nCons, 0);
	scilab_setDoubleArray(env, out[7], fLambda);
	
	out[8] = scilab_createDoubleMatrix2d(env, 1, nVars, 0);
	scilab_setDoubleArray(env, out[8], fZl);

	out[9] = scilab_createDoubleMatrix2d(env, 1, nVars, 0);
	scilab_setDoubleArray(env, out[9], fZu);

	out[10] = scilab_createDoubleMatrix2d(env, 1, nVars, 0);
	scilab_setDoubleArray(env, out[10], fGrad);

	out[11] = scilab_createDoubleMatrix2d(env, 1, nVars*nVars, 0);
	scilab_setDoubleArray(env, out[11], fHess);

	
	// As the SmartPtrs go out of scope, the reference count
	// will be decremented and the objects will automatically
	// be deleted.*/
	
	return 0;
}
}
