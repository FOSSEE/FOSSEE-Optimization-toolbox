// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// Author: Harpreet Singh, Pranav Deshpande and Akshay Miterani
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

#include "minuncTMINLP.hpp"

extern "C"
{
#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <localization.h>
#include <sciprint.h>
#include <wchar.h>
#include <assert.h>
}

#define LOCAL_DEBUG 1

using namespace Ipopt;
using namespace Bonmin;

minuncTMINLP::~minuncTMINLP()
{
	#if LOCAL_DEBUG
		printf("Calling destructor\n");
	#endif
	if(finalX_) delete[] finalX_;
}

// Set the type of every variable - CONTINUOUS or INTEGER
bool minuncTMINLP::get_variables_types(Index n, VariableType* var_types)
{

  #if LOCAL_DEBUG
		printf("Calling get_variables_types\n");
	#endif	
  n = numVars_;
  for(int i=0; i < n; i++)
    var_types[i] = CONTINUOUS;
  for(int i=0 ; i < intconSize_ ; ++i)
  	var_types[(int)(intcon_[i]-1)] = INTEGER;
  return true;
}

// The linearity of the variables - LINEAR or NON_LINEAR
bool minuncTMINLP::get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types)
{
	/*
	n = numVars_;
	for(int i = 0; i < n; i++)
		var_types[i] = Ipopt::TNLP::LINEAR;
  */
  return true;
}

// The linearity of the constraints - LINEAR or NON_LINEAR
bool minuncTMINLP::get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types)
{
 /* m = numConstr_;
  for(int i = 0; i < m; i++)
    const_types[i] = Ipopt::TNLP::LINEAR;
  */
  return true;
}

//get NLP info such as number of variables,constraints,no.of elements in jacobian and hessian to allocate memory
bool minuncTMINLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
{

	#if LOCAL_DEBUG
		printf("Calling get_nlp_info\n");
	#endif
	n=numVars_; // Number of variables
	m=0; // Number of constraints
	nnz_jac_g = 0; // No. of elements in Jacobian of constraints 
	nnz_h_lag = n*(n+1)/2; // No. of elements in lower traingle of Hessian of the Lagrangian.
	index_style=TNLP::C_STYLE; // Index style of matrices
	return true;
}

//get variable and constraint bound info
bool minuncTMINLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{

	#if LOCAL_DEBUG
		printf("Calling get_bounds_info\n");
	#endif
	unsigned int i;
	for(i=0;i<n;i++)
	{
		x_l[i]=-1.0e19;
		x_u[i]=1.0e19;
	}

        g_l=NULL;
        g_u=NULL;
	return true;
}

// return the value of the constraints: g(x)
bool minuncTMINLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  	// return the value of the constraints: g(x)
  	g=NULL;
  	return true;
}

// return the structure or values of the jacobian
bool minuncTMINLP::eval_jac_g(Index n, const Number* x, bool new_x,Index m, Index nele_jac, Index* iRow, Index *jCol,Number* values)
{

	#if LOCAL_DEBUG
		printf("Calling eval_jac_g\n");
	#endif
 	if (values == NULL) 
 	{
    		// return the structure of the jacobian of the constraints
    		iRow=NULL; 
    		jCol=NULL;
  	}
  	else 
	{
		values=NULL;
  	}

  	return true;
}

//get value of objective function at vector x
bool minuncTMINLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
	double obj=0;

	scilabVar* out = (scilabVar*)malloc(sizeof(double) * (numVars_) * 2);
	#if LOCAL_DEBUG
		printf("Calling eval_f\n");
	#endif	
  	double check;

	scilabVar* funcIn = (scilabVar*)malloc(sizeof(double) * (numVars_) * 2);
	funcIn[0] = scilab_createDoubleMatrix2d(env_, 1, numVars_, 0);
	scilab_setDoubleArray(env_, funcIn[0], x);

	scilab_call(env_, L"_f", 1, funcIn, 2, out);

	if (check==1)
	{
		return true;
	}	
	else
	{ 

		if (scilab_isDouble(env_, out[0]) == 0 || scilab_isScalar(env_, out[0]) == 0)
		{
			sciprint("No obj value\n");
			return 1;
		}

		scilab_getDouble(env_, out[0], &obj);
  		obj_value=obj;  
	}
	#if LOCAL_DEBUG
		printf("Obj value obtained\n");
	#endif
  	return true;
}

//get value of gradient of objective function at vector x.
bool minuncTMINLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
	char name[20]="_gradf";
  	Number *resg;

	#if LOCAL_DEBUG
		printf("eval_grad_f started\n");
	#endif	

	scilabVar* out = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 2);
  	double check = 0;

	scilabVar* funcIn = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 2);
	funcIn[0] = scilab_createDoubleMatrix2d(env_, 1, numVars_, 0);
	scilab_setDoubleArray(env_, funcIn[0], x);

	#if LOCAL_DEBUG
		printf("eval_grad_f scilab_setDoubleArray\n");
	#endif	

	scilab_call(env_, L"_gradf", 1, funcIn, 2, out);


	if (scilab_isDouble(env_, out[1]) == 0 || scilab_isScalar(env_, out[1]) == 0)
	{
    	Scierror(999, "Wrong type for input argument #%d: An int expected.\n", 2);
    	return 1;
	}
	
	#if LOCAL_DEBUG
		printf("eval_grad_f check\n");
	#endif
	scilab_getDouble(env_, out[1], &check);
	

	if (check==1)
	{
		return true;
	}	
	else
	{ 

  		if (scilab_isDouble(env_, out[0]) == 0 || scilab_isMatrix2d(env_, out[0]) == 0)
		{
			Scierror(999, "Wrong type for input argument #%d: An int expected.\n", 2);
			return 1;
		}
	
		scilab_getDoubleArray(env_, out[0], &resg);
  		Index i;
  		for(i=0;i<numVars_;i++)
  		{
			grad_f[i]=resg[i];
  		}
	}

	#if LOCAL_DEBUG
		printf("eval_grad_f finished\n");
	#endif
	return true;
}

// This method sets initial values for required vectors . For now we are assuming 0 to all values. 
bool minuncTMINLP::get_starting_point(Index n, bool init_x, Number* x,bool init_z, Number* z_L, Number* z_U,Index m, bool init_lambda,Number* lambda)
{
 	assert(init_x == true);
  	assert(init_z == false);
  	assert(init_lambda == false);
	if (init_x == true)
	{ //we need to set initial values for vector x
		for (Index var=0;var<n;var++)
			x[var]=varGuess_[var];//initialize with 0 or we can change.
	}

	return true;
}

/*
 * Return either the sparsity structure of the Hessian of the Lagrangian, 
 * or the values of the Hessian of the Lagrangian  for the given values for
 * x,lambda,obj_factor.
*/

bool minuncTMINLP::eval_h(Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,bool new_lambda, Index nele_hess, Index* iRow,Index* jCol, Number* values)
{
	#if LOCAL_DEBUG
		printf("eval_h started\n");
	#endif

	scilabVar* out = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 2);
	double check;
	if (values==NULL)
	{
		Index idx=0;
		for (Index row = 0; row < numVars_; row++) 
		{
			for (Index col = 0; col <= row; col++)
			{
				iRow[idx] = row;
				jCol[idx] = col;
				idx++;
		  	}
		}
	}

	else 
	{			
	  	Number *resh;

		#if LOCAL_DEBUG
			printf("in the gradhess block\n");
		#endif	
		scilabVar* out = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 2);
	  	double check = 0;

		scilabVar* funcIn = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 2);
		funcIn[0] = scilab_createDoubleMatrix2d(env_, 1, numVars_, 0);
		scilab_setDoubleArray(env_, funcIn[0], x);

	#if LOCAL_DEBUG
		printf("eval_h scilab_setDoubleArray\n");
	#endif	
		scilab_call(env_, L"_gradhess", 1, funcIn, 2, out);


		if (scilab_isDouble(env_, out[1]) == 0 || scilab_isScalar(env_, out[1]) == 0)
		{
			Scierror(999, "Wrong type for input argument #%d: An int expected.\n", 2);
			return 1;
		}
		

		scilab_getDouble(env_, out[1], &check);

		if (check==1)
		{
			
			return true;
		}	

		else
		{
			if (scilab_isDouble(env_, out[0]) == 0 || scilab_isMatrix2d(env_, out[0]) == 0)
			{
				Scierror(999, "Wrong type for input argument #%d: An int expected.\n", 2);
				return 1;
			}
		
			scilab_getDoubleArray(env_, out[0], &resh);

			Index index=0;
			for (Index row=0;row < numVars_ ;++row)
			{
				for (Index col=0; col <= row; ++col)
				{
					values[index++]=obj_factor*(resh[numVars_*row+col]);
				}
			}
		}
    }
    return true;
}


void minuncTMINLP::finalize_solution(SolverReturn status,Index n, const Number* x, Number obj_value)
{

	#if LOCAL_DEBUG
		printf("Calling finalize_solution\n");
	#endif	
	finalObjVal_ = obj_value;
	status_ = status;
	if(status==0 ||status== 3)
	{
		finalX_ = new double[n];
		for (Index i=0; i<numVars_; i++) 
		{
	    		 finalX_[i] = x[i];
		}
	}
	#if LOCAL_DEBUG
		printf("Finished finalize_solution\n");
	#endif	
	
}

const double * minuncTMINLP::getX()
{	
	
	return finalX_;
	
}

double minuncTMINLP::getObjVal()
{	

	return finalObjVal_;
	
}

int minuncTMINLP::returnStatus()
{	

	return status_;
	
}
