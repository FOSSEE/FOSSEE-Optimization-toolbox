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

#include "minconTMINLP.hpp"
#include "sci_iofunc.hpp"

extern "C"
{

#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <localization.h>
#include <wchar.h>
#include <sciprint.h>
#include <string.h>
#include <assert.h>
}

using namespace Ipopt;
using namespace Bonmin;



minconTMINLP::~minconTMINLP()
{
	if(finalX_) delete[] finalX_;
}



// Set the type of every variable - CONTINUOUS or INTEGER
bool minconTMINLP::get_variables_types(Index n, VariableType* var_types)
{
	#ifdef LOCAL_DEBUG
  		printf("Code is in get_variables_types\n");
	#endif
  n = numVars_;
  for(int i=0; i < n; i++)
    var_types[i] = CONTINUOUS;
  for(int i=0 ; i < intconSize_ ; ++i)
  	var_types[(int)(intcon_[i]-1)] = INTEGER;
  return true;
}

// The linearity of the variables - LINEAR or NON_LINEAR
bool minconTMINLP::get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types)
{  
	#ifdef LOCAL_DEBUG
  		printf("Code is in get_variables_linearity\n");
	#endif
	for(int i=0;i<n;i++)
	{
		  var_types[i] = Ipopt::TNLP::NON_LINEAR;
	}
	return true; }

// The linearity of the constraints - LINEAR or NON_LINEAR
bool minconTMINLP::get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types)
{	

	#ifdef LOCAL_DEBUG
  		printf("Code is in get_constraints_linearity\n");
	#endif
	for(int i=0;i<numLC_;i++)
	{
		const_types[i] = Ipopt::TNLP::LINEAR;
	}

	for(int i=numLC_;i<m;i++)
	{
		const_types[i] = Ipopt::TNLP::NON_LINEAR;
	}
	  return true;}

//get NLP info such as number of variables,constraints,no.of elements in jacobian and hessian to allocate memory
bool minconTMINLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
{
	#ifdef LOCAL_DEBUG
  		printf("Code is in get_nlp_info\n");
	#endif
	n=numVars_; // Number of variables
	m=numCons_; // Number of constraints
	nnz_jac_g = n*m; // No. of elements in Jacobian of constraints 
	nnz_h_lag = n*n; // No. of elements in Hessian of the Lagrangian.
	index_style=TNLP::C_STYLE; // Index style of matrices
	return true;
}

//get variable and constraint bound info
bool minconTMINLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
	#ifdef LOCAL_DEBUG
  		printf("Code is in get_bounds_info\n");
	#endif
	unsigned int i;
	for(i=0;i<n;i++)
	{
		x_l[i]=lb_[i];
		x_u[i]=ub_[i];
	}
	for(i=0;i<m;i++)
	{
		g_l[i]=conLb_[i];
        g_u[i]=conUb_[i];
	}
	return true;
}

// This method sets initial values for required vectors . For now we are assuming 0 to all values. 
bool minconTMINLP::get_starting_point(Index n, bool init_x, Number* x,bool init_z, Number* z_L, Number* z_U,Index m, bool init_lambda,Number* lambda)
{
 	assert(init_x == true);
  	assert(init_z == false);
  	assert(init_lambda == false);
	if (init_x == true)
	{ //we need to set initial values for vector x
		for (Index var=0;var<n;var++)
			{x[var]=x0_[var];}//initialize with 0.
	}
	return true;
}

//get value of objective function at vector x
bool minconTMINLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{	
	scilabVar* out = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 1);
	#if LOCAL_DEBUG
		printf("Calling eval_f\n");
	#endif	
  	double check;
	double obj=0;
	
	const Number *xNew=x;
	scilabVar* funcIn = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 1);
	funcIn[0] = scilab_createDoubleMatrix2d(env_, 1, numVars_, 0);
	scilab_setDoubleArray(env_, funcIn[0], x);

	scilab_call(env_, L"_f", 1, funcIn, 2, out);
	

	
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
		if (scilab_isDouble(env_, out[0]) == 0 || scilab_isScalar(env_, out[0]) == 0)
		{
			sciprint("No obj value\n");
			return 1;
		}

		scilab_getDouble(env_, out[0], &obj);
  		obj_value=obj;    
		
		return true;
	}
}

//get value of gradient of objective function at vector x.
bool minconTMINLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
	scilabVar* out = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) );
	const Number *xNew=x;
	#if LOCAL_DEBUG
		printf("grad_f obtained\n");
	#endif
	scilabVar* funcIn = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 1);
	funcIn[0] = scilab_createDoubleMatrix2d(env_, 1, numVars_, 0);
	scilab_setDoubleArray(env_, funcIn[0], x);

	scilab_call(env_, L"_gradf", 1, funcIn, 2, out);


	double* resg;
	double check;
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
	
		scilab_getDoubleArray(env_, out[0], &resg);


		Index i;
		for(i=0;i<numVars_;i++)
		{
			grad_f[i]=resg[i];

		}
	}		
	return true;
}

// return the value of the constraints: g(x)
bool minconTMINLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
	
 	// return the value of the constraints: g(x)
  	if(m==0)
  	{
  		g=NULL;
  	}  		
	else
	{
	  	unsigned int c=0;

		//value of non-linear constraints
		
		int* constr=NULL;  
		const Number *xNew=x;
		double check;


		scilabVar* out = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) );
		#if LOCAL_DEBUG
			printf("g obtained\n");
		#endif
		scilabVar* funcIn = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 1);
		funcIn[0] = scilab_createDoubleMatrix2d(env_, 1, numVars_, 0);
		scilab_setDoubleArray(env_, funcIn[0], x);

		scilab_call(env_, L"_addnlc", 1, funcIn, 2, out);

		double* resc;  
                     
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
	
			scilab_getDoubleArray(env_, out[0], &resc);



			for(int i=0;i<m;i++)
			{
				g[c]=resc[i];
				//printf("g element [%d] = %f\n",m, resc[i]);
				c++;
			}
		}
		
	}

  	return true;
}

// return the structure or values of the jacobian
bool minconTMINLP::eval_jac_g(Index n, const Number* x, bool new_x,Index m, Index nele_jac, Index* iRow, Index *jCol,Number* values)
{
 	if (values == NULL) 
 	{
    	if(m==0)// return the structure of the jacobian of the constraints
    	{
			iRow=NULL; 
    		jCol=NULL;
  		}
		else
		{
			unsigned int i,j,idx=0;
			for(i=0;i<m;i++)
				for(j=0;j<n;j++)
				{
					iRow[idx]=i;
					jCol[idx]=j;
					idx++;
				}
		}
	}
	else 
	{
		if(m==0)
		{
			values=NULL;
		}
		else
		{
			double* resj;


			scilabVar* out = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 2);
		  	double check = 0;

			scilabVar* funcIn = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 2);
			funcIn[0] = scilab_createDoubleMatrix2d(env_, 1, numVars_, 0);
			scilab_setDoubleArray(env_, funcIn[0], x);

			#if LOCAL_DEBUG
				printf("eval_jac_g scilab_setDoubleArray\n");
			#endif	

			scilab_call(env_, L"_gradnlc", 1, funcIn, 2, out);

			if (scilab_isDouble(env_, out[1]) == 0 || scilab_isScalar(env_, out[1]) == 0)
			{
				Scierror(999, "Wrong type for input argument #%d: An int expected.\n", 2);
				return 1;
			}
			
			#if LOCAL_DEBUG
				printf("eval_jac_g check\n");
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
			
				scilab_getDoubleArray(env_, out[0], &resj);

			
				int c = 0;
				for(int i=0;i<m;i++)
				{	for(int j=0;j<n;j++)
					{
						values[c] = resj[j*(int)m+i];
						//printf("jac_g element [%d, %d] = %f\n",m, n, resj[j*(int)m+i]);
						c++;
					}			
				}
			}
		}	
	}
  	return true;
}

/*
 * Return either the sparsity structure of the Hessian of the Lagrangian, 
 * or the values of the Hessian of the Lagrangian  for the given values for
 * x,lambda,obj_factor.
*/

bool minconTMINLP::eval_h(Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,bool new_lambda, Index nele_hess, Index* iRow,Index* jCol, Number* values)
{
	#ifdef DEBUG
  		printf("Code is in eval_h\n");
	#endif	
	if (values==NULL)
	{

		printf("values=NULL\n");
		Index idx=0;
		for (Index row = 0; row < numVars_; row++) 
		{
			for (Index col = 0; col < numVars_; col++)
			{
				iRow[idx] = row;
				jCol[idx] = col;
				idx++;
		  	}
		}
	}	
	else 
	{	

		for (int i = 0; i < n; i++)
		{
			//printf("x[%d] before enter getHess = %f\n",i,x[i]);
		}
	  	
		double check;
		scilabVar* out = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 2);

		scilabVar* funcIn = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_ * numCons_ + 1));

		funcIn[0] = scilab_createDoubleMatrix2d(env_, 1, n, 0);
		scilab_setDoubleArray(env_, funcIn[0], x);

		funcIn[1] =  scilab_createDouble(env_, obj_factor);

		funcIn[2] = scilab_createDoubleMatrix2d(env_, 1, numCons_, 0);
		scilab_setDoubleArray(env_, funcIn[2],lambda);	


	  	scilab_call(env_, L"_gradhess", 3, funcIn, 2, out);

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
			
			int r1; int c1;
			int size1 = scilab_getDim2d(env_, out[0], &r1, &c1);
			Number *resCh =(double*)malloc(sizeof(double)*r1*c1);
		
			scilab_getDoubleArray(env_, out[0], &resCh);

		
		
			Index index=0;
			for (Index row=0;row < numVars_ ;++row)
			{
				for (Index col=0; col < numVars_; ++col)
				{
					values[index++]=resCh[numVars_*row+col];
					//printf("hessian element [%d, %d] = %f\n",row, col, resCh[numVars_*row+col]);
				}
			}
		}
	}
	return true;
}

void minconTMINLP::finalize_solution(SolverReturn status,Index n, const Number* x, Number obj_value)
{
	#ifdef LOCAL_DEBUG
  		printf("Code is in finalize_solution\n");
  		printf("%d",status);
	#endif	
	finalObjVal_ = obj_value;
	status_ = status;
	if(status==0 ||status== 3)
	{
		finalX_ =  new double[n];
		for (Index i=0; i<numVars_; i++) 
		{
	    		 finalX_[i] = x[i];
		}
	}
}

const double * minconTMINLP::getX()
{	
	return finalX_;
}

double minconTMINLP::getObjVal()
{	
	return finalObjVal_;
}

int minconTMINLP::returnStatus()
{	
	return status_;
}
