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

#include "minconNLP.hpp"

extern "C"
{

#include <api_scilab.h>
#include <wchar.h>
#include <Scierror.h>
#include <BOOL.h>
#include <localization.h>
#include <sciprint.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>


#define LOCAL_DEBUG 0

	using namespace std;
	using namespace Ipopt;

minconNLP::~minconNLP()
{
	if(finalX_) delete[] finalX_;
	if(finalZl_) delete[] finalZl_;
	if(finalZu_) delete[] finalZu_;
	if(finalLambda_) delete[] finalLambda_;
}

//get NLP info such as number of variables,constraints,no.of elements in jacobian and hessian to allocate memory
bool minconNLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, IndexStyleEnum& index_style)
{
	finalGradient_ = (double*)malloc(sizeof(double) * numVars_ * 1);
	finalHessian_ = (double*)malloc(sizeof(double) * numVars_ * numVars_);
	
n=numVars_; // Number of variables
m=numConstr_; // Number of constraints

nnz_jac_g = n*m; // No. of elements in Jacobian of constraints 
nnz_h_lag = n*n; // No. of elements in lower traingle of Hessian of the Lagrangian.

index_style=C_STYLE; // Index style of matrices
return true;
}

//get variable and constraint bound info
bool minconNLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
	unsigned int i;
	
	//assigning bounds for the variables
	for(i=0;i<n;i++)
	{
		x_l[i]=varLB_[i];
		x_u[i]=varUB_[i];
	}
	
	if(m==0)
	{
		g_l=NULL;
		g_u=NULL;
	}

	else
	{
		unsigned int c=0;

		//bounds of non-linear inequality constraints
		for(i=0;i<nonlinIneqCon_;i++)
		{
			g_l[c]=-1.0e17;
			g_u[c]=0;
			c++;
		}
		
		//bounds of non-linear equality constraints
		for(i=0;i<nonlinCon_-nonlinIneqCon_;i++)
		{
			g_l[c]=g_u[c]=0;
			c++;
		}

		//bounds of linear equality constraints
		for(i=0;i<Aeqrows_;i++)
		{
			g_l[c]=g_u[c]=beq_[i];
			c++;
		}

		//bounds of linear inequality constraints
		for(i=0;i<Arows_;i++)
		{
			g_l[c]=-1.0e17;
			g_u[c]=b_[i];
			
			c++;
		}

	}

	return true;
}

// This method sets initial values for required vectors . For now we are assuming 0 to all values. 
bool minconNLP::get_starting_point(Index n, bool init_x, Number* x,bool init_z, Number* z_L, Number* z_U,Index m, bool init_lambda,Number* lambda)
{	
	assert(init_x == true);
	assert(init_z == false);
	assert(init_lambda == false);
	if (init_x == true)
	{ //we need to set initial values for vector x
		for (Index var=0;var<n;var++)
			x[var]=varGuess_[var];
	}

	return true;
}

//get value of objective function at vector x
bool minconNLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{	

	scilabVar* out = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 1);
	#ifdef LOCAL_DEBUG
		printf("Calling eval_f\n");
	#endif	
  	double check;
	double obj=0;
	
	const Number *xNew=x;
	scilabVar* funcIn = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 1);
	funcIn[0] = scilab_createDoubleMatrix2d(env_, 1, numVars_, 0);
	scilab_setDoubleArray(env_, funcIn[0], x);

	scilab_call(env_, L"f", 1, funcIn, 2, out);
	

	
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
bool minconNLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
	scilabVar* out = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) );
	const Number *xNew=x;
	#ifdef LOCAL_DEBUG
		printf("grad_f obtained\n");
	#endif
	scilabVar* funcIn = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 1);
	funcIn[0] = scilab_createDoubleMatrix2d(env_, 1, numVars_, 0);
	scilab_setDoubleArray(env_, funcIn[0], x);

	scilab_call(env_, L"fGrad1", 1, funcIn, 2, out);


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
			finalGradient_[i]=resg[i];
		}
	}		
	return true;
}

// return the value of the constraints: g(x)
bool minconNLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  	// return the value of the constraints: g(x)
	unsigned int i;
	unsigned int j;

	if(m==0)
		g=NULL;

	else
	{
		unsigned int c=0;

		//value of non-linear constraints
		if(nonlinCon_!=0)
		{
			int* constr=NULL;  
			const Number *xNew=x;
			double check;


			scilabVar* out = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) );

			#ifdef LOCAL_DEBUG
				printf("grad_f obtained\n");

			#endif
			scilabVar* funcIn = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 1);
			funcIn[0] = scilab_createDoubleMatrix2d(env_, 1, numVars_, 0);
			scilab_setDoubleArray(env_, funcIn[0], x);

			scilab_call(env_, L"addnlc1", 1, funcIn, 2, out);

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

				for(i=0;i<nonlinCon_;i++)
				{
					g[c]=resc[i];
					c++;
				}
			}
		}

		//value of linear equality constraints
		for(i=0;i<Aeqrows_;i++)
		{
			g[c]=0;
			for(j=0;j<Aeqcols_;j++)
				g[c] += Aeq_[j*Aeqrows_+i]*x[j];
			c++;			
		}

		//value of linear inequality constraints
		for(i=0;i<Arows_;i++)
		{
			g[c]=0;
			for(j=0;j<Acols_;j++)
				g[c] += A_[j*Arows_+i]*x[j];
			c++;
		}

	}

	return true;
}

// return the structure or values of the jacobian
bool minconNLP::eval_jac_g(Index n, const Number* x, bool new_x,Index m, Index nele_jac, Index* iRow, Index *jCol,Number* values)
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
    		for(int i=0;i<m;i++)
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
    			values=NULL;

    		else
    		{
    			unsigned int i,j,c=0;
    			double check;
			//jacobian of non-linear constraints
    			if(nonlinCon_!=0)
    			{

					scilabVar* out = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) );
  					double check = 0;
    				
    				const Number *xNew=x;
					scilabVar* funcIn = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 1);
					funcIn[0] = scilab_createDoubleMatrix2d(env_, 1, numVars_, 0);
					scilab_setDoubleArray(env_, funcIn[0], x);
    				
					scilab_call(env_, L"addcGrad1", 1, funcIn, 2, out);

    				
    				double* resj;  
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
					
						scilab_getDoubleArray(env_, out[0], &resj);
    					
    					for(i=0;i<nonlinCon_;i++)
    						for(j=0;j<n;j++)
    						{
    							values[c] = resj[j*(int)nonlinCon_+i];
    							c++;
    						}
    					}			
    				}

			//jacobian of linear equality constraints
    				for(i=0;i<Aeqrows_;i++)
    				{
    					for(j=0;j<Aeqcols_;j++)
    					{
    						values[c] = Aeq_[j*Aeqrows_+i];
    						c++;
    					}
    				}

			//jacobian of linear inequality constraints
    				for(i=0;i<Arows_;i++)
    				{
    					for(j=0;j<Acols_;j++)
    					{
    						values[c] = A_[j*Arows_+i];
    						c++;
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

bool minconNLP::eval_h(Index n, const Number* x, bool new_x,Number obj_factor, Index m, const Number* lambda,bool new_lambda, Index nele_hess, Index* iRow,Index* jCol, Number* values)
{

	scilabVar* out = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 1);
    #ifdef LOCAL_DEBUG
        printf("eval_h obtained\n");
    printf("%d %p %d %p %d %p",n,x,m,lambda,nele_hess,values);
    #endif
	if (values==NULL)
	{
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
		double check;


		const Number *xNew=x;
		#ifdef LOCAL_DEBUG
			printf("in the gradhess block\n");
		#endif	

		const Number *lambdaNew=lambda;
    	double objfac=obj_factor;

		scilabVar* funcIn = (scilabVar*)malloc(sizeof(scilabVar) * (numVars_) * 1);
		funcIn[0] = scilab_createDoubleMatrix2d(env_, 1, numVars_, 0);
		scilab_setDoubleArray(env_, funcIn[0], x);
        funcIn[1] = scilab_createDouble(env_, objfac);
		funcIn[2] = scilab_createDoubleMatrix2d(env_, 1, numConstr_, 0);
		scilab_setDoubleArray(env_, funcIn[2], lambdaNew);
        
        scilab_call(env_, L"lHess1", 3, funcIn, 2, out);
        
		double* resCh = (double*)malloc(sizeof(double)*numVars_*numVars_);
        
		if (scilab_isDouble(env_, out[1]) == 0 || scilab_isScalar(env_, out[1]) == 0)
		{
			Scierror(999, "Wrong type for input argument #%d: An int expected.\n", 2);
			return 1;
		}
		

		scilab_getDouble(env_, out[1], &check);
4
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
	
			scilab_getDoubleArray(env_, out[0], &resCh);
			
			Index index=0;
			for (Index row=0;row < numVars_ ;++row)
			{
				for (Index col=0; col < numVars_; ++col)
				{
					values[index++]=resCh[numVars_*row+col];
					printf("hessian element [%d, %d] = %f\n",row, col, values[numVars_*row+col]);
					
				}
			}
		}

		Index index=0;
		for (Index row=0;row < numVars_ ;++row)
		{
			for (Index col=0; col <= row; ++col)	
			{
				finalHessian_[n*row+col]=values[index++];
			}
		}

		index=0;
		for (Index col=0;col < numVars_ ;++col)
		{
			for (Index row=0; row <= col; ++row)	
			{
				finalHessian_[n*row+col]=values[index++];
			}
		}
		
	}
    #ifdef LOCAL_DEBUG
        printf("Reached end of eval_h\n");
    #endif
	return true;
}

//returning the results
    	void minconNLP::finalize_solution(SolverReturn status,Index n, const Number* x, const Number* z_L, const Number* z_U,Index m, const Number* g, const Number* lambda, Number obj_value,const IpoptData* ip_data,IpoptCalculatedQuantities* ip_cq)
    	{
            #ifdef LOCAL_DEBUG
                printf("Finalize start\n");
            #endif
            finalX_ = new double[n];
    		for (Index i=0; i<n; i++) 
    		{
    			finalX_[i] = x[i];
    		}
    		
    		finalZl_ = new double[n];
    		for (Index i=0; i<n; i++) 
    		{
    			finalZl_[i] = z_L[i];
    		}

    		finalZu_ = new double[n];
    		for (Index i=0; i<n; i++) 
    		{
    			finalZu_[i] = z_U[i];
    		}

    		finalLambda_ = new double[m];
    		for (Index i=0; i<m; i++) 
    		{
    			finalLambda_[i] = lambda[i];
    		}

    		finalObjVal_ = obj_value;
    		status_ = status;
            #ifdef LOCAL_DEBUG
                printf("Finalize end\n");
            #endif
    	}


    	const double * minconNLP::getX()
    	{	
    		return finalX_;
    	}

    	const double * minconNLP::getGrad()
    	{	
    		return finalGradient_;
    	}

    	const double * minconNLP::getHess()
    	{	
    		return finalHessian_;
    	}

    	const double * minconNLP::getZl()
    	{	
    		return finalZl_;
    	}

    	const double * minconNLP::getZu()
    	{	
    		return finalZu_;
    	}

    	const double * minconNLP::getLambda()
    	{	
    		return finalLambda_;
    	}

    	double minconNLP::getObjVal()
    	{	
    		return finalObjVal_;
    	}

    	int minconNLP::returnStatus()
    	{	
    		return status_;
    	}

    }

