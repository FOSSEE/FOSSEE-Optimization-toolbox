// Copyright (C) 2016 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Harpreet Singh, Pranav Deshpande and Akshay Miterani
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#include "QuadTMINLP.hpp"
#include "IpIpoptData.hpp"

extern "C"{
#include <sciprint.h>
}

// Go to http://coin-or.org/Ipopt and http://coin-or.org/Bonmin for the details of the below methods

// Set the type of every variable - CONTINUOUS or INTEGER
bool QuadTMINLP::get_variables_types(Index n, VariableType* var_types)
{
  n = numVars_;
  for(int i=0; i < n; i++)
    var_types[i] = CONTINUOUS;
  for(int i=0 ; i < intconSize_ ; ++i)
  	var_types[(int)(intcon_[i]-1)] = INTEGER;
  return true;
}

// The linearity of the variables - LINEAR or NON_LINEAR
bool QuadTMINLP::get_variables_linearity(Index n, Ipopt::TNLP::LinearityType* var_types)
{
	/*
	n = numVars_;
	for(int i = 0; i < n; i++)
		var_types[i] = Ipopt::TNLP::LINEAR;
  */
  return true;
}

// The linearity of the constraints - LINEAR or NON_LINEAR
bool QuadTMINLP::get_constraints_linearity(Index m, Ipopt::TNLP::LinearityType* const_types)
{
  m = numCons_;
  for(int i = 0; i < m; i++)
    const_types[i] = Ipopt::TNLP::LINEAR;
  return true;
}

// Get MINLP info such as the number of variables,constraints,no.of elements in jacobian and hessian to allocate memory
bool QuadTMINLP::get_nlp_info(Index& n, Index&m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
{
  n=numVars_; // Number of variables
	m=numCons_; // Number of constraints
	nnz_jac_g = n*m; // No. of elements in Jacobian of constraints 
	nnz_h_lag = (n*(n+1))/2; // No. of elements in lower traingle of Hessian of the Lagrangian.
	index_style=TNLP::C_STYLE; // Index style of matrices
	return true;
}

// Get the variables and constraints bound info
bool QuadTMINLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
{
  unsigned int i;
	for(i=0;i<n;i++){
		x_l[i]=varLB_[i];
		x_u[i]=varUB_[i];
		}

	for(i=0;i<m;i++){
		g_l[i]=conLB_[i];
		g_u[i]=conUB_[i];
		}
	return true;
}

// This method sets initial values for all the required vectors. We take 0 by default. 
bool QuadTMINLP::get_starting_point(Index n, bool init_x, Number* x,
                             bool init_z, Number* z_L, Number* z_U,
                             Index m, bool init_lambda,
                             Number* lambda)
{
  if (init_x == true){ //we need to set initial values for vector x
		for (Index var=0;var<n;var++)
			x[var]=varGuess_[var];//initialize with 0 or we can change.
		}

	if (init_z == true){ //we need to provide initial values for vector bound multipliers
		for (Index var=0;var<n;++var){
			z_L[var]=0.0; //initialize with 0 or we can change.
			z_U[var]=0.0;//initialize with 0 or we can change.
			}
		}
	
	if (init_lambda == true){ //we need to provide initial values for lambda values.
		for (Index var=0;var<m;++var){
			lambda[var]=0.0; //initialize with 0 or we can change.
			}
		}

	return true;
}

// Evaluate the objective function at a point
bool QuadTMINLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  unsigned int i,j;
	obj_value=0;
	for (i=0;i<n;i++){
		for (j=0;j<n;j++){
			obj_value+=0.5*x[i]*x[j]*qMatrix_[n*i+j];		
			}
  	obj_value+=x[i]*lMatrix_[i];
		}
	return true;
}

// Get the value of gradient of objective function at vector x.
bool QuadTMINLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  unsigned int i,j;
	for(i=0;i<n;i++)
	{
		grad_f[i]=lMatrix_[i];
		for(j=0;j<n;j++)
			{
				grad_f[i]+=(qMatrix_[n*i+j])*x[j];
			}
	}
	return true;
}

// Get the values of constraints at vector x.
bool QuadTMINLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  unsigned int i,j;
	for(i=0;i<m;i++)
	{
		g[i]=0;
		for(j=0;j<n;j++)
		{
			g[i]+=x[j]*conMatrix_[i+j*m];	
		}
	}
	return true;
}

// The Jacobian Matrix
bool QuadTMINLP::eval_jac_g(Index n, const Number* x, bool new_x,
                     Index m, Index nnz_jac, Index* iRow, Index *jCol,
                     Number* values)
{
	//It asks for the structure of the jacobian.
	if (values==NULL){ //Structure of jacobian (full structure)
		int index=0;
		for (int var=0;var<m;++var)//no. of constraints
			for (int flag=0;flag<n;++flag){//no. of variables
				iRow[index]=var;
				jCol[index]=flag;
				index++;
				}
		}
	//It asks for values
	else { 
		int index=0;
		for (int var=0;var<m;++var)
			for (int flag=0;flag<n;++flag)
				values[index++]=conMatrix_[var+flag*m];
		}
	return true;
}

/*
The structure of the Hessain matrix and the values
*/
bool QuadTMINLP::eval_h(Index n, const Number* x, bool new_x,
                 Number obj_factor, Index m, const Number* lambda,
                 bool new_lambda, Index nele_hess, Index* iRow,
                 Index* jCol, Number* values)
{
	if (values==NULL){
		Index idx=0;
		for (Index row = 0; row < n; row++) {
			for (Index col = 0; col <= row; col++) {
				iRow[idx] = row;
				jCol[idx] = col;
				idx++;
		  		}
			}
		}
	else {
		Index index=0;
		for (Index row=0;row < n;++row){
			for (Index col=0; col <= row; ++col){
				values[index++]=obj_factor*(qMatrix_[n*row+col]);
				}
			}
		}
	return true;
}

void QuadTMINLP::finalize_solution(TMINLP::SolverReturn status, Index n, const Number* x,Number obj_value)
{
	
	finalX_ = (double*)malloc(sizeof(double) * numVars_ * 1);
	for (Index i=0; i<n; i++) 
	{
    		 finalX_[i] = x[i];
	}

	finalObjVal_ = obj_value;
	status_ = status;
	sciprint(" Status %d",status);
}

const double * QuadTMINLP::getX()
{	
	return finalX_;
}

double QuadTMINLP::getObjVal()
{	
	return finalObjVal_;
}

int QuadTMINLP::returnStatus()
{	
	return status_;
}
