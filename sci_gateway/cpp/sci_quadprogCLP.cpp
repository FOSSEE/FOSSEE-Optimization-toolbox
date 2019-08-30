// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Adarsh Shah
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"
#include "CoinPackedMatrix.hpp"
#include "CoinPackedVector.hpp"
#include "ClpInterior.hpp"
#include "ClpSimplex.hpp"
#include "ClpModel.hpp"
#include "ClpCholeskyBase.hpp"
#include "ClpQuadraticObjective.hpp"
#include <cassert>
#include <iostream>

extern "C"{

#include <api_scilab.h>
#include <Scierror.h>
#include <localization.h>
#include <sciprint.h>
//Solver function
 const char fname[] = "quadprog_CLP";
/* ==================================================================== */
int sci_quadprog_CLP(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out) 
{

	//Quad Matrix H
	double* H = NULL;
	//Linear Matrix C
	double* C = NULL;
	//Constant d
	double d = 0;
	//Linear Inequality Matrix A
	double* A = NULL;
	//Linear Inequality Matrix b
	double* b = NULL;
	//Linear Equality Matrix N
	double* N = NULL;
	//Linear Equality Matrix e	
	double* e = NULL;
    //variables lower and upper bound
    double * lb = NULL;
    double * ub = NULL;
	//No of decision variables n
	int n = 0;
	int size = 0;
    int temp = 0;
	int temp1 = 0;
	//No of equality constraints p
	int p = 0;
	//No of inequality constraints m
	int m = 0;
	//isQuad variable
	int isQuad = 0;
    int i = 0;
    int j = 0;
	
	if (nin != 9) //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname, 9);
        	return STATUS_ERROR; 
	}
    
	if (nout !=6) //Checking the output arguments
	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname, 6);
		return 1;
	}
	////////// Manage the input argument //////////
	
	//Checking and Retriving Quadratic Matrix H 
	if (scilab_isDouble(env, in[0]) == 0 || scilab_isMatrix2d(env, in[0]) == 0)
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: A square matrix expected.\n", fname, 1);
    	return 1;
	}
	
	size = scilab_getDim2d(env,in[0], &n, &temp);
	
	if(n>0){
		isQuad = 1;
	}
	else if(n!=temp)
	{
		Scierror(999, "%s: Wrong type for input argument #%d: A square matrix expected.\n", fname, 1);
    	return 1;	
	}
	
	scilab_getDoubleArray(env, in[0], &H);

	//Checking and Retriving Linear Matrix C 
	if (scilab_isDouble(env, in[1]) == 0)
	{
    		Scierror(999, "%s: Wrong type for input argument #%d: An column vector expected.\n", fname, 2);
    		return 1;
	}
	
	if(isQuad){
		size = scilab_getDim2d(env,in[1], &temp, &temp1);
		if(n!=temp || temp1 != 1)
		{
			Scierror(999, "%s: Wrong type for input argument #%d: dimension of matrix is invalid.\n", fname, 2);
    		return 1;	
		}
	}
	else{
		size = scilab_getDim2d(env,in[1], &n, &temp1);
		if(temp1 != 1)
		{
			Scierror(999, "%s: Wrong type for input argument #%d: dimension of matrix is invalid.\n", fname, 2);
    		return 1;	
		}
	}
	scilab_getDoubleArray(env, in[1], &C);

	//Checking and Retriving Constant d
	if (scilab_isDouble(env, in[2]) == 0 || scilab_isScalar(env, in[2]) == 0)
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: An int expected.\n", fname, 3);
    	return 1;
	}

	scilab_getDouble(env, in[2], &d);

	//Checking and Retriving Linear Inequality Matrix A
	if (scilab_isDouble(env, in[3]) == 0 || scilab_isMatrix2d(env, in[3]) == 0)
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: A Matrix expected.\n", fname, 4);
    	return 1;
	}
	
	size = scilab_getDim2d(env,in[3], &m, &temp);
	if(m!=0){
		if(n!=temp)
		{
			Scierror(999, "%s: Wrong type for input argument #%d: dimension of matrix is invalid.\n", fname, 4);
    		return 1;	
		}
	
		scilab_getDoubleArray(env, in[3], &A);

	//Checking and Retriving Linear Matrix b 
		if (scilab_isDouble(env, in[4]) == 0)
		{
    		Scierror(999, "%s: Wrong type for input argument #%d: An column vector expected.\n", fname, 5);
    		return 1;
		}
	
		size = scilab_getDim2d(env,in[4], &temp, &temp1);
		if(m!=temp || temp1 != 1)
		{
			Scierror(999, "%s: Wrong type for input argument #%d: dimension of matrix is invalid.\n", fname, 5);
    		return 1;	
		}
	
		scilab_getDoubleArray(env, in[4], &b);
	}
	
	
	//Checking and Retriving Linear Equality Matrix N
	if (scilab_isDouble(env, in[5]) == 0 || scilab_isMatrix2d(env, in[5]) == 0)
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: A Matrix expected.\n", fname, 6);
    	return 1;
	}
	
	size = scilab_getDim2d(env,in[5], &p, &temp);
	if(p!=0){
		if(n!=temp)
		{
			Scierror(999, "%s: Wrong type for input argument #%d: dimension of matrix is invalid.\n", fname, 6);
    		return 1;	
		}
	
		scilab_getDoubleArray(env, in[5], &N);

		//Checking and Retriving Linear Equality Matrix e 
		if (scilab_isDouble(env, in[6]) == 0)
		{
    		Scierror(999, "%s: Wrong type for input argument #%d: An column vector expected.\n", fname, 7);
    		return 1;
		}
	
		size = scilab_getDim2d(env,in[6], &temp, &temp1);

	    if(p!=temp || temp1 != 1)
		{
			Scierror(999, "%s: Wrong type for input argument #%d: dimension of matrix is invalid.\n", fname, 7);
 		   	return 1;	
		}
	
		scilab_getDoubleArray(env, in[6], &e);
	}
    //Checking and retriving Lower Bounds Matrix lb
    if (scilab_isDouble(env, in[7]) == 0)
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: An column vector expected.\n", fname, 8);
    	return 1;
	}

    size = scilab_getDim2d(env,in[7], &temp, &temp1);

    if(n!=temp || temp1 != 1)
	{
	    Scierror(999, "%s: Wrong type for input argument #%d: dimension of matrix is invalid.\n", fname, 8);
    	return 1;	
	}
	
	scilab_getDoubleArray(env, in[7], &lb);

    //Checking and retriving Upper Bounds Matrix ub
    if (scilab_isDouble(env, in[8]) == 0)
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: An column vector expected.\n", fname, 9);
    	return 1;
	}

    size = scilab_getDim2d(env,in[8], &temp, &temp1);

    if(n!=temp || temp1 != 1)
	{
	    Scierror(999, "%s: Wrong type for input argument #%d: dimension of matrix is invalid.\n", fname, 9);
    	return 1;	
	}
	
	scilab_getDoubleArray(env, in[8], &ub);

	ClpSimplex model;
	try
		{
			/* code */
	
	
    CoinPackedMatrix * matrix;
    CoinPackedVector * col;
	//Quadratic Objective Function
    if(isQuad){
        matrix = new CoinPackedMatrix(false,0,0);
        matrix->setDimensions(n,0);
        temp = 0;
        for(i = 0 ; i < n ; i++){
            col = new CoinPackedVector();
            for(j = 0 ; j < n ; j++){
                col->insert(j,H[temp++]);
            }
            matrix->appendCol(*col);
        }
    }
	
	//Linear Contraints Function
    CoinPackedMatrix * constraints = new CoinPackedMatrix(false,0,0);
    constraints->setDimensions(m+p,0);
    temp=0;
    temp1=0;
    for(i = 0 ;i < n ; i++){
        col = new CoinPackedVector();
        for(j = 0 ; j < m ; j++){
            col->insert(j,A[temp++]);
        }
        for(j = 0 ; j < p ; j++){
            col->insert(j+m,N[temp1++]); 
        }
        constraints->appendCol(*col);
    }  

	//row lower bounds
    //row upper bounds
    double * row_lb = new double[m+p];
    double * row_ub = new double[m+p];
    for(i = 0 ; i < m ; i++){
        row_lb[i] = -COIN_DBL_MAX;
        row_ub[i] = b[i];
    }

    for(i = 0 ; i < p ; i++){
        row_lb[i+m] = e[i];
        row_ub[i+m] = e[i];
    }
	
	
	

	model.loadProblem(*constraints,lb,ub,C,row_lb,row_ub);	
	if(isQuad){
		model.loadQuadraticObjective(*matrix);
	}
    model.initialSolve();
	
	}
	catch(const std::exception& e)
	{
		std::cerr << e.what() << '\n';
	}
	//Output the solution to Scilab
	//get solution for x
	const double* xValue = NULL;
	xValue = model.getColSolution();
	
	//get objective value
	double objValue = model.getObjValue();	
	
	//get Status value
	double status_ = 0;
	if(model.isProvenOptimal())
			status_=0;
	else if(model.isProvenPrimalInfeasible())
			status_=1;
	else if(model.isProvenDualInfeasible())
			status_=2;
	else if(model.isIterationLimitReached())
			status_=3;
	else if(model.isAbandoned())
			status_=4;
	else if(model.isPrimalObjectiveLimitReached())
			status_=5;
	else if(model.isDualObjectiveLimitReached())
			status_=6;
		
	//get number of iterations
	double iterations  = model.getIterationCount(); 	
	
	//get reduced cost
	const double* Zl = model.getReducedCost();
	
	//get dual vector
	const double* dual = model.getRowPrice();
	
	//Create Output matrices
	out[0] = scilab_createDoubleMatrix2d(env, n, 1, 0);
	out[4] = scilab_createDoubleMatrix2d(env, n, 1, 0);
	out[5] = scilab_createDoubleMatrix2d(env, m+p, 1, 0);

	
	scilab_setDoubleArray(env, out[0], xValue);
	out[1] = scilab_createDouble(env, objValue);
	out[2] = scilab_createDouble(env, status_);
	out[3] = scilab_createDouble(env, iterations);
	scilab_setDoubleArray(env, out[4], Zl);
	scilab_setDoubleArray(env, out[5], dual);

	return 0;	

	}

}


  	
   
