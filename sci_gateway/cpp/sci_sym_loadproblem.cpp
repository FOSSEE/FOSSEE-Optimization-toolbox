
// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Keyur Joshi
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#include "symphony.h"


extern sym_environment* global_sym_env; //defined in globals.cpp

extern "C" {
#include "api_scilab.h"
#include "api_stack_sparse.h"
#include "Scierror.h"
#include "sciprint.h"
#include "BOOL.h"
#include <localization.h>
#include <string.h>

//error management variables

static int iRet;
const char fname[] = "sym_loadProblemBasic";
//data declarations
static int *varAddress=NULL,numVars,numConstr,*conMatrixColStart=NULL,*conMatrixRowIndex=NULL,*isIntVarBool=NULL,colIter,rowIter,inputMatrixCols,inputMatrixRows;
static double inputDouble,objSense,*objective=NULL,*lowerBounds=NULL,*upperBounds=NULL,*conLower=NULL,*conUpper=NULL,*conRange=NULL,*conRHS=NULL,*conMatrix=NULL;
static char *conType=NULL,*isIntVar=NULL;

//delete all allocd arrays before exit, and return output argument



//both basic and advanced loader use this code
static int commonCodePart1(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{
	
	//ensure that environment is active
	if(global_sym_env==NULL)
	{
		sciprint("Error: Symphony environment not initialized. Please run 'sym_open()' first.\n");
		return 1;
	}

	
	//code to check arguments and get them
	if (nin !=10)  //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname, 10);
        	return STATUS_ERROR; 
	}
	
	if (nout !=1) //Checking the output arguments

	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname, 1);
		return 1;
	}


	
	//get input 1: number of variables
	if (scilab_isInt32(env, in[0]) == 0 || scilab_isScalar(env, in[0]) == 0)
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: An int expected.\n", fname, 1);
    	return 1;
	}

	scilab_getInteger32(env, in[0], &numVars);
	
	//get input 2: number of constraints
	if (scilab_isInt32(env, in[1]) == 0 || scilab_isScalar(env, in[1]) == 0)
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: An int expected.\n", fname, 2);
    	return 1;
	}

	scilab_getInteger32(env, in[1], &numConstr);
	
	//allocate and prepare some arrays
	isIntVar=new char[numVars]; //is the variable constrained to be an integer?
	conType=new char[numConstr]; //char representing constraint type
	conRange=new double[numConstr]; //range of each constraint
	conRHS=new double[numConstr]; //RHS to be given to Symphony

	return 0;
}

//both basic and advanced loader use this code
static int commonCodePart2(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{
	//get input 3: lower bounds of variables
	
	if (scilab_isDouble(env, in[2]) == 0 || scilab_isMatrix2d(env, in[2]) == 0)
	{

		Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 3);
		return 1;
	}	
	
	scilab_getDoubleArray(env, in[2], &lowerBounds);
	
	//get input 4: upper bounds of variables

	if (scilab_isDouble(env, in[3]) == 0 || scilab_isMatrix2d(env, in[3]) == 0)
	{

		Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 4);
		return 1;
	}	
	
	scilab_getDoubleArray(env, in[3], &upperBounds);
	
	//get input 5: coefficients of variables in objective function to be minimized

	if (scilab_isDouble(env, in[4]) == 0 || scilab_isMatrix2d(env, in[4]) == 0)
	{

		Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 5);
		return 1;
	}	
	
	scilab_getDoubleArray(env, in[4], &objective);
	
	//get input 6: array that specifies wether a variable is constrained to be an integer
	
	if (scilab_isBoolean(env, in[5]) == 0 || scilab_isMatrix2d(env, in[5]) == 0)
	{

		Scierror(999, "%s: Wrong type for input argument #%d: A Boolean matrix expected.\n", fname, 6);
		return 1;
	}	
	
	scilab_getBooleanArray(env, in[5], &isIntVarBool);


	for(colIter=0;colIter<numVars;colIter++)
	{
		if(isIntVarBool[colIter])
			isIntVar[colIter]=TRUE;
		else
			isIntVar[colIter]=FALSE;
	}
	
	//get input 7: wether to minimize or maximize objective

	if (scilab_isDouble(env, in[6]) == 0 || scilab_isScalar(env, in[6]) == 0)
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: An int expected.\n", fname, 7);
    	return 1;
	}

	scilab_getDouble(env, in[6], &objSense);



	if(objSense!=-1 && objSense!=1)
	{
		Scierror(999, "Wrong type for input argument #7: Either 1 (sym_minimize) or -1 (sym_maximize) is expected.\n");
		return 1;
	}
	iRet=sym_set_obj_sense(global_sym_env,objSense);
	if(iRet==FUNCTION_TERMINATED_ABNORMALLY)
	{
		Scierror(999, "An error occured.\n");
		return 1;
	}
	
	//get input 9: constraint lower bound
	if(!(numConstr == 0))
	{	

		if (scilab_isDouble(env, in[8]) == 0 || scilab_isMatrix2d(env, in[8]) == 0)
		{

			Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 9);
			return 1;
		}	
		
		scilab_getDoubleArray(env, in[8], &conLower);
	}	
	//get input 10: constraint upper bound
	if(!(numConstr == 0))
	{
		
		if (scilab_isDouble(env, in[9]) == 0 || scilab_isMatrix2d(env, in[9]) == 0)
		{

			Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 10);
			return 1;
		}	
		
		scilab_getDoubleArray(env, in[9], &conUpper);
	}	
	//deduce type of constraint
	for(rowIter=0;rowIter<numConstr;rowIter++)
	{
		if(conLower[rowIter]>conUpper[rowIter])
		{
			Scierror(999, "Error: the lower bound of constraint %d is more than its upper bound.\n",rowIter);
			return 1;
		}
		//#ifdef _MSC_VER
  	      //  double INFINITY = sym_get_infinity();
		//#endif
		if(conLower[rowIter]==(-INFINITY) && conUpper[rowIter]==INFINITY){
			conType[rowIter]='N';
			conRange[rowIter]=0;
			conRHS[rowIter]=0;
		}else if(conLower[rowIter]==(-INFINITY)){
			conType[rowIter]='L';
			conRange[rowIter]=0;
			conRHS[rowIter]=conUpper[rowIter];
		}else if(conUpper[rowIter]==INFINITY){
			conType[rowIter]='G';
			conRange[rowIter]=0;
			conRHS[rowIter]=conLower[rowIter];
		}else if(conUpper[rowIter]==conLower[rowIter]){
			conType[rowIter]='E';
			conRange[rowIter]=0;
			conRHS[rowIter]=conLower[rowIter];
		}else{
			conType[rowIter]='R';
			conRange[rowIter]=conUpper[rowIter]-conLower[rowIter];
			conRHS[rowIter]=conUpper[rowIter];
		}
	}
	
	/*
	//for debug: show all data
	sciprint("Vars: %d Constr: %d ObjType: %lf\n",numVars,numConstr,objSense);
	for(colIter=0;colIter<numVars;colIter++)
		sciprint("Var %d: upper: %lf lower: %lf isInt: %d ObjCoeff: %lf\n",colIter,lowerBounds[colIter],upperBounds[colIter],isIntVar[colIter],objective[colIter]);
	for(rowIter=0;rowIter<numConstr;rowIter++)
		sciprint("Constr %d: type: %c lower: %lf upper: %lf range: %lf\n",rowIter,conType[rowIter],conLower[rowIter],conRange[rowIter]);
	*/
	
	//call problem loader
	sym_explicit_load_problem(global_sym_env,numVars,numConstr,conMatrixColStart,conMatrixRowIndex,conMatrix,lowerBounds,upperBounds,isIntVar,objective,NULL,conType,conRHS,conRange,TRUE);
	// sciprint("Problem loaded into environment.\n");
	
	//code to give output

	
	return 0;
}

//basic problem loader, expects normal matrix. Not suitable for large problems


int sci_sym_loadProblemBasic(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{
	
	if(commonCodePart1(env, nin, in, nopt, opt, nout, out))
		return 1;
	
	if(!(numConstr == 0))
	{
	   //get input 8: matrix of constraint equation coefficients

	   
		if (scilab_isDouble(env, in[7]) == 0 || scilab_isMatrix2d(env, in[7]) == 0)
		{
			Scierror(999, "%s: Wrong type for input argument #%d: A double matrix expected.\n", fname, 8);
			return 1;
		}	

		scilab_getDoubleArray(env, in[7], &conMatrix);
	}	 
	

	conMatrixColStart=new int[numVars+1]; //start of each column of constraint matrix, used internally
	conMatrixRowIndex=new int[numVars*numConstr]; //index of column elements in each column, used internally
	for(colIter=0;colIter<numVars;colIter++)
	{
		conMatrixColStart[colIter]=colIter*numConstr;
		for(rowIter=0;rowIter<numConstr;rowIter++)
		{
			conMatrixRowIndex[colIter*numConstr+rowIter]=rowIter;
		}
	}
	conMatrixColStart[numVars]=numVars*numConstr;
	
	if(commonCodePart2(env, nin, in, nopt, opt, nout, out))
		return 1;
	return 0;
}

}
