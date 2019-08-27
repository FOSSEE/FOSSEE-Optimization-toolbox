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
#include "Scierror.h"
#include "sciprint.h"
#include "BOOL.h"
#include <localization.h>
const char fname[] = "sym_getVarSoln";
int sci_sym_getVarSoln(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{
	
	//error management variable
	int iRet;
	
	//data declarations
	int numVars;
	double *solution;
	
	//ensure that environment is active
	if(global_sym_env==NULL){
		sciprint("Error: Symphony environment not initialized. Please run 'sym_open()' first.\n");
		return 1;
	}
	
	//code to check arguments and get them
	if (nin !=0) //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname, 0);
        	return STATUS_ERROR; 
	}
	
	if (nout !=1) //Checking the output arguments
	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname, 1);
		return 1;
	}
	
	//code to process input
	iRet=sym_get_num_cols(global_sym_env,&numVars);
	if(iRet==FUNCTION_TERMINATED_ABNORMALLY){
		Scierror(999, "An error occured. Has the problem been loaded and solved?\n");
		return 1;
	}
	solution=new double[numVars];
	iRet=sym_get_col_solution(global_sym_env,solution);
	if(iRet==FUNCTION_TERMINATED_ABNORMALLY){
		//Scierror(999, "An error occured. Has the problem been solved? Is the problem feasible?\n");
		delete[] solution;
		return 1;
	}
	
	//code to give output
	out[0] = scilab_createDoubleMatrix2d(env, 1, numVars, 0);
	scilab_setDoubleArray(env, out[0], solution);
	
	delete[] solution;
	
	return 0;
}

const char fname2[] = "sym_getObjVal";
int sci_sym_getObjVal(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{
	
	//error management variable
	SciErr sciErr;
	int iRet;
	
	//data declarations
	double solution;
	
	//ensure that environment is active
	if(global_sym_env==NULL)
	{
		sciprint("Error: Symphony environment not initialized. Please run 'sym_open()' first.\n");
		return 1;
	}
	
	//code to check arguments and get them
	if (nin !=0) //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname2, 0);
        	return STATUS_ERROR; 
	}
	
	if (nout !=1) //Checking the output arguments
	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname2, 1);
		return 1;
	}
	
	//code to process input
	iRet=sym_get_obj_val(global_sym_env,&solution);
	if(iRet==FUNCTION_TERMINATED_ABNORMALLY){
		Scierror(999, "An error occured. Has the problem been loaded and solved? Is the problem feasible?\n");
		return 1;
	}
	
	//code to give output

	out[0] = scilab_createDouble(env, solution);
	
	return 0;
}

}
