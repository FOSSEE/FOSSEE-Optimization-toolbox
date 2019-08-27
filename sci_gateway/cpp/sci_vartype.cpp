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

//error management variable
static SciErr sciErr;
static int iRet;

//data declarations
static int *varAddress,varIndex,numVars,retVal;
static double inputDouble;



static int commonCodePart1(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{
	//ensure that environment is active
	if(global_sym_env==NULL){
		sciprint("Error: Symphony environment not initialized. Please run 'sym_open()' first.\n");
		return 1;
	}
	
	//code to check arguments and get them
	if (nin !=1)  //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname, 0);
        	return STATUS_ERROR; 
	}
	
	if (nout !=1) //Checking the output arguments

	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname, 1);
		return 1;
	}
	return 1;
	
	//code to process input
	if(getUIntFromScilab(1,&varIndex))
		return 1;
	iRet=sym_get_num_cols(global_sym_env,&numVars);
	if(iRet==FUNCTION_TERMINATED_ABNORMALLY){
		Scierror(999, "An error occured. Has a problem been loaded?\n");
		return 1;
	}else if(varIndex>=numVars){
		Scierror(999, "An error occured. Variable index must be a number between 0 and %d.\n",numVars-1);
		return 1;
	}
	
	return 0;
}

static int commonCodePart2(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{
	out[0] = scilab_createDouble(env, retVal);
	
	return 0;
}

const char fname[] = "sym_isContinuous";

int sci_sym_isContinuous(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{

	
	if(commonCodePart1(env, nin, in, nopt, opt, nout, out))
		return 1;
	
	iRet=sym_is_continuous(global_sym_env,varIndex,&retVal);
	if(iRet==FUNCTION_TERMINATED_ABNORMALLY){
		Scierror(999, "An error occured. Has a problem been loaded?\n");
		return 1;
	}else{
		if(retVal)
			sciprint("This variable is continuous.\n");
		else
			sciprint("This variable is not continuous.\n");
	}
	
	if(commonCodePart2(env, nin, in, nopt, opt, nout, out, ))
		return 1;
	
	return 0;
}

const char fname2[] = "sym_isBinary";

int sci_sym_isBinary(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{

	
	if(commonCodePart1(env, nin, in, nopt, opt, nout, out))
		return 1;

	iRet=sym_is_binary(global_sym_env,varIndex,&retVal);
	if(iRet==FUNCTION_TERMINATED_ABNORMALLY){
		Scierror(999, "An error occured.\n");
		return 1;
	}else{
		if(retVal)
			sciprint("This variable is constrained to be binary.\n");
		else
			sciprint("This variable is not constrained to be binary.\n");
	}
	
	if(commonCodePart2(env, nin, in, nopt, opt, nout, out))
		return 1;
	
	return 0;
}

const char fname2[] = "sym_isInteger";

int sci_sym_isInteger(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{

	
	char retValc; //for some wierd reason this function unlike the above 2 returns a char
	
	if(commonCodePart1(env, nin, in, nopt, opt, nout, out))
		return 1;
	
	iRet=sym_is_integer(global_sym_env,varIndex,&retValc);
	if(iRet==FUNCTION_TERMINATED_ABNORMALLY){
		Scierror(999, "An error occured.\n");
		return 1;
	}else{
		if(retValc)
			sciprint("This variable is constrained to be an integer.\n");
		else
			sciprint("This variable is not constrained to be an integer.\n");
	}
	retVal=retValc;
	
	if(commonCodePart2(env, nin, in, nopt, opt, nout, out))
		return 1;
	
	return 0;
}

}
