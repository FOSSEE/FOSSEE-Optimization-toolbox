// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Sai Kiran
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#include "symphony.h"

extern sym_environment* global_sym_env;//defined in globals.cpp

extern "C" {
#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <stdlib.h>
#include <malloc.h>
#include <localization.h>
#include <sciprint.h>

#include <string.h>

/*
 * This function is used to get iteration count after solving a problem
*/

const char fname[] = "sym_get_iteration_count";

int sci_sym_get_iteration_count(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{
	
	//check whether we have no input and one output argument or not
	if (nin !=0)  //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname, 0);
        	return STATUS_ERROR; 
	}
	
	if (nout !=1) //Checking the output arguments
	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname, 1);
		return 1;
	}
	
	int iteration_count=0; // return value to the caller
	if(global_sym_env==NULL) //There is no environment opened.
		sciprint("Error: Symphony environment is not initialized.\n");
	else { //There is an environment opened
		 //Call symphony function
		int status=sym_get_iteration_count(global_sym_env,&iteration_count);
		if (status == FUNCTION_TERMINATED_ABNORMALLY) {
			sciprint("\nHave you solved a problem ?\n");
			iteration_count = 0;
			}
		}
	// Write the result to scilab

	out[0] = scilab_createDouble(env, iteration_count);
	return 0;
	}

}
