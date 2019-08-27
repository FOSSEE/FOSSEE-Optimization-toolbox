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

extern "C" 
{
#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <localization.h>
#include <sciprint.h>
#include <stdio.h>
int process_ret_val(int);



const char fname[] = "sym_solve";
int sci_sym_solve(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{

	
	int status=0; 
  
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

	// Check environment
	if(global_sym_env==NULL)
		sciprint("Error: Symphony environment is not initialized.\n");
	else 
	{// There is an environment opened
		double time_limit = -1.0;
		status = sym_get_dbl_param(global_sym_env,"time_limit",&time_limit);

		if (status == FUNCTION_TERMINATED_NORMALLY) {
			if ( time_limit >= 0.0 )
				// sciprint("\nNote: There is no limit on time.\n");
			sciprint("\nNote: Time limit has been set to %lf.\n",time_limit);
			status=process_ret_val(sym_solve(global_sym_env));// Call function	
			}
		else 
		{
			sciprint("\nUnable to read time limit.\n");
			status = 1; //Error state
		}
	}
	// Return result to scilab
	out[0] = scilab_createDouble(env, status);

	return 0;
	}
}
