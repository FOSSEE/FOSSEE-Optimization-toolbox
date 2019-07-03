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
#include "sci_iofunc.hpp"
extern sym_environment* global_sym_env;//defined in globals.cpp

extern "C" {
#include "api_scilab.h"
#include "Scierror.h"
#include "BOOL.h"
#include "localization.h"
#include "sciprint.h"
#include "string.h"

int process_ret_val(int);

/*
 * This function returns the status of problem that has been solved.
 */

const char fname[] = "sym_get_status";

int sci_sym_get_status(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
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
	else // There is an environment opened
		status=sym_get_status(global_sym_env);// Call function
	// Return result to scilab
	out[0] = scilab_createDouble(env, status);
	return 0;

}


/* This is a generelized function for 
 * sym_isOptimal,sym_isInfeasible,sym_isAbandoned,
 * sym_isIterLimitReached,sym_isTimeLimitReached,
 * and sym_isTargetGapAchieved.
 * All the above functions have same return value and input argument.
 * 
 * It returns (to scilab) 
 * 1 if the function is proved true.
 * 0 if the function is proved false.
 * -1 if there is an error.
 */

const char fname[] = "sym_get_solver";

int sci_sym_get_solver(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{
	int result= -1 ;// Result to caller. Set to error.
  
	// Check whether we have no input and one output argument or not
	if (nin !=0)  //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname, 0);
        	return STATUS_ERROR; 
	}

	// One output argument (For scilab 1 o/p argument is fixed)	
	if (nout !=1) //Checking the output arguments

	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname, 1);
		return 1;
	}

	/* Array of possible callers of this function */
	char *arr_caller[]={"sym_isOptimal","sym_isInfeasible","sym_isAbandoned",
						"sym_isIterLimitReached","sym_isTimeLimitReached",
						"sym_isTargetGapAchieved" };
	
	/* Array of functions to be called */
	int (*fun[])(sym_environment *)= { sym_is_proven_optimal,
										sym_is_proven_primal_infeasible,
										sym_is_abandoned,sym_is_iteration_limit_reached,
										sym_is_time_limit_reached,sym_is_target_gap_achieved
											};
	/* Output values if functions return TRUE */
	char *output_true[] = {"The problem is solved to optimality.",
							"The problem is proved to be infeasible.",
							"The problem is abandoned.",
							"Iteration limit is reached.",
							"Time limit is reached.","Target gap is reached."};

	/* Output values if functions return FALSE */
	char *output_false[] = {"The problem is not solved to optimality.",
							"The problem is not proved to be infeasible.",
							"The problem is not abandoned.",
							"Iteration limit is not reached.",
							"Time limit is not reached.","Target gap is not reached."};
	// Check environment
	if(global_sym_env==NULL)
		sciprint("Error: Symphony environment is not initialized.\n");
	else 
	{//there is an environment opened
		int iter = 0, length= sizeof(arr_caller) / sizeof(char *),found_at= -1;

		for (;iter < length ;++iter)
			if (!strcmp(fname,arr_caller[iter])) //Find caller
				found_at=iter;
		if (found_at != -1 ) 
		{
			result = fun[found_at](global_sym_env);
			sciprint("\n");
			switch (result) 
			{
				case TRUE: // TRUE = 1
					sciprint(output_true[found_at]);
					break;
				case FALSE: // FALSE = 0
					sciprint(output_false[found_at]);
					break;
				default:
					sciprint("Undefined return value.");
					result = -1;
				}
			sciprint("\n");
			}
		else // Very rare case
			sciprint("\nError in function mapping in scilab script\n");
		}

	out[0] = scilab_createDouble(env, result);

	return 0;
	}
}
