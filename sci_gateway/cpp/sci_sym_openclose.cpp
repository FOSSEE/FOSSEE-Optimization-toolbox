// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Keyur Joshi, Iswarya
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#include <symphony.h>


extern sym_environment* global_sym_env;//defined in globals.cpp

extern "C" {
#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <localization.h>
#include <sciprint.h>

/* Function that initializes the symphony environment
 * Returns 1 on success , 0 on failure
 */
const char fname[] = "sym_open";
int sci_sym_open(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{

	// Error management variable
	SciErr sciErr;
	double status=0;
	
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

	//check environment
	if(global_sym_env==NULL){
	// 	sciprint("Warning: Symphony environment is already initialized.\n");
	// }else{
		global_sym_env = sym_open_environment();//open an environment
		if (!global_sym_env)
			sciprint("Error: Unable to create symphony environment.\n");
		else{
			status=1;


		}
	}

	/*write satus of function (success-1 or failure-0) as output argument to scilab*/
	out[0] = scilab_createDouble(env, status);
	
	return 0;
}

/*Function that closes symphony environment
 * Returns 1 on success , 0 on failure
*/

const char fname2[] = "sym_close";
int sci_sym_close(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{
	
	// Error management variable
	SciErr sciErr;
	double status=0;
	int output;//output parameter for closing the environment
    

  	//check whether we have no input and one output argument or not
	if (nin !=0)  //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname2, 0);
        	return STATUS_ERROR; 
	}
	
	if (nout !=1) //Checking the output arguments

	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname2, 1);
		return 1;
	}

	if (global_sym_env==NULL){//check for environment
		sciprint("Error: symphony environment is not initialized.\n");
	}else{
		output=sym_close_environment(global_sym_env);//close environment
		if(output==ERROR__USER){	
			status=0;//User error detected in user_free_master() function or when function invoked unsuccessfully
			sciprint("Error in user_free_master()\n");
		}else if(output==FUNCTION_TERMINATED_ABNORMALLY){
			status=0;//function invoked unsuccessfully
			sciprint("Symphony environment could not be closed.\n");
		}else if(output==FUNCTION_TERMINATED_NORMALLY){			
			status=1;//function invoked successfully and no error
			global_sym_env=NULL;//important to set to NULL, so that other functions can detect that environment is not open.
			//sciprint("Symphony environement closed successfully. Please run 'sym_open()' to restart.\n");
			//delete the sym_ variables

		}
	}

	/*write satus of function (success-1 or failure-0) as output argument to scilab*/
	out[0] = scilab_createDouble(env, status);
	
	return 0;	
}

}
