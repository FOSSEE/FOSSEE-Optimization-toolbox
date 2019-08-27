// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Iswarya
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#include <symphony.h>

extern sym_environment* global_sym_env;//defined in globals.cpp

extern "C" {
#include <stdlib.h>
#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <localization.h>
#include <sciprint.h>
//This function is for loading a mps file to symphony

const char fname[] = "sym_set_defaults";
int sci_sym_set_defaults(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{
	
	double status=1.0;//assume error status
	int output=0;//out parameter for the setting of default values  function

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

	//ensure that environment is active
	if(global_sym_env==NULL)
	{
		sciprint("Error: Symphony environment not initialized. Please run 'sym_open()' first.\n");
	}
	else 
	{
		output=sym_set_defaults(global_sym_env);//setting all environment variables and parameters in this symphony environment passed to their default values
		if(output==FUNCTION_TERMINATED_ABNORMALLY)
		{
			status=1.0;//function did not invoke successfully
			sciprint("Function terminated abnormally, didnot execute");
		}
		else if(output==FUNCTION_TERMINATED_NORMALLY)
		{
			status=0.0;//no error in executing the function
			sciprint("Function executed successfully");
		}
		
		
	}
	
	out[0] = scilab_createDouble(env, status);

	return 0;
}



const char fname2[] = "sym_set_int_param";
int sci_sym_set_int_param(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{
	
	// Error management variable

	double status=1.0;//assume error status
	double num;//to store and check the value obtained is pure unsigned integer
	int output;//output variable to store the return value of symphony set integer function
	int value;//to store the value of integer to be set
	char variable_name[100];//string to hold the name of variable's value to be set
	char *ptr=variable_name;//pointer to point to address of the variable name
	wchar_t *sciPtr; 

	int ret;


	if (nin !=2) //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname2, 2);
        	return STATUS_ERROR; 
	}
	
	if (nout !=1) //Checking the output arguments
	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname2, 1);
		return 1;
	}


	if (scilab_isString(env, in[0]) == 0 || scilab_isScalar(env, in[0]) == 0) //Get the first input variable i.e. the option name
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: A function expected.\n", fname2, 1);
   		return 1;
	}	

	scilab_getString(env, in[0], &sciPtr);

	ret = wcstombs(ptr, sciPtr, 100);

	if (scilab_isDouble(env, in[1]) == 0 || scilab_isScalar(env, in[1]) == 0) //Get the second input variable i.e. the double option value
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: An int expected.\n", fname2, 2);
    	return 1;
	}

	scilab_getDouble(env, in[1], &num);

	

	//check for the integrity of integer value obtained
	if((double)(unsigned int)num!=(double)num)
		return 0;
	else
		value=(unsigned int)num;//if the value passed is an integer ,store it as an unsigned integer in value variable 

	//ensure that environment is active
	if(global_sym_env==NULL)
	{
		sciprint("Error: Symphony environment not initialized. Please run 'sym_open()' first.\n");
	}
	else 
	{
		output=sym_set_int_param(global_sym_env,ptr,value);//symphony function to set the variable name pointed by the ptr pointer to the integer value stored in value variable.
		if(output==FUNCTION_TERMINATED_NORMALLY)
		{
			sciprint("setting of integer parameter function executed successfully\n");
			status=0.0;	
		}
		else if(output==FUNCTION_TERMINATED_ABNORMALLY)
		{
			sciprint("setting of integer parameter was unsuccessful.....check your parameter and value\n");
			status=1.0;
		}
		else
			sciprint("\nerror while executing the setting integer function...check your parameter and value!!\n");
		
	}
	
	out[0] = scilab_createDouble(env, status);
	//ReturnArguments(pvApiCtx);

	return 0;
	}

const char fname3[] = "sym_get_int_param";
int sci_sym_get_int_param(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{

	
	// Error management variable
	double status=1.0;//assume error status
	char variable_name[100];//string to hold the name of variable's value to be retrieved
	char *ptr=variable_name;//pointer to point to address of the variable name
	int output;//output parameter for the symphony get_int_param function

	wchar_t* sciPtr;
	int ret;

	if (nin !=1) //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname3, 1);
        	return STATUS_ERROR; 
	}
	
	if (nout !=1) //Checking the output arguments
	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname3, 1);
		return 1;
	}


	if (scilab_isString(env, in[0]) == 0 || scilab_isScalar(env, in[0]) == 0) //Get the first input variable i.e. the option name
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: A function expected.\n", fname3, 1);
   		return 1;
	}	

	scilab_getString(env, in[0], &sciPtr);

	ret = wcstombs(ptr, sciPtr, 100);
	

	//ensure that environment is active
	if(global_sym_env==NULL){
		sciprint("Error: Symphony environment not initialized. Please run 'sym_open()' first.\n");
		}
	else {
		int a;//local variable to store the value of variable name we want to retrieve
		output=sym_get_int_param(global_sym_env,ptr,&a);//symphony function to get the value of integer parameter pointed by ptr pointer and store it in 'a' variable
		if(output==FUNCTION_TERMINATED_NORMALLY){			
			sciprint("value of integer parameter %s is :: %d\n",ptr,a);
			status=0.0;
		}
		else{
			sciprint("Unable to get the value of the parameter...check the input values!!\n");
			status=1.0;
		} 
				
		}
	
	out[0] = scilab_createDouble(env, status);
	//ReturnArguments(pvApiCtx);

	return 0;
	}


const char fname4[] = "sym_set_dbl_param";
int sci_sym_set_dbl_param(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{

	
	// Error management variable
	double status=1.0;//assume error status
	double num;//to store the value of the double parameter to be set
	int output;//output parameter for the symphomy setting double parameter function
	char variable_name[100];//string to hold the name of variable's value to be set
	char *ptr=variable_name;//pointer to point to address of the variable name

	wchar_t* sciPtr;
	int ret;

	if (nin !=2) //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname4, 2);
        	return STATUS_ERROR; 
	}
	
	if (nout !=1) //Checking the output arguments
	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname4, 1);
		return 1;
	}

	if (scilab_isString(env, in[0]) == 0 || scilab_isScalar(env, in[0]) == 0) //Get the first input variable i.e. the option name
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: A function expected.\n", fname4, 1);
   		return 1;
	}	

	scilab_getString(env, in[0], &sciPtr);

	ret = wcstombs(ptr, sciPtr, 100);

	if (scilab_isDouble(env, in[1]) == 0 || scilab_isScalar(env, in[1]) == 0) //Get the second input variable i.e. the double option value
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: An int expected.\n", fname4, 2);
    	return 1;
	}

	scilab_getDouble(env, in[1], &num);

	//ensure that environment is active
	if(global_sym_env==NULL){
		sciprint("Error: Symphony environment not initialized. Please run 'sym_open()' first.\n");
		}
	else {
		output=sym_set_dbl_param(global_sym_env,ptr,num);//symphony function to set the variable name pointed by the ptr pointer to the double value stored in 'value' variable.
		if(output==FUNCTION_TERMINATED_NORMALLY){		
			// sciprint("setting of double parameter function executed successfully\n");
			status=0.0;
		}
		else
			sciprint("Function did not execute successfully...check your inputs!!!\n");
		
		}
	
	out[0] = scilab_createDouble(env, status);
	//ReturnArguments(pvApiCtx);

	return 0;
	}

const char fname5[] = "sym_get_dbl_param";
int sci_sym_get_dbl_param(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{

	
	// Error management variable
	double status=1.0;//assume error status
	char variable_name[100];//string to hold the name of variable's value to be retrieved
	char *ptr=variable_name;//pointer to point to address of the variable name
	int output;//output parameter for the symphony get_dbl_param function

	wchar_t* sciPtr;
	int ret;
	
	if (nin !=1) //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname5, 2);
        	return STATUS_ERROR; 
	}
	
	if (nout !=1) //Checking the output arguments
	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname5, 1);
		return 1;
	}

	if (scilab_isString(env, in[0]) == 0 || scilab_isScalar(env, in[0]) == 0) //Get the first input variable i.e. the option name
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: A function expected.\n", fname5, 1);
   		return 1;
	}	

	scilab_getString(env, in[0], &sciPtr);

	ret = wcstombs(ptr, sciPtr, 100);
	

	//ensure that environment is active
	if(global_sym_env==NULL){
		sciprint("Error: Symphony environment not initialized. Please run 'sym_open()' first.\n");
		}
	else {
		double a;//local variable to store the value of variable name we want to retrieve
		output=sym_get_dbl_param(global_sym_env,ptr,&a);//symphony function to get the value of double parameter pointed by ptr pointer and store it in 'a' variable
		if(output==FUNCTION_TERMINATED_NORMALLY){
			sciprint("value of double parameter %s is :: %lf\n",ptr,a);
			status=1.0;
		}
		else{
			sciprint("Unable to get the value of the parameter...check the input values!!\n");
			status=1.0;
		} 


		}
	
	out[0] = scilab_createDouble(env, status);
	//ReturnArguments(pvApiCtx);

	return 0;
	}


const char fname6[] = "sym_set_str_param";
int sci_sym_set_str_param(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{

	
	// Error management variable
	double status=1.0;//assume error status
	double num;//to store the value of the double parameter to be set
	int output;//output return value of the setting of symphony string parameter function
	char variable_name[100],value[100];//string to hold the name of variable's value to be set and the value to be set is stored in 'value' string
	char *ptr=variable_name,*valptr=value;//pointer-'ptr' to point to address of the variable name and 'valptr' points to the address of the value to be set to the string parameter

	wchar_t* sciPtr;
	wchar_t* sciValPtr;
	int ret;

	if (nin !=2) //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname6, 2);
        	return STATUS_ERROR; 
	}
	
	if (nout !=1) //Checking the output arguments
	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname6, 1);
		return 1;
	}


	if (scilab_isString(env, in[0]) == 0 || scilab_isScalar(env, in[0]) == 0) //Get the first input variable i.e. the option name
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: A function expected.\n", fname6, 1);
   		return 1;
	}	

	scilab_getString(env, in[0], &sciPtr);

	ret = wcstombs(ptr, sciPtr, 100);

	if (scilab_isString(env, in[1]) == 0 || scilab_isScalar(env, in[1]) == 0) //Get the second input variable i.e. the double option value
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: An int expected.\n", fname6, 2);
    	return 1;
	}

	scilab_getString(env, in[1], &sciValPtr);
	
	ret = wcstombs(valptr, sciValPtr, 100);

	//ensure that environment is active
	if(global_sym_env==NULL){
		sciprint("Error: Symphony environment not initialized. Please run 'sym_open()' first.\n");
		}
	else {
		output=sym_set_str_param(global_sym_env,ptr,valptr);//symphony function to set the variable name pointed by the ptr pointer to the double value stored in 'value' variable.
		if(output==FUNCTION_TERMINATED_NORMALLY){
			sciprint("setting of string parameter function executed successfully\n");
			status=0.0;
		}
		else
			sciprint("Setting of the string parameter was unsuccessfull...check the input values!!\n");
		
		}
	

	out[0] = scilab_createDouble(env, status);
	//ReturnArguments(pvApiCtx);

	return 0;
	}

const char fname7[] = "sym_get_str_param";
int sci_sym_get_str_param(scilabEnv env, int nin, scilabVar* in, int nopt, scilabOpt opt, int nout, scilabVar* out)
{

	
	// Error management variable
	SciErr sciErr1;
	double status=1.0;//assume error status
	int *piAddressVarOne = NULL;//pointer used to access first argument of the function
	char variable_name[100];//string to hold the name of variable's value to be retrieved
	char *ptr=variable_name;//pointer to point to address of the variable name
	int output;//output parameter for the symphony get_dbl_param function

	wchar_t* sciPtr;
	int ret;
	

	if (nin !=1) //Checking the input arguments
	{
        	Scierror(999, "%s: Wrong number of input arguments: %d expected.\n", fname7, 1);
        	return STATUS_ERROR; 
	}
	
	if (nout !=1) //Checking the output arguments
	{
		Scierror(999, "%s: Wrong number of output argument(s): %d expected.\n", fname7, 1);
		return 1;
	}

	if (scilab_isString(env, in[0]) == 0 || scilab_isScalar(env, in[0]) == 0) //Get the first input variable i.e. the option name
	{
    	Scierror(999, "%s: Wrong type for input argument #%d: A function expected.\n", fname7, 1);
   		return 1;
	}	

	scilab_getString(env, in[0], &sciPtr);

	ret = wcstombs(ptr, sciPtr, 100);
	

	//ensure that environment is active
	if(global_sym_env==NULL){
		sciprint("Error: Symphony environment not initialized. Please run 'sym_open()' first.\n");
		}
	else {
		char value[100];//local variable to store the value of variable name we want to retrieve
		char *p=value;//pointer to store the address of the character array that contains the value of the string parameter to be retrieved
		output=sym_get_str_param(global_sym_env,ptr,&p);//symphony function to get the value of string parameter pointed by ptr pointer and store it in 'p' pointer variable
		if(output==FUNCTION_TERMINATED_NORMALLY){
			sciprint("value of string parameter %s is :: %s\n",ptr,p);
			status=0.0;
		}
		else
			sciprint("The string parameter value could not be retrieved successfully...check the input values!!\n");

				
		}
	
	out[0] = scilab_createDouble(env, status);
	//ReturnArguments(pvApiCtx);

	return 0;
	}


}
