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

#include <iomanip>
#include <fstream>
#include <iostream>
#include "CoinPragma.hpp"
#include "CoinTime.hpp"
#include "CoinError.hpp"

#include "BonOsiTMINLPInterface.hpp"
#include "BonIpoptSolver.hpp"
#include "QuadTMINLP.hpp"
#include "BonCbc.hpp"
#include "BonBonminSetup.hpp"

#include "BonOACutGenerator2.hpp"
#include "BonEcpCuts.hpp"
#include "BonOaNlpOptim.hpp"

#include "sci_iofunc.hpp"
extern  "C"
{
#include <api_scilab.h>
#include <Scierror.h>
#include <BOOL.h>
#include <localization.h>
#include <sciprint.h>

int cpp_intqpipopt(char *fname)
{
  using namespace Ipopt;
  using namespace Bonmin;

	CheckInputArgument(pvApiCtx, 15, 15); // We need total 15 input arguments.
	CheckOutputArgument(pvApiCtx, 3, 3);  // 3 output arguments
	
	// Error management variable
	SciErr sciErr;

	// Input arguments
	double *QItems=NULL,*PItems=NULL, *intcon = NULL, *ConItems=NULL,*conUB=NULL,*conLB=NULL;
	double *varUB=NULL,*varLB=NULL,*init_guess=NULL,*options=NULL, *ifval=NULL;
	static unsigned int nVars = 0,nCons = 0, intconSize = 0;
	unsigned int temp1 = 0,temp2 = 0;
	char *bonmin_options_file = NULL;
	// Output arguments
	double *fX = NULL, ObjVal = 0,iteration=0;
	int rstatus = 0;
	
	//Number of Variables
	if(getIntFromScilab(1,&nVars))
	{
		return 1;
	}

	//Number of Constraints
	if (getIntFromScilab(2,&nCons))
	{
		return 1;
	}
	
	//Number of variables constrained to be integers
	if (getIntFromScilab(3,&intconSize))
	{
	  return 1;
	}

	//Q matrix from scilab
	temp1 = nVars;
	temp2 = nVars;
	if (getFixedSizeDoubleMatrixFromScilab(4,temp1,temp1,&QItems))
	{
		return 1;
	}
	
	//P matrix from scilab
	temp1 = 1;
	temp2 = nVars; 
	if (getFixedSizeDoubleMatrixFromScilab(5,temp1,temp2,&PItems))
	{
		return 1;
	}

	temp1 = 1;
	temp2 = intconSize;
	// Getting intcon
	if (getDoubleMatrixFromScilab(6,&temp1,&temp2,&intcon))
	{
		return 1;
	}

	if (nCons!=0)
	{
		//conMatrix matrix from scilab
		temp1 = nCons;
		temp2 = nVars;

		if (getFixedSizeDoubleMatrixFromScilab(7,temp1,temp2,&ConItems))
		{
			return 1;
		}

		//conLB matrix from scilab
		temp1 = 1;
		temp2 = nCons;
		if (getFixedSizeDoubleMatrixFromScilab(8,temp1,temp2,&conLB))
		{
			return 1;
		}

		//conUB matrix from scilab
		if (getFixedSizeDoubleMatrixFromScilab(9,temp1,temp2,&conUB))
		{
			return 1;
		}
	}

	//varLB matrix from scilab
	temp1 = 1;
	temp2 = nVars;
	if (getFixedSizeDoubleMatrixFromScilab(10,temp1,temp2,&varLB))
	{
		return 1;
	}

	//varUB matrix from scilab
	if (getFixedSizeDoubleMatrixFromScilab(11,temp1,temp2,&varUB))
	{
		return 1;
	}

	//Initial Value of variables from scilab
	if (getFixedSizeDoubleMatrixFromScilab( 12,temp1,temp2,&init_guess))
	{
		return 1;
	}

	temp1=1;
	temp2=5;
	if (getFixedSizeDoubleMatrixFromScilab(13,temp1,temp2,&options))
	{
		return 1;
	}
	
	temp1=1;
	temp2=5;
	if (getFixedSizeDoubleMatrixFromScilab(14,temp1,temp2,&ifval))
	{
		return 1;
	}
	
	if (getStringFromScilab(15, &bonmin_options_file))
	{
		return 1;
	}
  

  
		
  SmartPtr<QuadTMINLP> tminlp = new QuadTMINLP(nVars,nCons,intconSize,QItems, PItems, intcon,ConItems,conLB,conUB,varLB,varUB,init_guess);

  BonminSetup bonmin;
  bonmin.initializeOptionsAndJournalist();
            
  // Here we can change the default value of some Bonmin or Ipopt option
  bonmin.options()->SetStringValue("mu_oracle","loqo");
  

  //Register an additional option
	if((int)ifval[0])
            bonmin.options()->SetNumericValue("bonmin.integer_tolerance", (options[0]));
    if((int)ifval[1])
            bonmin.options()->SetIntegerValue("bonmin.node_limit", (int)(options[1]));
    if((int)ifval[2])
            bonmin.options()->SetNumericValue("bonmin.time_limit", (options[2]));
    if((int)ifval[3])
            bonmin.options()->SetNumericValue("bonmin.allowable_gap", (options[3]));
    if((int)ifval[4])
            bonmin.options()->SetIntegerValue("bonmin.iteration_limit", (int)(options[4]));


  //Here we read the option file
  if ( bonmin_options_file!=NULL )
    bonmin.readOptionsFile(bonmin_options_file);

  //Now initialize from tminlp
  bonmin.initialize(GetRawPtr(tminlp));
  
  //Set up done, now let's branch and bound
  try {
    Bab bb;
    bb(bonmin);//process parameter file using Ipopt and do branch and bound using Cbc
  }
  catch(TNLPSolver::UnsolvedError *E) {
    //There has been a failure to solve a problem with Ipopt.
    std::cerr<<"Ipopt has failed to solve a problem!"<<std::endl;
    sciprint(999, "\nIpopt has failed to solve the problem!\n");
  }
  catch(OsiTMINLPInterface::SimpleError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
	     <<std::endl
	     <<E.message()<<std::endl;
	  sciprint(999, "\nFailed to solve a problem!\n");
	}
  catch(CoinError &E) {
    std::cerr<<E.className()<<"::"<<E.methodName()
	     <<std::endl
	     <<E.message()<<std::endl;
	  sciprint(999, "\nFailed to solve a problem!\n");
	}
	rstatus=tminlp->returnStatus();
	if (rstatus >= 0 | rstatus <= 5){
		fX = tminlp->getX();
		ObjVal = tminlp->getObjVal();
		if (returnDoubleMatrixToScilab(1, 1, nVars, fX))
		{
			return 1;
		}

		if (returnDoubleMatrixToScilab(2, 1, 1, &ObjVal))
		{
			return 1;
		}

		if (returnIntegerMatrixToScilab(3, 1, 1, &rstatus))
		{
			return 1;
		}
		
	}
	else
	{
		if (returnDoubleMatrixToScilab(1, 0, 0, fX))
		{
			return 1;
		}

		if (returnDoubleMatrixToScilab(2, 1, 1, &ObjVal))
		{
			return 1;
		}

		if (returnIntegerMatrixToScilab(3, 1, 1, &rstatus))
		{
			return 1;
		}
		
		sciprint(999, "\nThe problem could not be solved!\n");
  }

	// As the SmartPtrs go out of scope, the reference count
	// will be decremented and the objects will automatically
	// be deleted(No memory leakage). 
	
  return 0;
}
}

