// Copyright (C) 2016 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Georgey John
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#include "sci_iofunc.hpp"

#include "ecos.h"
#include "ecos_bb.h"

extern "C"
{
	#include <api_scilab.h>
	#include <Scierror.h>
	#include <BOOL.h>
	#include <localization.h>
	#include <sciprint.h>

	// function to convert int to idxint
	idxint * int2idxint(int * sci_int, int n) {
		// int *sci_int;
		idxint *sci_idxint;
		sci_idxint = (idxint *) malloc(n * sizeof(idxint));
		for (int i = 0; i < n; i++) {
			sci_idxint[i] = sci_int[i];
			}
		return(sci_idxint);
	}

	int sci_ecos(char *fname){

		CheckInputArgument(pvApiCtx, 17, 17);
		CheckOutputArgument(pvApiCtx, 5, 5);
		// Error management variable
		SciErr sciErr;
		
		double *cptr=NULL, *Gprptr=NULL, *hptr=NULL, *Aprptr=NULL, *bptr=NULL;
		int  *Girptr=NULL,*Gjcptr=NULL, *Airptr=NULL, *Ajcptr=NULL, *qptr=NULL, mptr=0,nptr=0,pptr=0,nconesptr=0;
		static int  lptr=0, eptr=0;
		unsigned int temp1 = 0,temp2 = 0, iret = 0;
		double *maxit=NULL,*feastol=NULL,*reltol=NULL,*abstol=NULL,*feastol_inacc=NULL,*reltol_inacc=NULL,*abstol_inacc=NULL,*verbose=NULL,*mi_max_iter=NULL,*mi_int_tol=NULL,*mi_abs_eps=NULL,*mi_rel_eps=NULL;
		int c_rows=0, Gpr_rows=0, Gir_rows=0, Gjc_rows=0, h_rows=0, Apr_rows=0, Air_rows=0, Ajc_rows=0, b_rows=0, l_rows=0, q_rows=0, e_rows=0;
		int c_cols=0, Gpr_cols=0, Gir_cols=0, Gjc_cols=0, h_cols=0, Apr_cols=0, Air_cols=0, Ajc_cols=0, b_cols=0, l_cols=0, q_cols=0, e_cols=0;
	
		idxint exitflag = ECOS_FATAL;

		
		//Getting objective matrix 
		if(getDoubleMatrixFromScilab(1, &c_rows, &c_cols, &cptr))
		{
			return 1;
		}

		//Getting Gpr matrix linear inequality constraints 
		
		if(getDoubleMatrixFromScilab(2, &Gpr_rows, &Gpr_cols, &Gprptr))
		{
			return 1;
		}

		//Getting Gjc integer matrix representing coloumn indices of linear inequality constraints 
		if(getIntMatrixFromScilab(3, &Gjc_rows, &Gjc_cols, &Gjcptr))
		{
			return 1;
		}

		//Getting Gir integer matrix representing row indices of linear inequality constraints 
		if(getIntMatrixFromScilab(4, &Gir_rows, &Gir_cols, &Girptr))
		{
			return 1;
		}

		//Getting matrix representing RHS of linear inequality constraints 
		if(getDoubleMatrixFromScilab(5, &h_rows, &h_cols, &hptr))
		{
			return 1;
		}

		//Getting Apr matrix representing linear equality constraints 
		if(getDoubleMatrixFromScilab(6, &Apr_rows, &Apr_cols, &Aprptr))
		{
			return 1;
		}

		//Getting Ajc integer matrix representing coloumn indices linear inequality constraints 
		if(getIntMatrixFromScilab(7, &Ajc_rows, &Ajc_cols, &Ajcptr))
		{
			return 1;
		}

		//Getting Air integer matrix representing row indices of linear inequality constraints 
		if(getIntMatrixFromScilab(8, &Air_rows, &Air_cols, &Airptr))
		{
			return 1;
		}

		//Getting matrix representing RHS of linear inequality constraints 
		if(getDoubleMatrixFromScilab(9, &b_rows, &b_cols, &bptr))
		{
			return 1;
		}

		//Getting Integer representing dimensions of positive orthant
		if(getIntFromScilab(10, &lptr))
		{
			return 1;
		}

		//Getting integer matrix representing dimensions of each cone  
		if(getIntMatrixFromScilab(11, &q_rows, &q_cols, &qptr))
		{
			return 1;
		}

		//Getting integer representing number of exponential cone
		if(getIntFromScilab(12, &eptr))
		{
			return 1;
		}
		
		//Getting matrix representing maximum iteration
		if(getFixedSizeDoubleMatrixInList(13,2,temp1,temp2,&maxit))
		{
			return 1;
		}

		//Getting matrix representing feasible tolerance
		if(getFixedSizeDoubleMatrixInList(13,4,temp1,temp2,&feastol))
		{
			return 1;
		}

		//Getting matrix representing relative tolerance
		if(getFixedSizeDoubleMatrixInList(13,6,temp1,temp2,&reltol))
		{
			return 1;
		}

		//Getting matrix representing absolute tolerance
		if(getFixedSizeDoubleMatrixInList(13,8,temp1,temp2,&abstol))
		{
			return 1;
		}

		//Getting matrix representing feastol_inacc
		if(getFixedSizeDoubleMatrixInList(13,10,temp1,temp2,&feastol_inacc))
		{
			return 1;
		}

		//Getting matrix representing reltol_inacc
		if(getFixedSizeDoubleMatrixInList(13,12,temp1,temp2,&reltol_inacc))
		{
			return 1;
		}

		//Getting matrix representing abstol_inacc
		if(getFixedSizeDoubleMatrixInList(13,14,temp1,temp2,&abstol_inacc))
		{
			return 1;
		}
		//Getting matrix representing verbose mode
		if(getFixedSizeDoubleMatrixInList(13,16,temp1,temp2,&verbose))
		{
			return 1;
		}

		// /ecos_bb options
		// if(getFixedSizeDoubleMatrixInList(13,18,temp1,temp2,&mi_max_iter))
		// {
		// 	return 1;
		// }

		// if(getFixedSizeDoubleMatrixInList(13,20,temp1,temp2,&mi_int_tol))
		// {
		// 	return 1;
		// }

		// if(getFixedSizeDoubleMatrixInList(13,22,temp1,temp2,&mi_rel_eps))
		// {
		// 	return 1;
		// }

		// if(getFixedSizeDoubleMatrixInList(13,24,temp1,temp2,&mi_abs_eps))
		// {
		// 	return 1;
		// }

		// //Getting integer representing number of inequality constraints
		if(getIntFromScilab(14, &mptr))
		{
			return 1;
		}
		
		//Getting integer representing the number of primal variable 
		if(getIntFromScilab(15, &nptr))
		{
			return 1;
		}

		//Getting integer representing the number of equality constraints
		if(getIntFromScilab(16, &pptr))
		{
			return 1;
		}

		//Getting integer representing the number of second order cones
		if(getIntFromScilab(17, &nconesptr))
		{
			return 1;
		}

		// Intializing ECOS variables
		pfloat *Gpr=NULL,*c=NULL,*h=NULL,*Apr=NULL,*b=NULL;
		idxint *Gjc=NULL,*Gir=NULL,*Ajc=NULL,*Air=NULL, n=0, m=0, p=0, l=0, *q=NULL,ncones=0, e=0;
		double *x = NULL, *y = NULL, *s = NULL, *z = NULL;
		const char* infostring;
		
		c = (pfloat *)cptr;

		if (Gpr_rows != 0)
		{
			Gpr = (pfloat *) Gprptr;
			Gjc = int2idxint(Gjcptr,Gjc_rows);
			Gir = int2idxint(Girptr,Gir_rows);
			h = (pfloat *) hptr;
		}		

		if (Apr_rows != 0)
		{
			Apr = (pfloat *)Aprptr;
			Ajc = int2idxint(Ajcptr,Ajc_rows);
			Air = int2idxint(Airptr,Air_rows);
			b = (pfloat *)bptr;
		}
		
		n = (idxint)nptr;
		m = (idxint)mptr;
		p = (idxint)pptr;
		l = (idxint)lptr;
		ncones = (idxint)nconesptr;

		if (ncones != 0)
		{
			q = int2idxint(qptr,ncones);
		}

		e = (idxint)eptr;


		pwork* input_struct;

		// for (int i = 0; i < Gpr_rows; i++)
		// {
		// 	sciprint("%f\t",Gpr[i]);
		// }
		// sciprint("\n");

		// for (int i = 0; i < Gjc_rows; i++)
		// {
		// 	sciprint("%d\t",Gjc[i]);
		// }
		// sciprint("\n");

		// for (int i = 0; i < Gir_rows; i++)
		// {
		// 	sciprint("%d\t",Gir[i]);
		// }
		// sciprint("\n");
		// sciprint("%d\t",n);sciprint("%d\t",m);sciprint("%d\t",p);sciprint("%d\t",l);sciprint("%d\t",ncones);

		// setup ecos 
		input_struct = ECOS_setup(n, m, p, l, ncones, q, e, Gpr, Gjc, Gir, Apr, Ajc, Air,c, h, b);

		// if setup is successful
		if( input_struct != NULL ){
			// passing options to ecos
			input_struct->stgs->feastol 			= (pfloat) *feastol;
			input_struct->stgs->reltol 				= (pfloat) *reltol;
			input_struct->stgs->abstol 				= (pfloat) *abstol;
			input_struct->stgs->feastol_inacc 		= (pfloat) *feastol_inacc;
			input_struct->stgs->abstol_inacc 		= (pfloat) *abstol_inacc;
			input_struct->stgs->reltol_inacc 		= (pfloat) *reltol_inacc;
			input_struct->stgs->maxit 				= (idxint) maxit;
			input_struct->stgs->verbose 			= (idxint) verbose;
       		
       		// solve
	        exitflag = ECOS_solve(input_struct); 

	        // output
	        x = input_struct->x;
			y = input_struct->y;
			s = input_struct->s;
			z = input_struct->z;

			double pcost = (double)input_struct->info->pcost;
			double dcost = (double)input_struct->info->dcost;
			double pres = (double)input_struct->info->pres;
			double dres = (double)input_struct->info->dres;
			double pinf = (double)input_struct->info->pinf;
			double dinf = (double)input_struct->info->dinf;
			double pinfres = (double)input_struct->info->pinfres;
			double dinfres = (double)input_struct->info->dinfres;
			double gap = (double)input_struct->info->gap;
			double rel_gap = (double)input_struct->info->relgap;
			double iter = (double)input_struct->info->iter;

			// #if PROFILING > 1
			// 	double tkktcreate = (double)input_struct->info->tkktcreate;
			// 	double tkktsolve = (double)input_struct->info->tkktsolve;
			// 	double tkktfactor = (double)input_struct->info->tfactor;
			// 	double torder = (double)input_struct->info->torder;
			// 	double ttranspose = (double)input_struct->info->ttranspose;
			// #else
				double tsetup = (double)input_struct->info->tsetup;
				double tsolve = (double)input_struct->info->tsolve;
			// #endif

			// #if PROFILING > 1
			// 	double info[19]={tsetup,tsolve,pcost,dcost,pres,dres,pinf,dinf,pinfres,dinfres,gap,rel_gap,iter,(double) exitflag,tkktcreate,tkktsolve,tfactor,torder,ttranspose};
			// 	int o = 19;
			// #else
				double info[14]={tsetup,tsolve,pcost,dcost,pres,dres,pinf,dinf,pinfres,dinfres,gap,rel_gap,iter,(double) exitflag};
			// 	int o = 14;
			// #endif

			if (returnDoubleMatrixToScilab(1, n, 1, x))
			{
				return 1;
			}

			if (returnDoubleMatrixToScilab(2, p, 1, y))
			{
				return 1;
			}

			if (returnDoubleMatrixToScilab(3, 14, 1, info))
			{
				return 1;
			}

			if (returnDoubleMatrixToScilab(4, m, 1, s))
			{
				return 1;
			}

			if (returnDoubleMatrixToScilab(5, m, 1, z))
			{
				return 1;
			}

			/* clean up memory */
			ECOS_cleanup(input_struct, 0);
	    }
	    else{
	    	// if setup fails
	    	exitflag =-8;
			
			// #if PROFILING > 1
			// 	double info[19]={0,0,0,0,0,0,0,0,0,0,0,0,0,(double) exitflag,0,0,0,0,0};
			// 	int o = 19;
			// #else
				double info[14]={0,0,0,0,0,0,0,0,0,0,0,0,0,(double) exitflag};
			// 	int o = 14;
			// #endif

			if (returnDoubleMatrixToScilab(1, 0, 0, x))
			{
				return 1;
			}

			if (returnDoubleMatrixToScilab(2, 0, 0, y))
			{
				return 1;
			}

			if (returnDoubleMatrixToScilab(3, 14, 1, info))
			{
				return 1;
			}

			if (returnDoubleMatrixToScilab(4, 0, 0, s))
			{
				return 1;
			}

			if (returnDoubleMatrixToScilab(5, 0, 0, z))
			{
				return 1;
			}
	    }	

	}
}