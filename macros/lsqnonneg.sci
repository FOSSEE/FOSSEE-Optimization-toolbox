// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Harpreet Singh
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in


function [xopt,resnorm,residual,exitflag,output,lambda] = lsqnonneg (varargin)
	// Solves nonnegative least-squares curve fitting problems.
	//
	//   Calling Sequence
	//   xopt = lsqnonneg(C,d)
	//   xopt = lsqnonneg(C,d,param)
	//   [xopt,resnorm,residual,exitflag,output,lambda] = lsqnonneg( ... )
	//   
	//   Parameters
	//   C : a matrix of double, represents the multiplier of the solution x in the expression C⋅x - d. Number of columns in C is equal to the number of elements in x.
	//   d : a vector of double, represents the additive constant term in the expression C⋅x - d. Number of elements in d is equal to the number of rows in C matrix.
	//   xopt : a vector of double, the computed solution of the optimization problem.
	//   resnorm : a double, objective value returned as the scalar value norm(C⋅x-d)^2.
	//   residual : a vector of double, solution residuals returned as the vector d-C⋅x.
	//   exitflag : The exit status. See below for details.
	//   output : The structure consist of statistics about the optimization. See below for details.
	//   lambda : The structure consist of the Lagrange multipliers at the solution of problem. See below for details.
	//   
	//   Description
	//   Solves nonnegative least-squares curve fitting problems specified by :
	//
	//   <latex>
	//    \begin{eqnarray}
	//    &\mbox{min}_{x}
	//    & 1/2||C⋅x - d||_2^2  \\
	//    & & x \geq 0 \\
	//    \end{eqnarray}
	//   </latex>
	//   
	//   The routine calls Ipopt for solving the nonnegative least-squares curve fitting problems, Ipopt is a library written in C++.
	//
  	// The options allows the user to set various parameters of the Optimization problem. 
  	// It should be defined as type "list" and contains the following fields.
	// <itemizedlist>
	//   <listitem>Syntax : options= list("MaxIter", [---], "CpuTime", [---]);</listitem>
	//   <listitem>MaxIter : a Scalar, containing the Maximum Number of Iteration that the solver should take.</listitem>
	//   <listitem>CpuTime : a Scalar, containing the Maximum amount of CPU Time that the solver should take.</listitem>
	//   <listitem>Default Values : options = list("MaxIter", [3000], "CpuTime", [600]);</listitem>
	// </itemizedlist>
	//
	// The exitflag allows to know the status of the optimization which is given back by Ipopt.
	// <itemizedlist>
	//   <listitem>exitflag=0 : Optimal Solution Found </listitem>
	//   <listitem>exitflag=1 : Maximum Number of Iterations Exceeded. Output may not be optimal.</listitem>
	//   <listitem>exitflag=2 : Maximum CPU Time exceeded. Output may not be optimal.</listitem>
	//   <listitem>exitflag=3 : Stop at Tiny Step.</listitem>
	//   <listitem>exitflag=4 : Solved To Acceptable Level.</listitem>
	//   <listitem>exitflag=5 : Converged to a point of local infeasibility.</listitem>
	// </itemizedlist>
	//
	// For more details on exitflag see the ipopt documentation, go to http://www.coin-or.org/Ipopt/documentation/
	//
	// The output data structure contains detailed informations about the optimization process. 
	// It has type "struct" and contains the following fields.
	// <itemizedlist>
	//   <listitem>output.iterations: The number of iterations performed during the search</listitem>
	//   <listitem>output.constrviolation: The max-norm of the constraint violation.</listitem>
	// </itemizedlist>
	//
	// The lambda data structure contains the Lagrange multipliers at the end 
	// of optimization. In the current version the values are returned only when the the solution is optimal. 
	// It has type "struct" and contains the following fields.
	// <itemizedlist>
	//   <listitem>lambda.lower: The Lagrange multipliers for the lower bound constraints.</listitem>
	//   <listitem>lambda.upper: The Lagrange multipliers for the upper bound constraints.</listitem>
	// </itemizedlist>
	//    
	// Examples 
	// // A basic lsqnonneg problem
	//		C = [1 1 1;
	//			1 1 0;
	//			0 1 1;
	//			1 0 0;
	//			0 0 1]
	//		d = [89;
	//			67;
	//			53;
	//			35;
	//			20;]
	// [xopt,resnorm,residual,exitflag,output,lambda] = lsqnonneg(C,d)
	// Authors
	// Harpreet Singh


	//To check the number of input and output argument
	[lhs , rhs] = argn();

	//To check the number of argument given by user
	if ( rhs < 2 | rhs > 3 ) then
		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be in the set of [2 3]"), "lsqnonneg", rhs);
		error(errmsg)
	end

	C = varargin(1);
	d = varargin(2);
	nbVar = size(C,2);

	if(nbVar == 0) then
		errmsg = msprintf(gettext("%s: Cannot determine the number of variables because input objective coefficients is empty"), "lsqnonneg");
		error(errmsg);
	end

	if ( rhs<3 | size(varargin(3)) ==0 ) then
		param = list();
	else
		param =varargin(3);
	end

	if (type(param) ~= 15) then
		errmsg = msprintf(gettext("%s: param should be a list "), "lsqnonneg");
		error(errmsg);
	end

	//Check type of variables
	Checktype("lsqnonneg", C, "C", 1, "constant")
	Checktype("lsqnonneg", d, "d", 2, "constant")

	if (modulo(size(param),2)) then
		errmsg = msprintf(gettext("%s: Size of parameters should be even"), "lsqnonneg");
		error(errmsg);
	end

	options = list(	"MaxIter"   , [3000], ...
					"CpuTime"   , [600] ...
	);

	for i = 1:(size(param))/2

		select convstr(param(2*i-1),'l')
			case "maxiter" then
				options(2*i) = param(2*i);
			case "cputime" then
				options(2*i) = param(2*i);
			else
				errmsg = msprintf(gettext("%s: Unrecognized parameter name ''%s''."), "lsqlin", param(2*i-1));
				error(errmsg)
		end
	end

	// Check if the user gives row vector 
	// and Changing it to a column matrix


	if (size(d,2)== [nbVar]) then
		d=d';
	end

	//Check the size of f which should equal to the number of variable
	if ( size(d,1) ~= size(C,1)) then
		errmsg = msprintf(gettext("%s: The number of rows in C must be equal the number of elements of d"), "lsqnonneg");
		error(errmsg);
	end

	//Converting it into Quadratic Programming Problem

	Q = C'*C;
	p = [-C'*d]';
	op_add = d'*d;
	lb = repmat(0,1,nbVar);
	ub = repmat(%inf,1,nbVar);	
	x0 = repmat(0,1,nbVar);;
	conMatrix = [];
	nbCon = size(conMatrix,1);
	conLB = [];
	conUB = [] ; 
	[xopt,fopt,status,iter,Zl,Zu,lmbda] = solveqp(int32(nbVar),int32(nbCon),Q,p,conMatrix,conLB,conUB,lb,ub,x0,options);

	xopt = xopt';
	residual = -1*(C*xopt-d);
	resnorm = residual'*residual;
	exitflag = status;
	output = struct("Iterations"      , [], ..
					"ConstrViolation" ,[]);
	output.Iterations = iter;
	output.ConstrViolation = max([0;(lb'-xopt);(xopt-ub')]);

	lambda = struct("lower"           , [], ..
		           "upper"           , []);
	lambda.lower = Zl;
	lambda.upper = Zu;

	select status 
		case 0 then
			printf("\nOptimal Solution Found.\n");
		case 1 then
			printf("\nMaximum Number of Iterations Exceeded. Output may not be optimal.\n");
		case 2 then
			printf("\nMaximum CPU Time exceeded. Output may not be optimal.\n");
		case 3 then
			printf("\nStop at Tiny Step\n");
		case 4 then
			printf("\nSolved To Acceptable Level\n");
		case 5 then
			printf("\nConverged to a point of local infeasibility.\n");
		case 6 then
			printf("\nStopping optimization at current point as requested by user.\n");
		case 7 then
			printf("\nFeasible point for square problem found.\n");
		case 8 then 
			printf("\nIterates diverging; problem might be unbounded.\n");
		case 9 then
			printf("\nRestoration Failed!\n");
		case 10 then
			printf("\nError in step computation (regularization becomes too large?)!\n");
		case 12 then
			printf("\nProblem has too few degrees of freedom.\n");
		case 13 then
			printf("\nInvalid option thrown back by Ipopt\n");
		case 14 then
			printf("\nNot enough memory.\n");
		case 15 then
			printf("\nINTERNAL ERROR: Unknown SolverReturn value - Notify Ipopt Authors.\n");
		else
			printf("\nInvalid status returned. Notify the Toolbox authors\n");
		break;
	end

endfunction
