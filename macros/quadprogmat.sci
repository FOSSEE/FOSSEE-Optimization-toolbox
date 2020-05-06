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


function [xopt,fopt,exitflag,output,lambda] = quadprogmat (varargin)
	// Solves a linear quadratic problem.
	//
	//   Calling Sequence
	//   xopt = quadprogmat(H,f)
	//   xopt = quadprogmat(H,f,A,b)
	//   xopt = quadprogmat(H,f,A,b,Aeq,beq)
	//   xopt = quadprogmat(H,f,A,b,Aeq,beq,lb,ub)
	//   xopt = quadprogmat(H,f,A,b,Aeq,beq,lb,ub,x0)
	//   xopt = quadprogmat(H,f,A,b,Aeq,beq,lb,ub,x0,param)
	//   [xopt,fopt,exitflag,output,lamda] = quadprogmat( ... )
	//   
	//   Parameters
	//   H : a symmetric matrix of double, represents coefficients of quadratic in the quadratic problem.
	//   f : a vector of double, represents coefficients of linear in the quadratic problem
	//   A : a matrix of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b. 
	//   b : a vector of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b.
	//   Aeq : a matrix of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   beq : a vector of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   lb : a vector of double, contains lower bounds of the variables.
	//   ub : a vector of double, contains upper bounds of the variables.
	//   x0 : a vector of double, contains initial guess of variables.
	//   param : a list containing the parameters to be set.
	//   xopt : a vector of double, the computed solution of the optimization problem.
	//   fopt : a double, the value of the function at x.
	//   exitflag : The exit status. See below for details.
	//   output : The structure consist of statistics about the optimization. See below for details.
	//   lambda : The structure consist of the Lagrange multipliers at the solution of problem. See below for details.
	//   
	//   Description
	//   Search the minimum of a constrained linear quadratic optimization problem specified by :
	//
	//   <latex>
	//    \begin{eqnarray}
	//    &\mbox{min}_{x}
	//    & 1/2⋅x^T⋅H⋅x + f^T⋅x  \\
	//    & \text{subject to} & A⋅x \leq b \\
	//    & & Aeq⋅x = beq \\
	//    & & lb \leq x \leq ub \\
	//    \end{eqnarray}
	//   </latex>
	//   
	//   The routine calls Ipopt for solving the quadratic problem, Ipopt is a library written in C++.
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
	//   <listitem>lambda.eqlin: The Lagrange multipliers for the linear equality constraints.</listitem>
	//   <listitem>lambda.ineqlin: The Lagrange multipliers for the linear inequality constraints.</listitem>
	// </itemizedlist>
	//
	// Examples
	//		//Ref : example 14 :
	//		//https://www.me.utexas.edu/~jensen/ORMM/supplements/methods/nlpmethod/S2_quadratic.pdf
	//		// min. -8*x1*x1 -16*x2*x2 + x1 + 4*x2
	//		// such that
	//		//	x1 + x2 <= 5,
	//		//	x1 <= 3,
	//		//	x1 >= 0,
	//		//	x2 >= 0
	//	H = [2 0
	//		 0 8]; 
	//	f = [-8; -16];
	//  A = [1 1;1 0];
	//	b = [5;3];
	//	lb = [0; 0];
	//	ub = [%inf; %inf];
	//	[xopt,fopt,exitflag,output,lambda] = quadprogmat(H,f,A,b,[],[],lb,ub)
	// // Press ENTER to continue 
	//
	// Examples 
	//  //Find x in R^6 such that:
	//    Aeq= [1,-1,1,0,3,1;
	//         -1,0,-3,-4,5,6;
	//          2,5,3,0,1,0];
	//    beq=[1; 2; 3];
	//    A= [0,1,0,1,2,-1;
	//       -1,0,2,1,1,0];
	//    b = [-1; 2.5];
	//    lb=[-1000; -10000; 0; -1000; -1000; -1000];
	//    ub=[10000; 100; 1.5; 100; 100; 1000];
	//    x0 = repmat(0,6,1);
	//	  param = list("MaxIter", 300, "CpuTime", 100);
	//    //and minimize 0.5*x'*H*x + f'*x with
	//    f=[1; 2; 3; 4; 5; 6]; H=eye(6,6);
	//    [xopt,fopt,exitflag,output,lambda]=quadprogmat(H,f,A,b,Aeq,beq,lb,ub,x0,param)
	// Authors
	// Keyur Joshi, Saikiran, Iswarya, Harpreet Singh
    
    
	//To check the number of input and output argument
	[lhs , rhs] = argn();

	//To check the number of argument given by user
	if ( rhs < 2 | rhs == 3 | rhs == 5 | rhs == 7 | rhs > 10 ) then
		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be in the set of [2 4 6 8 9 10]"), "quadprogmat", rhs);
		error(errmsg)
	end

	H = [];
	f = [];
	A = [];
	b = [];
	Aeq = [];
	beq = []; 
	lb = [];
	ub = [];

	H = varargin(1);
	f = varargin(2);
	nbVar = size(H,1);

	if(nbVar == 0) then
		errmsg = msprintf(gettext("%s: Cannot determine the number of variables because input objective coefficients is empty"), "quadprogmat");
		error(errmsg);
	end

	if ( rhs<3 ) then
	  A = []
	  b = []
	else
	  A = varargin(3);
	  b = varargin(4);
	end
	  
	if ( rhs<5 ) then
	  Aeq = []
	  beq = []
	else
	  Aeq = varargin(5);
	  beq = varargin(6);
	end

	if ( rhs<7 ) then
		lb = repmat(-%inf,nbVar,1);
		ub = repmat(%inf,nbVar,1);
	else
		lb = varargin(7);
		ub = varargin(8);
	end

	if ( rhs<9 | size(varargin(9)) ==0 ) then
		x0 = repmat(0,nbVar,1)
	else
		x0 = varargin(9);
	end

	if ( rhs<10 | size(varargin(10)) ==0 ) then
		param = list();
	else
		param =varargin(10);
	end

	if (size(lb,2)==0) then
		lb = repmat(-%inf,nbVar,1);
	end

	if (size(ub,2)==0) then
		ub = repmat(%inf,nbVar,1);
	end

	if (size(f,2)==0) then
		f = repmat(0,nbVar,1);
	end

	if (modulo(size(param),2)) then
		errmsg = msprintf(gettext("%s: Size of parameters should be even"), "quadprogmat");
		error(errmsg);
	end

	//Check type of variables
	fot_Checktype("quadprogmat", H, "H", 1, "constant")
	fot_Checktype("quadprogmat", f, "f", 2, "constant")
	fot_Checktype("quadprogmat", A, "A", 3, "constant")
	fot_Checktype("quadprogmat", b, "b", 4, "constant")
	fot_Checktype("quadprogmat", Aeq, "Aeq", 5, "constant")
	fot_Checktype("quadprogmat", beq, "beq", 6, "constant")
	fot_Checktype("quadprogmat", lb, "lb", 7, "constant")
	fot_Checktype("quadprogmat", ub, "ub", 8, "constant")
	fot_Checktype("quadprogmat", x0, "x0", 9, "constant")
	fot_Checktype("quadprogmat", param, "param", 10, "list")

	options = list("MaxIter"     , [3000], ...
				  "CpuTime"   , [600]);

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

	nbConInEq = size(A,1);
	nbConEq = size(Aeq,1);

	// Check if the user gives row vector 
	// and Changing it to a column matrix

	if (size(f,2)== [nbVar]) then
		f=f';
	end

	if (size(lb,2)== [nbVar]) then
		lb = lb';
	end

	if (size(ub,2)== [nbVar]) then
		ub = ub';
	end

	if (size(b,2)==nbConInEq) then
		b = b';
	end

	if (size(beq,2)== nbConEq) then
		beq = beq';
	end

	if (size(x0,2)== [nbVar]) then
		x0=x0';
	end

	//Checking the H matrix which needs to be a symmetric matrix
	if ( ~isequal(H,H')) then
		errmsg = msprintf(gettext("%s: H is not a symmetric matrix"), "quadprogmat");
		error(errmsg);
	end

	//Check the size of f which should equal to the number of variable
	if ( size(f,1) ~= [nbVar]) then
		errmsg = msprintf(gettext("%s: The number of rows and columns in H must be equal the number of elements of f"), "quadprogmat");
		error(errmsg);
	end

	//Check the size of inequality constraint which should be equal to the number of variables
	if ( size(A,2) ~= nbVar & size(A,2) ~= 0) then
		errmsg = msprintf(gettext("%s: The number of columns in A must be the same as the number of elements of f"), "quadprogmat");
		error(errmsg);
	end

	//Check the size of equality constraint which should be equal to the number of variables
	if ( size(Aeq,2) ~= nbVar & size(Aeq,2) ~= 0 ) then
		errmsg = msprintf(gettext("%s: The number of columns in Aeq must be the same as the number of elements of f"), "quadprogmat");
		error(errmsg);
	end

	//Check the size of Lower Bound which should be equal to the number of variables
	if ( size(lb,1) ~= nbVar) then
		errmsg = msprintf(gettext("%s: The Lower Bound is not equal to the number of variables"), "quadprogmat");
		error(errmsg);
	end

	//Check the size of Upper Bound which should equal to the number of variables
	if ( size(ub,1) ~= nbVar) then
		errmsg = msprintf(gettext("%s: The Upper Bound is not equal to the number of variables"), "quadprogmat");
		error(errmsg);
	end

	//Check the size of constraints of Lower Bound which should equal to the number of constraints
	if ( size(b,1) ~= nbConInEq & size(b,1) ~= 0) then
		errmsg = msprintf(gettext("%s: The number of rows in A must be the same as the number of elements of b"), "quadprogmat");
		error(errmsg);
	end

	//Check the size of constraints of Upper Bound which should equal to the number of constraints
	if ( size(beq,1) ~= nbConEq & size(beq,1) ~= 0) then
		errmsg = msprintf(gettext("%s: The number of rows in Aeq must be the same as the number of elements of beq"), "quadprogmat");
		error(errmsg);
	end

	//Check the size of initial of variables which should equal to the number of variables
	if ( size(x0,1) ~= nbVar) then
		warnmsg = msprintf(gettext("%s: Ignoring initial guess of variables as it is not equal to the number of variables"), "quadprogmat");
		warning(warnmsg);
		x0 = repmat(0,nbVar,1);
	end

	//Check if the user gives a matrix instead of a vector

	if ((size(f,1)~=1)& (size(f,2)~=1)) then
		errmsg = msprintf(gettext("%s: f should be a vector"), "quadprogmat");
		error(errmsg); 
	end

	if (size(lb,1)~=1)& (size(ub,2)~=1) then
		errmsg = msprintf(gettext("%s: Lower Bound should be a vector"), "quadprogmat");
		error(errmsg); 
	end

	if (size(ub,1)~=1)& (size(ub,2)~=1) then
		errmsg = msprintf(gettext("%s: Upper Bound should be a vector"), "quadprogmat");
		error(errmsg); 
	end

	if (nbConInEq) then
		if ((size(b,1)~=1)& (size(b,2)~=1)) then
		    errmsg = msprintf(gettext("%s: Constraint Lower Bound should be a vector"), "quadprogmat");
		    error(errmsg); 
		end
	end

	if (nbConEq) then
		if (size(beq,1)~=1)& (size(beq,2)~=1) then
		    errmsg = msprintf(gettext("%s: Constraint should be a vector"), "quadprogmat");
		    error(errmsg); 
		end
	end

	for i = 1:nbConInEq
		if (b(i) == -%inf)
		   	errmsg = msprintf(gettext("%s: Value of b can not be negative infinity"), "quadprogmat");
		    error(errmsg); 
		end	
	end

	for i = 1:nbConEq
		if (beq(i) == -%inf)
		   	errmsg = msprintf(gettext("%s: Value of beq can not be negative infinity"), "quadprogmat");
		    error(errmsg); 
		end	
	end

	for i = 1:nbVar
		if(lb(i)>ub(i))
			errmsg = msprintf(gettext("%s: Problem has inconsistent variable bounds"), "lsqlin");
			error(errmsg);
		end
	end
   
	//Converting it into ipopt format
	f = f';
	lb = lb';
	ub = ub';
	x0 = x0';
	conMatrix = [Aeq;A];
	nbCon = size(conMatrix,1);
	conLB = [beq; repmat(-%inf,nbConInEq,1)]';
	conUB = [beq;b]' ; 

	[xopt,fopt,status,iter,Zl,Zu,lmbda] = solveqp(int32(nbVar),int32(nbCon),H,f,conMatrix,conLB,conUB,lb,ub,x0,options);
	xopt = xopt';
	exitflag = status;
	output = struct("Iterations"      , [], ..
					"ConstrViolation" ,[]);
	output.Iterations = iter;
	output.ConstrViolation = max([0;norm(Aeq*xopt-beq, 'inf');(lb'-xopt);(xopt-ub');(A*xopt-b)]);
	lambda = struct("lower"           , [], ..
		           "upper"           , [], ..
		           "eqlin"           , [], ..
				   "ineqlin"         , []);

	lambda.lower = Zl;
	lambda.upper = Zu;
	lambda.eqlin = lmbda(1:nbConEq);
	lambda.ineqlin = lmbda(nbConEq+1:nbCon);
   
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
