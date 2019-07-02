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


function [xopt,fopt,exitflag,output,lambda] = quadprog (varargin)
	// Solves a linear quadratic problem.
	//
	//   Calling Sequence
	//   xopt = quadprog(nbVar,nbCon,H,f,lb,ub,A,conLB,conUB)
	//   xopt = quadprog(nbVar,nbCon,H,f,lb,ub,A,conLB,conUB,x0)
	//   xopt = quadprog(nbVar,nbCon,H,f,lb,ub,A,conLB,conUB,x0,param)
	//   [xopt,fopt,exitflag,output,lamda] = quadprog( ... )
	//   
	//   Parameters
	//   nbVar : a double, number of variables
	//   nbCon : a double, number of constraints
	//   H : a symmetric matrix of double, represents coefficients of quadratic in the quadratic problem.
	//   f : a vector of double, represents coefficients of linear in the quadratic problem
	//   lb : a vector of double, contains lower bounds of the variables.
	//   ub : a vector of double, contains upper bounds of the variables.
	//   A : a matrix of double, contains the constraint matrix conLB ≤ A⋅x ≤ conUB.
	//   conLB : a vector of double, contains lower bounds of the constraints conLB ≤ A⋅x ≤ conUB. 
	//   conUB : a vector of double, contains upper bounds of the constraints conLB ≤ A⋅x ≤ conUB. 
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
	//    & \text{subject to} & conLB \leq A⋅x \leq conUB \\
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
	//	conUB = [5;3];
	//	conLB = [-%inf; -%inf];
	//	lb = [0; 0];
	//	ub = [%inf; %inf];
	//	nbVar = 2;
	//	nbCon = 2;
	//	[xopt,fopt,exitflag,output,lambda] = quadprog(nbVar,nbCon,H,f,lb,ub,A,conLB,conUB)
	//  //Press ENTER to continue
	//    
	// Examples 
	//      //Find x in R^6 such that:
	//      A= [1,-1,1,0,3,1;
	//         -1,0,-3,-4,5,6;
	//          2,5,3,0,1,0
	//          0,1,0,1,2,-1;
	//         -1,0,2,1,1,0];
	//      conLB=[1;2;3;-%inf;-%inf];
	//      conUB = [1;2;3;-1;2.5];
	//      lb=[-1000;-10000; 0; -1000; -1000; -1000];
	//      ub=[10000; 100; 1.5; 100; 100; 1000];
	//      //and minimize 0.5*x'⋅H⋅x + f'⋅x with
	//      f=[1; 2; 3; 4; 5; 6]; H=eye(6,6);
	//      nbVar = 6;
	//      nbCon = 5;
	//      x0 = repmat(0,nbVar,1);
	//	  param = list("MaxIter", 300, "CpuTime", 100);
	//      [xopt,fopt,exitflag,output,lambda]=quadprog(nbVar,nbCon,H,f,lb,ub,A,conLB,conUB,x0,param)
	// Authors
	// Keyur Joshi, Saikiran, Iswarya, Harpreet Singh
    
    
	//To check the number of input and output argument
    [lhs , rhs] = argn();
	
	//To check the number of argument given by user
	if ( rhs < 9 | rhs > 11 ) then
		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be 9, 10 or 11"), "quadprog", rhs);
		error(errmsg)
	end

	nbVar = [];
	nbCon = [];
	H = [];
	f = [];
	A = [];
	conLB = [];
	conUB = [];
	lb = [];
	ub = [];
   
	nbVar = varargin(1);
	nbCon = varargin(2);
	H = varargin(3);
	f = varargin(4);
	lb = varargin(5);
	ub = varargin(6);
	A = varargin(7);
	conLB = varargin(8);
	conUB = varargin(9);

    if (size(lb,2)==0) then
        lb = repmat(-%inf,nbVar,1);
    end
    
    if (size(ub,2)==0) then
        ub = repmat(%inf,nbVar,1);
    end

    if (size(f,2)==0) then
        f = repmat(0,nbVar,1);
    end
    
	if ( rhs<10 | size(varargin(10)) ==0 ) then
		x0 = repmat(0,nbVar,1);
	else
		x0 = varargin(10);
	end
   
   if ( rhs<11 | size(varargin(11)) ==0 ) then
      param = list(); 
   else
      param =varargin(11);
   end
   
   if (type(param) ~= 15) then
      errmsg = msprintf(gettext("%s: param should be a list "), "quadprog");
      error(errmsg);
   end
   
	if (modulo(size(param),2)) then
		errmsg = msprintf(gettext("%s: Size of parameters should be even"), "quadprog");
		error(errmsg);
	end

	//Check type of variables
	Checktype("quadprog", nbVar, "nbVar", 1, "constant")
	Checktype("quadprog", nbCon, "nbCon", 2, "constant")
	Checktype("quadprog", H, "H", 3, "constant")
	Checktype("quadprog", f, "f", 4, "constant")
	Checktype("quadprog", lb, "lb", 5, "constant")
	Checktype("quadprog", ub, "lb", 6, "constant")
	Checktype("quadprog", A, "A", 7, "constant")
	Checktype("quadprog", conLB, "conlb", 8, "constant")
	Checktype("quadprog", conUB, "conub", 9, "constant")
	Checktype("quadprog", x0, "x0", 10, "constant")

   options = list(..
      "MaxIter"     , [3000], ...
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

	if (size(f,2)== [nbVar]) then
		f=f';
	end

	if (size(lb,2)== [nbVar]) then
		lb = lb';
	end

   if (size(ub,2)== [nbVar]) then
      ub = ub';
   end

   if (size(conUB,2)== [nbCon]) then
      conUB = conUB';
   end

   if (size(conLB,2)== [nbCon]) then
      conLB = conLB';
   end

   if (size(x0,2)== [nbVar]) then
	x0=x0';
   end

   //IPOpt wants it in row matrix form
   f = f';
   lb = lb';
   ub = ub';
   conLB = conLB';
   conUB = conUB';
   x0 = x0';
   
   //Checking the H matrix which needs to be a symmetric matrix
   if ( ~isequal(H,H') ) then
      errmsg = msprintf(gettext("%s: H is not a symmetric matrix"), "quadprog");
      error(errmsg);
   end

   //Check the size of H which should equal to the number of variable
   if ( size(H) ~= [nbVar nbVar]) then
      errmsg = msprintf(gettext("%s: The Size of H is not equal to the number of variables"), "quadprog");
      error(errmsg);
   end
   
   //Check the size of p which should equal to the number of variable
   if ( size(f,2) ~= [nbVar]) then
      errmsg = msprintf(gettext("%s: The Size of f is not equal to the number of variables"), "quadprog");
      error(errmsg);
   end
   
   if (nbCon) then
          //Check the size of constraint which should equal to the number of variables
       if ( size(A,2) ~= nbVar) then
          errmsg = msprintf(gettext("%s: The size of constraints is not equal to the number of variables"), "quadprog");
          error(errmsg);
       end
   end

   //Check the number of constraint
   if ( size(A,1) ~= nbCon) then
      errmsg = msprintf(gettext("%s: The size of constraint matrix is not equal to the number of constraint given i.e. %d"), "quadprog", nbCon);
      error(errmsg);
   end

   //Check the size of Lower Bound which should equal to the number of variables
   if ( size(lb,2) ~= nbVar) then
      errmsg = msprintf(gettext("%s: The size of Lower Bound is not equal to the number of variables"), "quadprog");
      error(errmsg);
   end

   //Check the size of Upper Bound which should equal to the number of variables
   if ( size(ub,2) ~= nbVar) then
      errmsg = msprintf(gettext("%s: The size of Upper Bound is not equal to the number of variables"), "quadprog");
      error(errmsg);
   end

   //Check the size of constraints of Lower Bound which should equal to the number of constraints
   if ( size(conLB,2) ~= nbCon) then
      errmsg = msprintf(gettext("%s: The size of Lower Bound of constraints is not equal to the number of constraints"), "quadprog");
      error(errmsg);
   end

   //Check the size of constraints of Upper Bound which should equal to the number of constraints
   if ( size(conUB,2) ~= nbCon) then
      errmsg = msprintf(gettext("%s: The size of Upper Bound of constraints is not equal to the number of constraints"), "quadprog");
      error(errmsg);
   end
    
   //Check the size of initial of variables which should equal to the number of variables
   if ( size(x0,2) ~= nbVar | size(x0,"*")>nbVar) then
      warnmsg = msprintf(gettext("%s: Ignoring initial guess of variables as it is not equal to the number of variables"), "quadprog");
      warning(warnmsg);
	  x0 = repmat(0,1,nbVar);
   end
   
   //Check if the user gives a matrix instead of a vector
   
   if ((size(f,1)~=1)& (size(f,2)~=1)) then
      errmsg = msprintf(gettext("%s: f should be a vector"), "quadprog");
      error(errmsg); 
   end
   
   if (size(lb,1)~=1)& (size(lb,2)~=1) then
      errmsg = msprintf(gettext("%s: Lower Bound should be a vector"), "quadprog");
      error(errmsg); 
   end
   
   if (size(ub,1)~=1)& (size(ub,2)~=1) then
      errmsg = msprintf(gettext("%s: Upper Bound should be a vector"), "quadprog");
      error(errmsg); 
   end
   
   if (nbCon) then
        if ((size(conLB,1)~=1)& (size(b,2)~=1)) then
            errmsg = msprintf(gettext("%s: Constraint Lower Bound should be a vector"), "quadprog");
            error(errmsg); 
        end

        if (size(conUB,1)~=1)& (size(beq,2)~=1) then
            errmsg = msprintf(gettext("%s: Constraint should be a vector"), "quadprog");
            error(errmsg); 
        end
   end
   
    // Check if the user gives infinity or negative infinity in conLB or conUB
	for i = 1:nbCon
		if (conLB(i) == %inf)
		   	errmsg = msprintf(gettext("%s: Value of Lower Bound can not be infinity"), "quadprog");
    		error(errmsg); 
  		end	

		if (conUB(i) == -%inf)
		   	errmsg = msprintf(gettext("%s: Value of Upper Bound can not be negative infinity"), "quadprog");
    		error(errmsg); 
		end	
	end

	for i = 1:nbVar
		if(lb(i)>ub(i))
			errmsg = msprintf(gettext("%s: Problem has inconsistent variable bounds"), "lsqlin");
			error(errmsg);
		end
	end

	[xopt,fopt,status,iter,Zl,Zu,lmbda] = solveqp(int32(nbVar), int32(nbCon),H,f,A,conLB,conUB,lb,ub,x0,options);
  
   xopt = xopt';
   exitflag = status;
	output = struct("Iterations"      , [], ..
					"ConstrViolation" ,[]);
	output.Iterations = iter;
	output.ConstrViolation = max([0;(conLB'-A*xopt);(A*xopt - conUB');(lb'-xopt);(xopt-ub')]);
   lambda = struct("lower"           , [], ..
                   "upper"           , [], ..
                   "constraint"      , []);
   
   lambda.lower = Zl;
   lambda.upper = Zu;
   lambda.constraint = lmbda;

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
