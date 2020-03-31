// Copyright (C) 2016 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Pranav Deshpande and Akshay Miterani
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

function varargout = intquadprog (varargin)
	// Solves a linear quadratic problem.
	//
	//   Calling Sequence
	//   xopt = intquadprog(H,f)
	//   xopt = intquadprog(H,f,intcon)
	//   xopt = intquadprog(H,f,intcon,A,b)
	//   xopt = intquadprog(H,f,intcon,A,b,Aeq,beq)
	//   xopt = intquadprog(H,f,intcon,A,b,Aeq,beq,lb,ub)
	//   xopt = intquadprog(H,f,intcon,A,b,Aeq,beq,lb,ub,x0)
	//   xopt = intquadprog(H,f,intcon,A,b,Aeq,beq,lb,ub,x0,options)
	//   xopt = intquadprog(H,f,intcon,A,b,Aeq,beq,lb,ub,x0,options,"file_path")
	//   [xopt,fopt,exitflag,output] = intquadprog( ... )
	//   
	//   Parameters
	//   H : a symmetric matrix of double, represents coefficients of quadratic in the quadratic problem.
	//   f : a vector of double, represents coefficients of linear in the quadratic problem
	//   intcon : a vector of integers, represents which variables are constrained to be integers
	//   A : a matrix of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b. 
	//   b : a vector of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b.
	//   Aeq : a matrix of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   beq : a vector of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   lb : a vector of double, contains lower bounds of the variables.
	//   ub : a vector of double, contains upper bounds of the variables.
	//   x0 : a vector of double, contains initial guess of variables.
	//   options : a list containing the parameters to be set.
	//	 file_path : path to bonmin opt file if used.
	//   xopt : a vector of double, the computed solution of the optimization problem.
	//   fopt : a double, the value of the function at x.
	//   exitflag : The exit status. See below for details.
	//   output : The structure consist of statistics about the optimization. See below for details.
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
	//	  & & x_i \in \!\, \mathbb{Z}, i \in \!\, intcon\\
	//    \end{eqnarray}
	//   </latex>
	//   
	//   The routine calls Bonmin for solving the quadratic problem, Bonmin is a library written in C++.
	//
	// The exitflag allows to know the status of the optimization which is given back by Bonmin.
	// <itemizedlist>
	//   <listitem>exitflag=0 : Optimal Solution Found. </listitem>
	//   <listitem>exitflag=1 : InFeasible Solution.</listitem>
	//   <listitem>exitflag=2 : Output is Continuous Unbounded.</listitem>
	//   <listitem>exitflag=3 : Limit Exceeded.</listitem>
	//   <listitem>exitflag=4 : User Interrupt.</listitem>
	//   <listitem>exitflag=5 : MINLP Error.</listitem>
	// </itemizedlist>
	// 
	// For more details on exitflag see the Bonmin page, go to http://www.coin-or.org/Bonmin
	//
	// The output data structure contains detailed informations about the optimization process. 
	// It has type "struct" and contains the following fields.
	// <itemizedlist>
	//   <listitem>output.constrviolation: The max-norm of the constraint violation.</listitem>
	// </itemizedlist>
	//
	// Examples
	// H = [1 -1; -1 2]; 
	// f = [-2; -6];
	// A = [1 1; -1 2; 2 1];
	// b = [2; 2; 3];
	// lb=[0,0];
	// ub=[%inf, %inf];
	// intcon = [1 2];
	//[xopt,fopt,status,output]=intquadprog(H,f,intcon,A,b,[],[],lb,ub)
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
	//	  intcon = [2 4];
	//    [xopt,fopt,exitflag,output]=intquadprog(H,f,intcon,A,b,Aeq,beq,lb,ub,x0,param)
	// Authors
	// Akshay Miterani and Pranav Deshpande

	//To check the number of input and output argument
	[lhs , rhs] = argn();

	//To check the number of argument given by user
	if ( rhs < 2 | rhs == 4 | rhs == 6 | rhs == 8 | rhs > 12 ) then
		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be in the set of [2 3 5 7 9 10 11 12]"), "intquadprog", rhs);
		error(errmsg)
	end
	
	//To check the number of output arguments
	if lhs > 4 then
		errmsg = msprintf(gettext("%s: Unexpected number of output arguments: %d provided while should be in the set of [1 2 3 4]"), "intquadprog", lhs);
	end

	H = [];
	f = [];
	intcon = [];
	A = [];
	b = [];
	Aeq = [];
	beq = []; 
	lb = [];
	ub = [];
	options = list();
	bonmin_options_file = '';

	H = varargin(1);
	f = varargin(2);
	nbVar = size(H,1);

	if(nbVar == 0) then
		errmsg = msprintf(gettext("%s: Cannot determine the number of variables because input objective coefficients is empty"), "intquadprog");
		error(errmsg);
	end
	
	if ( rhs>=3 ) then
		intcon=varargin(3);
	end
	
	if ( rhs<=3 ) then
	  A = []
	  b = []
	else
	  A = varargin(4);
	  b = varargin(5);
	end
	  
  if ( rhs<6 ) then
	  Aeq = []
	  beq = []
	else
	  Aeq = varargin(6);
	  beq = varargin(7);
	end

	if ( rhs<8 ) then
		lb = repmat(-%inf,nbVar,1);
		ub = repmat(%inf,nbVar,1);
	else
		lb = varargin(8);
		ub = varargin(9);
	end

	
	if ( rhs<10 | size(varargin(10)) ==0 ) then
		x0 = repmat(0,nbVar,1)
	  else
	    x0 = varargin(10);
	end
	
	if (rhs<11|size(varargin(11))==0) then
      options = list();
	else
      options = varargin(11);
	end
   
	if(rhs==12) then
		bonmin_options_file=varargin(12);
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

	//Check type of variables
	Checktype("intquadprog", H, "H", 1, "constant")
	Checktype("intquadprog", f, "f", 2, "constant")
	Checktype("intquadprog", intcon, "intcon", 3, "constant")
	Checktype("intquadprog", A, "A", 4, "constant")
	Checktype("intquadprog", b, "b", 5, "constant")
	Checktype("intquadprog", Aeq, "Aeq", 6, "constant")
	Checktype("intquadprog", beq, "beq", 7, "constant")
	Checktype("intquadprog", lb, "lb", 8, "constant")
	Checktype("intquadprog", ub, "ub", 9, "constant")
	Checktype("intquadprog", x0, "x0", 10, "constant")
	Checktype("intquadprog", options, "options", 11, "list")
	Checktype("intquadprog", bonmin_options_file, "bonmin_options_file", 12, "string")
	
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
		errmsg = msprintf(gettext("%s: H is not a symmetric matrix"), "intquadprog");
		error(errmsg);
	end

	//Check the size of f which should equal to the number of variable
	if ( size(f,1) ~= [nbVar]) then
		errmsg = msprintf(gettext("%s: The number of rows and columns in H must be equal the number of elements of f"), "intquadprog");
		error(errmsg);
	end
	
	//Error Checks for intcon
	 
	for i=1:size(intcon,2)
      if(intcon(i)>nbVar) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be less than the number of variables"), "intquadprog");
        error(errmsg);
      end

      if (intcon(i)<0) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be greater than 0 "), "intquadprog");
        error(errmsg);
      end

      if(modulo(intcon(i),1)) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be an integer "), "intquadprog");
        error(errmsg);
      end
	end

	//Check the size of inequality constraint which should be equal to the number of variables
	if ( size(A,2) ~= nbVar & size(A,2) ~= 0) then
		errmsg = msprintf(gettext("%s: The number of columns in A must be the same as the number of elements of f"), "intquadprog");
		error(errmsg);
	end

	//Check the size of equality constraint which should be equal to the number of variables
	if ( size(Aeq,2) ~= nbVar & size(Aeq,2) ~= 0 ) then
		errmsg = msprintf(gettext("%s: The number of columns in Aeq must be the same as the number of elements of f"), "intquadprog");
		error(errmsg);
	end

	//Check the size of Lower Bound which should be equal to the number of variables
	if ( size(lb,1) ~= nbVar) then
		errmsg = msprintf(gettext("%s: The Lower Bound is not equal to the number of variables"), "intquadprog");
		error(errmsg);
	end

	//Check the size of Upper Bound which should equal to the number of variables
	if ( size(ub,1) ~= nbVar) then
		errmsg = msprintf(gettext("%s: The Upper Bound is not equal to the number of variables"), "intquadprog");
		error(errmsg);
	end

	//Check the size of constraints of Lower Bound which should equal to the number of constraints
	if ( size(b,1) ~= nbConInEq & size(b,1) ~= 0) then
		errmsg = msprintf(gettext("%s: The number of rows in A must be the same as the number of elements of b"), "intquadprog");
		error(errmsg);
	end

	//Check the size of constraints of Upper Bound which should equal to the number of constraints
	if ( size(beq,1) ~= nbConEq & size(beq,1) ~= 0) then
		errmsg = msprintf(gettext("%s: The number of rows in Aeq must be the same as the number of elements of beq"), "intquadprog");
		error(errmsg);
	end

	//Check the size of initial of variables which should equal to the number of variables
	if ( size(x0,1) ~= nbVar) then
		warnmsg = msprintf(gettext("%s: Ignoring initial guess of variables as it is not equal to the number of variables"), "intquadprog");
		warning(warnmsg);
		x0 = repmat(0,nbVar,1);
	end

	//Check if the user gives a matrix instead of a vector

	if ((size(f,1)~=1)& (size(f,2)~=1)) then
		errmsg = msprintf(gettext("%s: f should be a vector"), "intquadprog");
		error(errmsg); 
	end

	if (size(lb,1)~=1)& (size(ub,2)~=1) then
		errmsg = msprintf(gettext("%s: Lower Bound should be a vector"), "intquadprog");
		error(errmsg); 
	end

	if (size(ub,1)~=1)& (size(ub,2)~=1) then
		errmsg = msprintf(gettext("%s: Upper Bound should be a vector"), "intquadprog");
		error(errmsg); 
	end

	if (nbConInEq) then
		if ((size(b,1)~=1)& (size(b,2)~=1)) then
		    errmsg = msprintf(gettext("%s: Constraint Lower Bound should be a vector"), "intquadprog");
		    error(errmsg); 
		end
	end

	if (nbConEq) then
		if (size(beq,1)~=1)& (size(beq,2)~=1) then
		    errmsg = msprintf(gettext("%s: Constraint should be a vector"), "intquadprog");
		    error(errmsg); 
		end
	end

	for i = 1:nbConInEq
		if (b(i) == -%inf) then
		   	errmsg = msprintf(gettext("%s: Value of b can not be negative infinity"), "intquadprog");
		    error(errmsg); 
		end	
	end

	for i = 1:nbConEq
		if (beq(i) == -%inf) then
		   	errmsg = msprintf(gettext("%s: Value of beq can not be negative infinity"), "intquadprog");
		    error(errmsg); 
		end	
	end

	for i = 1:nbVar
		if(lb(i)>ub(i)) then
			errmsg = msprintf(gettext("%s: Problem has inconsistent variable bounds"), "intquadprog");
			error(errmsg);
		end
	end
	
	// Checking if the specified options file exists or not

	
  	if (~isfile(bonmin_options_file) & ~(bonmin_options_file==  "") ) then
   		error(999, "The specified options file does not exist!");
 	end
	
 	
 	// Validating Options
 	if (type(options) ~= 15) then
		errmsg = msprintf(gettext("%s: Options should be a list "), "intlinprog");
		error(errmsg);
	end
   

	if (modulo(size(options),2)) then
		errmsg = msprintf(gettext("%s: Size of parameters should be even"), "intlinprog");
		error(errmsg);
	end
   
	//Pusing options as required to a double array
    optval = [];
    ifval =[];
    if length(options) == 0 then
        optval = [0 0 0 0 0];
        ifval = [0 0 0 0 0];
    else
        optval = [0 0 0 0 0];
        ifval = [0 0 0 0 0];
        for i=1:2:length(options)
            select convstr(options(i),'l')
                case 'integertolerance' then
					ifval(1) = 1;
                    optval(1) = options(i+1);
                case 'maxnodes' then
					ifval(2) = 1;
                    optval(2) = options(i+1);
                case 'cputime' then
					ifval(3) = 1;
                    optval(3) = options(i+1);
                case 'allowablegap' then
					ifval(4) = 1;
                    optval(4) = options(i+1);
                case 'maxiter' then
					ifval(5) = 1;
                    optval(5) = options(i+1);
                else
                    error(999, 'Unknown string argument passed.');
                end
        end
    end

	//Converting it into bonmin format
	f = f';
	lb = lb';
	ub = ub';
	x0 = x0';
	conMatrix = [Aeq;A];
	nbCon = size(conMatrix,1);
	conLB = [beq; repmat(-%inf,nbConInEq,1)]';
	conUB = [beq;b]';
	intcon = intcon'
	intconSize = length(intcon);
	xopt=[]
	fopt=[]
	status=[]
	[xopt,fopt,status] = solveintqp(int32(nbVar), int32(nbCon), int32(intconSize), H, f, intcon, conMatrix, conLB, conUB, lb, ub, x0, optval,ifval,bonmin_options_file);
	xopt = xopt';

	output.ConstrViolation = max([0;norm(Aeq*xopt-beq, 'inf');(lb'-xopt);(xopt-ub');(A*xopt-b)]);
	
	varargout(1) = xopt
	varargout(2) = fopt
	varargout(3) = status
	varargout(4) = output
	
	select status
    case 0 then
        printf("\nOptimal Solution Found.\n");
    case 1 then
        printf("\nInFeasible Solution.\n");
    case 2 then
        printf("\nOutput is Continuous Unbounded.s\n");
    case 3 then
        printf("\nStop at Tiny Step.\n");
    case 4 then
        printf("\nUser Interrupt.\n");
    case 5 then
        printf("\nMINLP Error.\n");
    else
        printf("\nInvalid status returned. Notify the Toolbox authors\n");
        break;
    end

endfunction
