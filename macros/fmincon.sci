// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: R.Vidyadhar & Vignesh Kannan
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

function [xopt,fopt,exitflag,output,lambda,gradient,hessian] = fmincon (varargin)
  	// Solves a multi-variable constrainted optimization problem
  	//
  	//   Calling Sequence
  	//   xopt = fmincon(f,x0,A,b)
  	//   xopt = fmincon(f,x0,A,b,Aeq,beq)
  	//   xopt = fmincon(f,x0,A,b,Aeq,beq,lb,ub)
  	//   xopt = fmincon(f,x0,A,b,Aeq,beq,lb,ub,nlc)
  	//   xopt = fmincon(f,x0,A,b,Aeq,beq,lb,ub,nlc,options)
  	//   [xopt,fopt] = fmincon(.....)
  	//   [xopt,fopt,exitflag]= fmincon(.....)
  	//   [xopt,fopt,exitflag,output]= fmincon(.....)
  	//   [xopt,fopt,exitflag,output,lambda]=fmincon(.....)
  	//   [xopt,fopt,exitflag,output,lambda,gradient]=fmincon(.....)
  	//   [xopt,fopt,exitflag,output,lambda,gradient,hessian]=fmincon(.....)
  	//
  	//   Parameters
  	//   f : a function, representing the objective function of the problem 
  	//   x0 : a vector of doubles, containing the starting values of variables of size (1 X n) or (n X 1) where 'n' is the number of Variables
  	//   A : a matrix of doubles, containing the coefficients of linear inequality constraints of size (m X n) where 'm' is the number of linear inequality constraints
  	//   b : a vector of doubles, related to 'A' and containing the the Right hand side equation of the linear inequality constraints of size (m X 1)
  	//   Aeq : a matrix of doubles, containing the coefficients of linear equality constraints of size (m1 X n) where 'm1' is the number of linear equality constraints
  	//   beq : a vector of doubles, related to 'Aeq' and containing the the Right hand side equation of the linear equality constraints of size (m1 X 1)
  	//   lb : a vector of doubles, containing the lower bounds of the variables of size (1 X n) or (n X 1) where 'n' is the number of Variables
  	//   ub : a vector of doubles, containing the upper bounds of the variables of size (1 X n) or (n X 1) where 'n' is the number of Variables
  	//   nlc : a function, representing the Non-linear Constraints functions(both Equality and Inequality) of the problem. It is declared in such a way that non-linear inequality constraints are defined first as a single row vector (c), followed by non-linear equality constraints as another single row vector (ceq). Refer Example for definition of Constraint function.
  	//   options : a list, containing the option for user to specify. See below for details. 
  	//   xopt : a vector of doubles, cointating the computed solution of the optimization problem
  	//   fopt : a scalar of double, containing the the function value at x
  	//   exitflag : a scalar of integer, containing the flag which denotes the reason for termination of algorithm. See below for details.
  	//   output : a structure, containing the information about the optimization. See below for details.
  	//   lambda : a structure, containing the Lagrange multipliers of lower bound, upper bound and constraints at the optimized point. See below for details.
  	//   gradient : a vector of doubles, containing the Objective's gradient of the solution.
  	//   hessian  : a matrix of doubles, containing the Lagrangian's hessian of the solution.
  	//
  	//   Description
  	//   Search the minimum of a constrained optimization problem specified by :
  	//   Find the minimum of f(x) such that 
  	//
  	//   <latex>
  	//    \begin{eqnarray}
  	//    &\mbox{min}_{x}
  	//    & f(x) \\
  	//    & \text{subject to} & A*x \leq b \\
  	//    & & Aeq*x \ = beq\\
  	//	& & c(x) \leq  0\\
  	//    & & ceq(x) \ =  0\\
  	//    & & lb \leq x \leq ub \\
  	//    \end{eqnarray}
  	//   </latex>
  	//
  	//   The routine calls Ipopt for solving the Constrained Optimization problem, Ipopt is a library written in C++.
  	//
  	// The options allows the user to set various parameters of the Optimization problem. 
  	// It should be defined as type "list" and contains the following fields.
	// <itemizedlist>
	//   <listitem>Syntax : options= list("MaxIter", [---], "CpuTime", [---], "GradObj", ---, "Hessian", ---, "GradCon", ---);</listitem>
	//   <listitem>MaxIter : a Scalar, containing the Maximum Number of Iteration that the solver should take.</listitem>
	//   <listitem>CpuTime : a Scalar, containing the Maximum amount of CPU Time that the solver should take.</listitem>
	//   <listitem>GradObj : a function, representing the gradient function of the Objective in Vector Form.</listitem>
	//   <listitem>Hessian : a function, representing the hessian function of the Lagrange in Symmetric Matrix Form with Input parameters x, Objective factor and Lambda. Refer Example for definition of Lagrangian Hessian function.</listitem>
	//   <listitem>GradCon : a function, representing the gradient of the Non-Linear Constraints (both Equality and Inequality) of the problem. It is declared in such a way that gradient of non-linear inequality constraints are defined first as a separate Matrix (cg of size m2 X n or as an empty), followed by gradient of non-linear equality constraints as a separate Matrix (ceqg of size m2 X n or as an empty) where m2 & m3 are number of non-linear inequality and equality constraints respectively.</listitem>
	//   <listitem>Default Values : options = list("MaxIter", [3000], "CpuTime", [600]);</listitem>
	// </itemizedlist>
	//
	// The exitflag allows to know the status of the optimization which is given back by Ipopt.
	// <itemizedlist>
	//   <listitem>exitflag=0 : Optimal Solution Found </listitem>
	//   <listitem>exitflag=1 : Maximum Number of Iterations Exceeded. Output may not be optimal.</listitem>
	//   <listitem>exitflag=2 : Maximum amount of CPU Time exceeded. Output may not be optimal.</listitem>
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
	//   <listitem>output.Iterations: The number of iterations performed during the search</listitem>
	//   <listitem>output.Cpu_Time: The total cpu-time spend during the search</listitem>
	//   <listitem>output.Objective_Evaluation: The number of Objective Evaluations performed during the search</listitem>
	//   <listitem>output.Dual_Infeasibility: The Dual Infeasiblity of the final soution</listitem>
	//	 <listitem>output.Message: The output message for the problem</listitem>
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
	//   <listitem>lambda.eqnonlin: The Lagrange multipliers for the non-linear equality constraints.</listitem>
	//   <listitem>lambda.ineqnonlin: The Lagrange multipliers for the non-linear inequality constraints.</listitem>
	// </itemizedlist>
	//
  	// Examples
  	//	//Find x in R^2 such that it minimizes:
  	//    //f(x)= -x1 -x2/3
  	//    //x0=[0,0]
  	//    //constraint-1 (c1): x1 + x2 <= 2
  	//    //constraint-2 (c2): x1 + x2/4 <= 1 
  	//    //constraint-3 (c3): x1 - x2 <= 2
  	//    //constraint-4 (c4): -x1/4 - x2 <= 1
  	//    //constraint-5 (c5): -x1 - x2 <= -1
  	//    //constraint-6 (c6): -x1 + x2 <= 2
  	//    //constraint-7 (c7): x1 + x2 = 2  
  	//    //Objective function to be minimised
  	//  function y=f(x)
  	//		y=-x(1)-x(2)/3;
  	//	endfunction
  	//	//Starting point, linear constraints and variable bounds  
  	//	x0=[0 , 0];
  	//	A=[1,1 ; 1,1/4 ; 1,-1 ; -1/4,-1 ; -1,-1 ; -1,1];
  	//	b=[2;1;2;1;-1;2];
  	//	Aeq=[1,1];
  	// 	beq=[2];
  	//	lb=[];
  	//	ub=[];
  	//    nlc=[];
  	//	//Gradient of objective function
  	//	function y= fGrad(x)
  	//		y= [-1,-1/3];
  	// 	endfunction
  	//	//Hessian of lagrangian
  	//	function y= lHess(x,obj,lambda)
  	//		y= obj*[0,0;0,0] 
  	//  endfunction
  	//	//Options
  	//	options=list("GradObj", fGrad, "Hessian", lHess);
  	//    //Calling Ipopt
  	//	[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options)
	// // Press ENTER to continue
  	//
  	// Examples
  	//	//Find x in R^3 such that it minimizes:
  	//    //f(x)= x1*x2 + x2*x3
  	//    //x0=[0.1 , 0.1 , 0.1]
  	//    //constraint-1 (c1): x1^2 - x2^2 + x3^2 <= 2
  	//    //constraint-2 (c2): x1^2 + x2^2 + x3^2 <= 10  
  	//    //Objective function to be minimised
  	//    function y=f(x)
  	//		y=x(1)*x(2)+x(2)*x(3);
  	//	endfunction
  	//	//Starting point, linear constraints and variable bounds  
  	//	x0=[0.1 , 0.1 , 0.1];
  	//	A=[];
  	//	b=[];
  	//	Aeq=[];
  	// 	beq=[];
  	//	lb=[];
  	//	ub=[];
  	//	//Nonlinear constraints  
  	//	function [c,ceq]=nlc(x)
  	//		c = [x(1)^2 - x(2)^2 + x(3)^2 - 2 , x(1)^2 + x(2)^2 + x(3)^2 - 10];
  	//  		ceq = [];
  	//	endfunction
  	//	//Gradient of objective function
  	//	function y= fGrad(x)
  	//		y= [x(2),x(1)+x(3),x(2)];
  	//	endfunction
  	//    //Hessian of the Lagrange Function
  	//	function y= lHess(x,obj,lambda)
  	//		y= obj*[0,1,0;1,0,1;0,1,0] + lambda(1)*[2,0,0;0,-2,0;0,0,2] + lambda(2)*[2,0,0;0,2,0;0,0,2]
  	//	endfunction
  	//    //Gradient of Non-Linear Constraints
  	//	function [cg,ceqg] = cGrad(x)
  	//		cg=[2*x(1) , -2*x(2) , 2*x(3) ; 2*x(1) , 2*x(2) , 2*x(3)];
  	//		ceqg=[];
  	//	endfunction
  	//	//Options  
  	//	options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", fGrad, "Hessian", lHess,"GradCon", cGrad);
  	//    //Calling Ipopt
  	//	[x,fval,exitflag,output] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options)
	// // Press ENTER to continue
  	//
  	// Examples
  	//    //The below problem is an unbounded problem:
  	//	//Find x in R^3 such that it minimizes:
  	//    //f(x)= -(x1^2 + x2^2 + x3^2)
  	//    //x0=[0.1 , 0.1 , 0.1]
  	//    //  x1 <= 0
  	//    //  x2 <= 0
  	//    //  x3 <= 0
  	//    //Objective function to be minimised
  	//    function y=f(x)
  	//		y=-(x(1)^2+x(2)^2+x(3)^2);
  	//	endfunction
  	//	//Starting point, linear constraints and variable bounds  
  	//	x0=[0.1 , 0.1 , 0.1];
  	//	A=[];
  	//	b=[];
  	//	Aeq=[];
  	// 	beq=[];
  	//	lb=[];
  	//	ub=[0,0,0];
  	//	//Options
  	//	options=list("MaxIter", [1500], "CpuTime", [500]);
  	//    //Calling Ipopt
  	//	[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,[],options)
	// // Press ENTER to continue
  	//
  	// Examples
  	//    //The below problem is an infeasible problem:
  	//	//Find x in R^3 such that in minimizes:
  	//    //f(x)=x1*x2 + x2*x3
  	//    //x0=[1,1,1]
  	//    //constraint-1 (c1): x1^2 <= 1
  	//    //constraint-2 (c2): x1^2 + x2^2 <= 1    
  	//    //constraint-3 (c3): x3^2 <= 1  
  	//    //constraint-4 (c4): x1^3 = 0.5  
  	//    //constraint-5 (c5): x2^2 + x3^2 = 0.75
  	//    // 0 <= x1 <=0.6
  	//    // 0.2 <= x2 <= inf
  	//    // -inf <= x3 <= 1
  	//    //Objective function to be minimised
  	//    function y=f(x)
  	//		y=x(1)*x(2)+x(2)*x(3);
  	//	endfunction
  	//	//Starting point, linear constraints and variable bounds  
  	//	x0=[1,1,1];
  	//	A=[];
  	//	b=[];
  	//	Aeq=[];
  	// 	beq=[];
  	//	lb=[0 0.2,-%inf];
  	//	ub=[0.6 %inf,1];
  	//	//Nonlinear constraints  
  	//	function [c,ceq]=nlc(x)
  	//		c=[x(1)^2-1,x(1)^2+x(2)^2-1,x(3)^2-1];
  	//		ceq=[x(1)^3-0.5,x(2)^2+x(3)^2-0.75];
  	//	endfunction
  	//	//Gradient of objective function
  	//	function y= fGrad(x)
  	//		y= [x(2),x(1)+x(3),x(2)];
  	//	endfunction
  	//    //Hessian of the Lagrange Function
  	//	function y= lHess(x,obj,lambda)
  	//		y= obj*[0,1,0;1,0,1;0,1,0] + lambda(1)*[2,0,0;0,0,0;0,0,0] + ..
	//				lambda(2)*[2,0,0;0,2,0;0,0,0] +lambda(3)*[0,0,0;0,0,0;0,0,2] + ..
	//				lambda(4)*[6*x(1),0,0;0,0,0;0,0,0] + lambda(5)*[0,0,0;0,2,0;0,0,2];
  	//	endfunction
  	//    //Gradient of Non-Linear Constraints
  	//	function [cg,ceqg] = cGrad(x)
  	//		cg = [2*x(1),0,0;2*x(1),2*x(2),0;0,0,2*x(3)];
  	//		ceqg = [3*x(1)^2,0,0;0,2*x(2),2*x(3)];
  	//	endfunction
   	//	//Options  
  	//	options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", fGrad, "Hessian", lHess,"GradCon", cGrad);
  	//    //Calling Ipopt
  	//	[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options)
 	// // Press ENTER to continue
  	// Authors
  	// 	R.Vidyadhar , Vignesh Kannan
 	
	
	//To check the number of input and output arguments
   	[lhs , rhs] = argn();
	
	//To check the number of arguments given by the user
   	if ( rhs<4 | rhs>10 ) then
    		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while it should be 4,6,8,9,10"), "fmincon", rhs);
    		error(errmsg)
   	end
    	
	if (rhs==5 | rhs==7) then
    	errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while it should be 4,6,8,9,10s"), "fmincon", rhs);
    	error(errmsg)
   	end
 
	//Storing the Input Parameters  
   	fun    	 = varargin(1);
   	x0   	 = varargin(2);
   	A    	 = varargin(3);
   	b    	 = varargin(4);
   	Aeq  	 = [];
   	beq  	 = [];
   	lb       = [];
   	ub       = [];
   	nlc      = [];
   	
   	if (rhs>4) then
   		Aeq  	 = varargin(5);
   		beq  	 = varargin(6);
   	end

   	if (rhs>6) then
   		lb       = varargin(7);
   		ub       = varargin(8);
   	end

   	if (rhs>8) then
   		nlc      = varargin(9);
	end
	
	//To check whether the 1st Input argument (fun) is a function or not
    Checktype("fmincon", fun, "f", 1, "function");
   
	//To check whether the 2nd Input argument (x0) is a vector/scalar
  	Checktype("fmincon", x0, "x0", 2, "constant");
  	
  	//To check and convert the 2nd Input argument (x0) to a row vector 
   	if((size(x0,1)~=1) & (size(x0,2)~=1)) then
   		errmsg = msprintf(gettext("%s: Expected Row Vector or Column Vector for x0 (Starting Point) or Starting Point cannot be Empty"), "fmincon");
   		error(errmsg);
    end

   	if(size(x0,2)==1) then
   		x0=x0';		//Converting x0 to a row vector, if it is a column vector
   	else 
   	 	x0=x0;		//Retaining the same, if it is already a row vector
   	end   	 	
    	s=size(x0);
  	
  	//To check the match between fun (1st Parameter) and x0 (2nd Parameter)
   	if(execstr('init=fun(x0)','errcatch')==21) then
		errmsg = msprintf(gettext("%s: Objective function and x0 did not match"), "fmincon");
   		error(errmsg);
	end
	
	//Converting the User defined Objective function into Required form (Error Detectable)
   	function [y,check] = f(x)
   		if(execstr('y=fun(x)','errcatch')==32 | execstr('y=fun(x)','errcatch')==27)
			y=0;
			check=1;
		else
			y=fun(x);
			if (isreal(y)==%F) then
				y=0;
				check=1;
  			else
				check=0;
			end
		end
	endfunction
   	
  	//To check whether the 3rd Input argument (A) is a Matrix/Vector
  	Checktype("fmincon", A, "A", 3, "constant");

	//To check for correct size of A(3rd paramter)
   	if(size(A,2)~=s(2) & size(A,2)~=0) then
   		errmsg = msprintf(gettext("%s: Expected Matrix of size (No of linear inequality constraints X No of Variables) or an Empty Matrix for Linear Inequality Constraint coefficient Matrix A"), "fmincon");
   		error(errmsg);
   	end

   	s1=size(A);
   	
	//To check whether the 4th Input argument (b) is a vector/scalar
  	Checktype("fmincon", b, "b", 4, "constant");
  	
  	//To check for the correct size of b (4th paramter) and convert it into a column vector which is required for Ipopt
    if(s1(2)==0) then
    	if(size(b,2)~=0) then
    		errmsg = msprintf(gettext("%s: As Linear Inequality Constraint coefficient Matrix A (3rd parameter) is empty, b (4th Parameter) should also be empty"), "fmincon");
   			error(errmsg);
   		end
	else
   		if((size(b,1)~=1) & (size(b,2)~=1)) then
   			errmsg = msprintf(gettext("%s: Expected Non empty Row/Column Vector for b (4th Parameter) for your Inputs "), "fmincon");
   			error(errmsg);
   		elseif(size(b,1)~=s1(1) & size(b,2)==1) then
   			errmsg = msprintf(gettext("%s: Expected Column Vector (number of linear inequality constraints X 1) for b (4th Parameter) "), "fmincon");
   			error(errmsg);
   		elseif(size(b,1)==s1(1) & size(b,2)==1) then 
   	 		b=b;
   		elseif(size(b,1)==1 & size(b,2)~=s1(1)) then
   			errmsg = msprintf(gettext("%s: Expected Row Vector (1 X number of linear inequality constraints) for b (4th Parameter) "), "fmincon");
   			error(errmsg);
   		elseif(size(b,1)==1 & size(b,2)==s1(1)) then
   			b=b';
   		end 
   	end
  	
  	//To check whether the 5th Input argument (Aeq) is a matrix/vector
  	Checktype("fmincon", Aeq, "Aeq", 5, "constant");
  	
  	//To check for the correct size of Aeq (5th paramter)
   	if(size(Aeq,2)~=s(2) & size(Aeq,2)~=0) then
   		errmsg = msprintf(gettext("%s: Expected Matrix of size (No of linear equality constraints X No of Variables) or an Empty Matrix for Linear Equality Constraint coefficient Matrix Aeq"), "fmincon");
   		error(errmsg);
   	end

   	s2=size(Aeq);

	//To check whether the 6th Input argument(beq) is a vector/scalar
	Checktype("fmincon", beq, "beq", 6, "constant");
  	
  	//To check for the correct size of beq(6th paramter) and convert it into a column vector which is required for Ipopt
    if(s2(2)==0) then
    	if(size(beq,2)~=0) then
    		errmsg = msprintf(gettext("%s: As Linear Equality Constraint coefficient Matrix Aeq (5th parameter) is empty, beq (6th Parameter) should also be empty"), "fmincon");
   			error(errmsg);
   		end
   	else
   		if((size(beq,1)~=1) & (size(beq,2)~=1)) then
   			errmsg = msprintf(gettext("%s: Expected Non empty Row/Column Vector for beq (6th Parameter)"), "fmincon");
   			error(errmsg);
   		elseif(size(beq,1)~=s2(1) & size(beq,2)==1) then
   			errmsg = msprintf(gettext("%s: Expected Column Vector (number of linear equality constraints X 1) for beq (6th Parameter) "), "fmincon");
   			error(errmsg);
   		elseif(size(beq,1)==s2(1) & size(beq,2)==1) then 
   	 		beq=beq;
   		elseif(size(beq,1)==1 & size(beq,2)~=s2(1)) then
   			errmsg = msprintf(gettext("%s: Expected Row Vector (1 X number of linear equality constraints) for beq (6th Parameter) "), "fmincon");
   			error(errmsg);
   		elseif(size(beq,1)==1 & size(beq,2)==s2(1)) then
   			beq=beq';
   		end 
   	end
   	
  	
  	//To check whether the 7th Input argument (lb) is a vector/scalar
	Checktype("fmincon", lb, "lb", 7, "constant");
  	
  	//To check for the correct size and data of lb (7th paramter) and convert it into a column vector as required by Ipopt
   	if (size(lb,2)==0) then
        	lb = repmat(-%inf,1,s(2));
    end
    
   	if (size(lb,1)~=1) & (size(lb,2)~=1) then
      errmsg = msprintf(gettext("%s: Lower Bound (7th Parameter) should be a vector"), "fmincon");
      error(errmsg); 
    elseif(size(lb,1)~=s(2) & size(lb,2)==1) then
   		errmsg = msprintf(gettext("%s: Expected Column Vector (number of Variables X 1) for lower bound (7th Parameter) "), "fmincon");
   		error(errmsg);
   	elseif(size(lb,1)==s(2) & size(lb,2)==1) then
   	 	lb=lb;
   	elseif(size(lb,1)==1 & size(lb,2)~=s(2)) then
   		errmsg = msprintf(gettext("%s: Expected Row Vector (1 X number of Variables) for lower bound (7th Parameter) "), "fmincon");
   		error(errmsg);
   	elseif(size(lb,1)==1 & size(lb,2)==s(2)) then
   		lb=lb';
   	end 
   	
   	//To check whether the 8th Input argument (ub) is a vector/scalar
   	Checktype("fmincon", ub, "ub", 8, "constant");
   	
   	//To check for the correct size and data of ub (8th paramter) and convert it into a column vector as required by Ipopt
    if (size(ub,2)==0) then
        ub = repmat(%inf,1,s(2));
    end
    
    if (size(ub,1)~=1)& (size(ub,2)~=1) then
      errmsg = msprintf(gettext("%s: Upper Bound (8th Parameter) should be a vector"), "fmincon");
      error(errmsg); 
    elseif(size(ub,1)~=s(2) & size(ub,2)==1) then
   		errmsg = msprintf(gettext("%s: Expected Column Vector (number of Variables X 1) for upper bound (8th Parameter) "), "fmincon");
   		error(errmsg);
   	elseif(size(ub,1)==s(2) & size(ub,2)==1) then
   	 	ub=ub;
   	elseif(size(ub,1)==1 & size(ub,2)~=s(2)) then
   		errmsg = msprintf(gettext("%s: Expected Row Vector (1 X number of Variables) for upper bound (8th Parameter) "), "fmincon");
   		error(errmsg);
   	elseif(size(ub,1)==1 & size(ub,2)==s(2)) then
   		ub=ub';
   	end 
    
    //To check the contents of lb & ub (7th & 8th Parameter)
    for i = 1:s(2)
		if (lb(i) == %inf) then
		   	errmsg = msprintf(gettext("%s: Value of Lower Bound can not be infinity"), "fmincon");
    		error(errmsg); 
  		end	

		if (ub(i) == -%inf) then
		   	errmsg = msprintf(gettext("%s: Value of Upper Bound can not be negative infinity"), "fmincon");
    		error(errmsg); 
		end	

		if(ub(i)-lb(i)<=1e-6) then
			errmsg = msprintf(gettext("%s: Difference between Upper Bound and Lower bound should be atleast > 10^-6 for variable number= %d "), "fmincon", i);
    		error(errmsg)
    	end
	end
	
	//To check whether the 10th Input argument (nlc) is a function or an empty matrix
	if (type(nlc) == 1 & size(nlc,2)==0 ) then
  		addnlc=[];
  		addnlc1=[];
  		no_nlc=0;
  		no_nlic=0;
		no_nlec=0;
  	elseif (type(nlc) == 13 | type(nlc) == 11) then
  		
  		if(execstr('[sample_c,sample_ceq] = nlc(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Non-Linear Constraint function(9th Parameter) and x0(2nd Parameter) did not match"), "fmincon");
   			error(errmsg);
		end
  		[sample_c,sample_ceq] = nlc(x0);
  		
  		if (size(sample_c,1)~=1 & size(sample_c,1)~=0) then
  			errmsg = msprintf(gettext("%s: Definition of c in Non-Linear Constraint function(9th Parameter) should be in the form of Row Vector or Empty Vector"), "fmincon");
    		error(errmsg)
    	end
  		
  		if (size(sample_ceq,1)~=1 & size(sample_ceq,1)~=0) then
  			errmsg = msprintf(gettext("%s: Definition of ceq in Non-Linear Constraint function(9th Parameter) should be in the form of Row Vector or Empty Vector"), "fmincon");
    		error(errmsg)
    	end
    	
  		no_nlic = size(sample_c,2);
  		no_nlec = size(sample_ceq,2);
  		no_nlc = no_nlic + no_nlec;
  		//Constructing a single output variable function for nlc
  		function y = addnlc(x)
  			[c,ceq] = nlc(x);
  			y = [c';ceq'];
  		endfunction
  		
  		//To check the addnlc function
  		if(execstr('sample_allcon = addnlc(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Non-Linear Constraint function(9th Parameter) and x0(2nd Parameter) did not match"), "fmincon");
   			error(errmsg);
		end
  		sample_allcon = addnlc(x0);
  		
  		if (size(sample_allcon,1)==0 & size(sample_allcon,2)==0) then
  		
  		elseif (size(sample_allcon,1)~=no_nlc | size(sample_allcon,2)~=1) then
  			errmsg = msprintf(gettext("%s: Please check the Non-Linear Constraint function (9th Parameter) function"), "fmincon");
    		error(errmsg)
    	end
    	
    	//Constructing a nlc function with error deduction
  		function [y,check] = addnlc1(x)
  			if(execstr('y = addnlc(x)','errcatch')==32 | execstr('y = addnlc(x)','errcatch')==27)
				y = 0;
				check=1;
			else
  				y = addnlc(x);
				if (isreal(y)==%F) then
					y = 0;
					check=1;
  				else
					check=0;
				end
			end
  		endfunction
    	
    else 
    	errmsg = msprintf(gettext("%s: Non Linear Constraint (9th Parameter) should be a function or an Empty Matrix"), "fmincon");
    	error(errmsg)
    end
   	
   	//To check whether options has been entered by the user   
   	if ( rhs<10 ) then
      		param = list();
       else
      		param =varargin(10); //Storing the 3rd Input Parameter in an intermediate list named 'param'
    end
   
	//If options has been entered, then check its type for 'list'   
   Checktype("fmincon", param, "options", 10, "list");
   
	//If options has been entered, then check whether an even number of entires has been entered   
   	if (modulo(size(param),2)) then
		errmsg = msprintf(gettext("%s: Size of Options (list) should be even"), "fmincon");
		error(errmsg);
   	end
   	
   	//To set default values for options, if user doesn't enter options
	options = list("MaxIter", [3000], "CpuTime", [600]);

	//Flags to check whether Gradient is "ON"/"OFF" and Hessian is "ON"/"OFF" 
   	flag1=0;
   	flag2=0;
   	flag3=0;
   	
   	//Function for Gradient and Hessian
   	fGrad=[];
   	fGrad1=[];
   	lHess=[];
   	lHess1=[];
   	cGrad=[];
	addcGrad=[];
	addcGrad1=[];
   	
	//To check the user entry for options and storing it
   	for i = 1:(size(param))/2
       	select convstr(param(2*i-1),'l')
          	case "maxiter" then
          			if (type(param(2*i))~=1) then
          				errmsg = msprintf(gettext("%s: Value for Maximum Iteration should be a Constant"), "fmincon");
    	      			error(errmsg);
          			else
          				options(2) = param(2*i);    //Setting the maximum number of iterations as per user entry
          			end
       		case "cputime" then
          			if (type(param(2*i))~=1) then
          				errmsg = msprintf(gettext("%s: Value for Maximum Cpu-time should be a Constant"), "fmincon");
    	      			error(errmsg);
          			else
          				options(4) = param(2*i);    //Setting the maximum CPU time as per user entry
          			end
        	case "gradobj" then
        			if (type(param(2*i))==10) then
        				if (convstr(param(2*i))=="off") then
        					flag1 =0;
        				else
        					errmsg = msprintf(gettext("%s: Unrecognized String [%s] entered for the option- %s."), "fmincon",param(2*i), param(2*i-1));
    	      				error(errmsg);
        				end
        			else
						flag1 = 1;
        				fGrad = param(2*i);        				      
        			end
        	case "hessian" then
        			if (type(param(2*i))==10) then
        				if (convstr(param(2*i))=="off") then
        					flag2 =0;
        				else
        					errmsg = msprintf(gettext("%s: Unrecognized String [%s] entered for the option- %s."), "fmincon",param(2*i), param(2*i-1));
    	      				error(errmsg);
        				end
        			else
						flag2 = 1;
        				lHess = param(2*i);        				      
        			end
        	case "gradcon" then
        			if (type(param(2*i))==10) then
        				if (convstr(param(2*i))=="off") then
        					flag3 =0;
        				else
        					errmsg = msprintf(gettext("%s: Unrecognized String [%s] entered for the option- %s."), "fmincon",param(2*i), param(2*i-1));
    	      				error(errmsg);
        				end
        			else
						flag3 = 1;
        				cGrad = param(2*i);        				      
        			end
       		else
    	     	 	errmsg = msprintf(gettext("%s: Unrecognized parameter name %s."), "fmincon", param(2*i-1));
    	      		error(errmsg);
        end        					
     end 	       				      
	
   //To check for correct input of Objective Gradient function from the user	     	
   if (flag1==1) then
   		Checktype("fmincon", fGrad, "fGrad", 10, "function");
   		
   		if(execstr('sample_fGrad=fGrad(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Gradient function of Objective and x0 did not match "), "fmincon");
   			error(errmsg);
		end
		
		sample_fGrad=fGrad(x0);
		
		if (size(sample_fGrad,1)==s(2) & size(sample_fGrad,2)==1) then
		elseif (size(sample_fGrad,1)==1 & size(sample_fGrad,2)==s(2)) then
		else
   			errmsg = msprintf(gettext("%s: Wrong Input for Objective Gradient function(3rd Parameter)---->Row Vector function of size [1 X %d] is Expected"), "fmincon",s(2));
   			error(errmsg);
   		end  		
   end
   
   //To check for correct input of Lagrangian Hessian function from the user
   if (flag2==1) then
   		Checktype("fmincon", lHess, "lHess", 10, "function");
   		
   		if(execstr('sample_lHess=lHess(x0,1,1:no_nlc)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Hessian function of Objective and x0 did not match "), "fmincon");
   			error(errmsg);
		end
		sample_lHess=lHess(x0,1,1:no_nlc);
   		if(size(sample_lHess,1)~=s(2) | size(sample_lHess,2)~=s(2)) then
   			errmsg = msprintf(gettext("%s: Wrong Input for Objective Hessian function(10th Parameter)---->Symmetric Matrix function is Expected "), "fmincon");
   			error(errmsg);
   		end 		
   	end
   	
   	//To check for correct input of Constraint Gradient function from the user
   	if (flag3==1) then
   		Checktype("fmincon", cGrad, "cGrad", 10, "function");
   		
   		if(execstr('[sample_cGrad,sample_ceqg]=cGrad(x0)','errcatch')==21)
			errmsg = msprintf(gettext("%s: Gradient function of Constraint and x0 did not match "), "fmincon");
   			error(errmsg);
		end
		[sample_cGrad,sample_ceqg]=cGrad(x0);
		
		if  (size(sample_cGrad,2)==0) then 		
		elseif (size(sample_cGrad,1)~=no_nlic | size(sample_cGrad,2)~=s(2)) then
   			errmsg = msprintf(gettext("%s: Definition of (cGrad) in Non-Linear Constraint function(10th Parameter) should be in the form of (m X n) or Empty Matrix where m is number of Non- linear inequality constraints and n is number of Variables"), "fmincon");
   			error(errmsg);		
		end
		
		if (size(sample_ceqg,2)==0) then		
		elseif (size(sample_ceqg,1)~=no_nlec | size(sample_ceqg,2)~=s(2)) then
   			errmsg = msprintf(gettext("%s: Definition of (ceqg) in Non-Linear Constraint function(10th Parameter) should be in the form of (m X n) or Empty Matrix where m is number of Non- linear equality constraints and n is number of Variables"), "fmincon");
   			error(errmsg);
		end

		function y = addcGrad(x)
  			[sample_cGrad,sample_ceqg] = cGrad(x);
  			y = [sample_cGrad;sample_ceqg];
  		endfunction
   	
   		sample_addcGrad=addcGrad(x0);
		if(size(sample_addcGrad,1)~=no_nlc | size(sample_addcGrad,2)~=s(2)) then
   			errmsg = msprintf(gettext("%s:  Wrong Input for Constraint Gradient function(10th Parameter) (Refer Help)"), "fmincon");
   			error(errmsg);
   		end
   	end
   	
   	//Defining an inbuilt Objective gradient function 
   	function [y,check] = fGrad1(x) 
   			if flag1==1 then
   			    if(execstr('y=fGrad(x)','errcatch')==32 | execstr('y=fGrad(x)','errcatch')==27)
			    	y = 0;
			    	check=1;
			    else
			    	y=fGrad(x);
			    	if (isreal(y)==%F) then
			    		y = 0;
			    		check=1;
  			    	else
			    		check=0;
			    	end
			    end
		    else
		        if(execstr('y=numderivative(fun,x)','errcatch')==10000)
			    	y=0;
			    	check=1;
			    else
			    	y=numderivative(fun,x);
			    	if (isreal(y)==%F) then
			        	y=0;
			    		check=1;
  			    	else
			    		check=0;
			    	end
			    end		
		    end
  	endfunction
  	
  	//Defining an inbuilt Lagrangian Hessian function 
  	function [y,check] = lHess1(x,obj,lambda)
   			if flag2==1 then
   				if(execstr('y=lHess(x,obj,lambda)','errcatch')==32 | execstr('y=lHess(x,obj,lambda)','errcatch')==27)
					y = 0;
					check=1;
				else
					y=lHess(x,obj,lambda);
					if (isreal(y)==%F) then
						y = 0;
						check=1;
  					else
						check=0;
					end
				end
			else
				if(execstr('[grad,y1]=numderivative(fun,x)','errcatch')==10000)
					check1=1;
				else
					[grad,y1]=numderivative(fun,x);
					if (isreal(y1)==%F) then
						check1=1;
  					else
						check1=0;
					end	
				end	
				
				if check1==0 then
					if no_nlc~=0 then
						if(execstr('[grad,y2]=numderivative(addnlc,x)','errcatch')==10000)
							check2=1;
						else
							[grad,y2]=numderivative(addnlc,x);
							if (isreal(y2)==%F) then
								check2=1;
  							else
								check2=0;
							end	
						end
						if check2==0 then
							y2=matrix(y2, no_nlc*s(2)*s(2),1)
							for i = 1:s(2)*s(2)
								y(i)=0;
								for j = 1:no_nlc
									y(i)= y(i) + lambda(j)*y2((i-1)*no_nlc+j);
								end
							end
					
							for i = 1:s(2)*s(2)
								y(i) = y(i)+ obj*y1(i);
							end
							check=0;
						else
							check=1;
						end
					else
						check=0;
						for i = 1:s(2)*s(2)
							y(i) = obj*y1(i);
						end
					end
				else
					check=1;
					y=[0];	
				end				
			end
		disp (y)
		disp(check)
  	endfunction

   	//Defining an inbuilt Constraint gradient function 
   	function [y,check] = addcGrad1(x)
   			if flag3==1 then
   			    if(execstr('y=addcGrad(x)','errcatch')==32 | execstr('y=addcGrad(x)','errcatch')==27)
			    	y = 0;
			    	check=1;
			    else
			    	y=addcGrad(x);
			    	if (isreal(y)==%F) then
			    		y = 0;
			    		check=1;
  			    	else
			    		check=0;
			    	end
			    end
			else
			    if(execstr('y=numderivative(addnlc,x)','errcatch')==10000)
			    	y=0;
			    	check=1;
			    else
			    	y=numderivative(addnlc,x);
			    	if (isreal(y)==%F) then
			    		y=0;
			    		check=1;
  			    	else
			    		check=0;
			       	end
			    end
			end	
  	endfunction

   	//Creating a Dummy Variable for IPopt use
   	empty=[0];

   	//Calling the Ipopt function for solving the above problem
    [xopt,fopt,status,iter,cpu,obj_eval,dual,lambda1,zl,zu,gradient,hessian1] = solveminconp("f",A,b,Aeq,beq,lb,ub,no_nlc,no_nlic,"addnlc1","fGrad1","lHess1","addcGrad1",x0,options,empty)	
   
	//Calculating the values for the output   	
   	xopt = xopt';
    exitflag = status;
    output = struct("Iterations", [],"Cpu_Time",[],"Objective_Evaluation",[],"Dual_Infeasibility",[],"Message","");
   	output.Iterations = iter;
    output.Cpu_Time = cpu;
    output.Objective_Evaluation = obj_eval;
    output.Dual_Infeasibility = dual;
	lambda = struct("lower", zl,"upper",zu,"ineqlin",[],"eqlin",[],"ineqnonlin",[],"eqnonlin",[]);

	if (no_nlic ~= 0) then
		for i = 1:no_nlic
			lambda.ineqnonlin (i) = lambda1(i)
		end
		lambda.ineqnonlin = lambda.ineqnonlin'
	end
	
	if (no_nlec ~= 0) then
		j=1;
		for i = no_nlic+1 : no_nlc
			lambda.eqnonlin (j) = lambda1(i)
			j= j+1;
		end
		lambda.eqnonlin = lambda.eqnonlin'
	end
	
	if (Aeq ~=[]) then
		j=1;
		for i = no_nlc+1 : no_nlc + size(Aeq,1)
			lambda.eqlin (j) = lambda1(i)
			j= j+1;
		end
		lambda.eqlin = lambda.eqlin'
	end
	
	if (A ~=[]) then
		j=1;
		for i = no_nlc+ size(Aeq,1)+ 1 : no_nlc + size(Aeq,1) + size(A,1)
			lambda.ineqlin (j) = lambda1(i)
			j= j+1;
		end
		lambda.ineqlin = lambda.ineqlin'
	end
	
    //Converting hessian of order (1 x (numberOfVariables)^2) received from Ipopt to order (numberOfVariables x numberOfVariables)
    s1=size(gradient)
    for i =1:s1(2)
    	for j =1:s1(2)
			hessian(i,j)= hessian1(j+((i-1)*s1(2)))
		end
    end
    
    //In the cases of the problem not being solved, return NULL to the output matrices
    if( status~=0 & status~=1 & status~=2 & status~=3 & status~=4 & status~=7 ) then
		xopt=[];
		fopt=[];
		output = struct("Iterations", [],"Cpu_Time",[],"Message","");
		output.Iterations = iter;
    	output.Cpu_Time = cpu;
		lambda = struct("lower",[],"upper",[],"ineqlin",[],"eqlin",[],"ineqnonlin",[],"eqnonlin",[]);
		gradient=[];
		hessian=[];
    end
		

    //To print output message
    select status
    
    	case 0 then
        	printf("\nOptimal Solution Found.\n");
        	output.Message="Optimal Solution Found";
    	case 1 then
        	printf("\nMaximum Number of Iterations Exceeded. Output may not be optimal.\n");
        	output.Message="Maximum Number of Iterations Exceeded. Output may not be optimal";
    	case 2 then
       		printf("\nMaximum CPU Time exceeded. Output may not be optimal.\n");
       		output.Message="Maximum CPU Time exceeded. Output may not be optimal";
    	case 3 then
        	printf("\nStop at Tiny Step\n");
        	output.Message="Stop at Tiny Step";
    	case 4 then
        	printf("\nSolved To Acceptable Level\n");
        	output.Message="Solved To Acceptable Level";
    	case 5 then
        	printf("\nConverged to a point of local infeasibility.\n");
        	output.Message="Converged to a point of local infeasibility";
    	case 6 then
        	printf("\nStopping optimization at current point as requested by user.\n");
        	output.Message="Stopping optimization at current point as requested by user";
    	case 7 then
        	printf("\nFeasible point for square problem found.\n");
        	output.Message="Feasible point for square problem found";
    	case 8 then 
        	printf("\nIterates diverging; problem might be unbounded.\n");
        	output.Message="Iterates diverging; problem might be unbounded";
    	case 9 then
        	printf("\nRestoration Failed!\n");
        	output.Message="Restoration Failed!";
    	case 10 then
        	printf("\nError in step computation (regularization becomes too large?)!\n");
        	output.Message="Error in step computation (regularization becomes too large?)!";
    	case 11 then
        	printf("\nProblem has too few degrees of freedom.\n");
        	output.Message="Problem has too few degrees of freedom";
    	case 12 then
        	printf("\nInvalid option thrown back by Ipopt\n");
        	output.Message="Invalid option thrown back by Ipopt";
    	case 13 then
        	printf("\nNot enough memory.\n");
        	output.Message="Not enough memory";
    	case 15 then
        	printf("\nINTERNAL ERROR: Unknown SolverReturn value - Notify Ipopt Authors.\n");
        	output.Message="INTERNAL ERROR: Unknown SolverReturn value - Notify Ipopt Authors";
    	else
        	printf("\nInvalid status returned. Notify the Toolbox authors\n");
        	output.Message="Invalid status returned. Notify the Toolbox authors";
        	break;
        end
endfunction

