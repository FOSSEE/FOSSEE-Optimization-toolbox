// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Harpreet Singh, Pranav Deshpande and Akshay Miterani
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

function [xopt,fopt,exitflag,gradient,hessian] = intfmincon (varargin)
	// Solves a constrainted multi-variable mixed integer non linear programming problem
	//
	//   Calling Sequence
	//   xopt = intfmincon(f,x0,intcon,A,b)
	//   xopt = intfmincon(f,x0,intcon,A,b,Aeq,beq)
	//   xopt = intfmincon(f,x0,intcon,A,b,Aeq,beq,lb,ub)
	//   xopt = intfmincon(f,x0,intcon,A,b,Aeq,beq,lb,ub,nlc)
	//   xopt = intfmincon(f,x0,intcon,A,b,Aeq,beq,lb,ub,nlc,options)
	//   [xopt,fopt] = intfmincon(.....)
	//   [xopt,fopt,exitflag]= intfmincon(.....)
	//   [xopt,fopt,exitflag,gradient]=intfmincon(.....)
	//   [xopt,fopt,exitflag,gradient,hessian]=intfmincon(.....)
	//
	//   Parameters
	//   f : a function, representing the objective function of the problem 
	//   x0 : a vector of doubles, containing the starting values of variables. 
	//   intcon : a vector of integers, represents which variables are constrained to be integers
	//   A : a matrix of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b. 
	//   b : a vector of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b.
	//   Aeq : a matrix of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   beq : a vector of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   lb : Lower bounds, specified as a vector or array of double. lb represents the lower bounds elementwise in lb ≤ x ≤ ub.
	//   ub : Upper bounds, specified as a vector or array of double. ub represents the upper bounds elementwise in lb ≤ x ≤ ub.
	//   nlc : a function, representing the Non-linear Constraints functions(both Equality and Inequality) of the problem. It is declared in such a way that non-linear inequality constraints are defined first as a single row vector (c), followed by non-linear equality constraints as another single row vector (ceq). Refer Example for definition of Constraint function.
	//   options : a list, containing the option for user to specify. See below for details. 
	//   xopt : a vector of doubles, containing the the computed solution of the optimization problem.
	//   fopt : a scalar of double, containing the the function value at x.
	//   exitflag : a scalar of integer, containing the flag which denotes the reason for termination of algorithm. See below for details.
	//   gradient : a vector of doubles, containing the Objective's gradient of the solution.
	//   hessian  : a matrix of doubles, containing the Objective's hessian of the solution.
	//
	//   Description
  	//   Search the minimum of a mixed integer constrained optimization problem specified by :
  	//   Find the minimum of f(x) such that 
  	//
  	//   <latex>
  	//    \begin{eqnarray}
  	//    &\mbox{min}_{x}
  	//    & f(x) \\
  	//    & \text{subject to} & A*x \leq b \\
  	//    & & Aeq*x \ = beq\\
  	//	  & & c(x) \leq  0\\
  	//    & & ceq(x) \ =  0\\
  	//    & & lb \leq x \leq ub \\
	//    & & x_i \in \!\, \mathbb{Z}, i \in \!\, I
  	//    \end{eqnarray}
  	//   </latex>
	//
	//   The routine calls Bonmin for solving the Bounded Optimization problem, Bonmin is a library written in C++.
	//
	// The options allows the user to set various parameters of the Optimization problem. 
	// It should be defined as type "list" and contains the following fields.
	// <itemizedlist>
	//   <listitem>Syntax : options= list("IntegerTolerance", [---], "MaxNodes",[---], "MaxIter", [---], "AllowableGap",[---] "CpuTime", [---],"gradobj", "off", "hessian", "off" );</listitem>
  //   <listitem>IntegerTolerance : a Scalar, a number with that value of an integer is considered integer..</listitem>
  //   <listitem>MaxNodes : a Scalar, containing the Maximum Number of Nodes that the solver should search.</listitem>
	//   <listitem>CpuTime : a Scalar, containing the Maximum amount of CPU Time that the solver should take.</listitem>
	//   <listitem>AllowableGap : a Scalar, to stop the tree search when the gap between the objective value of the best known solution is reached.</listitem>
  //   <listitem>MaxIter : a Scalar, containing the Maximum Number of Iteration that the solver should take.</listitem>
  //   <listitem>gradobj : a string, to turn on or off the user supplied objective gradient.</listitem>
  //   <listitem>hessian : a Scalar, to turn on or off the user supplied objective hessian.</listitem>
	//   <listitem>Default Values : options = list('integertolerance',1d-06,'maxnodes',2147483647,'cputime',1d10,'allowablegap',0,'maxiter',2147483647,'gradobj',"off",'hessian',"off")</listitem>
	// </itemizedlist>
	//
	// The exitflag allows to know the status of the optimization which is given back by Ipopt.
	// <itemizedlist>
	//   <listitem>exitflag=0 : Optimal Solution Found </listitem>
	//   <listitem>exitflag=1 : InFeasible Solution.</listitem>
	//   <listitem>exitflag=2 : Objective Function is Continuous Unbounded.</listitem>
	//   <listitem>exitflag=3 : Limit Exceeded.</listitem>
	//   <listitem>exitflag=4 : User Interrupt.</listitem>
	//   <listitem>exitflag=5 : MINLP Error.</listitem>
	// </itemizedlist>
	//
	// For more details on exitflag see the Bonmin documentation, go to http://www.coin-or.org/Bonmin
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
  	//  function [y,dy]=f(x)
  	//		y=-x(1)-x(2)/3;
	//		dy= [-1,-1/3];
  	//	endfunction
  	//	//Starting point, linear constraints and variable bounds  
  	//	x0=[0 , 0];
	//  intcon = [1]
  	//	A=[1,1 ; 1,1/4 ; 1,-1 ; -1/4,-1 ; -1,-1 ; -1,1];
  	//	b=[2;1;2;1;-1;2];
  	//	Aeq=[1,1];
  	// 	beq=[2];
  	//	lb=[];
  	//	ub=[];
  	//  nlc=[];
  	//	//Options
  	//	options=list("GradObj", "on");
  	//    //Calling Ipopt
  	//	[x,fval,exitflag,grad,hessian] =intfmincon(f, x0,intcon,A,b,Aeq,beq,lb,ub,nlc,options)
	// // Press ENTER to continue
  	//
  	// Examples
  	//	//Find x in R^3 such that it minimizes:
  	//    //f(x)= x1*x2 + x2*x3
  	//    //x0=[0.1 , 0.1 , 0.1]
  	//    //constraint-1 (c1): x1^2 - x2^2 + x3^2 <= 2
  	//    //constraint-2 (c2): x1^2 + x2^2 + x3^2 <= 10  
  	//    //Objective function to be minimised
  	//  function [y,dy]=f(x)
  	//		y=x(1)*x(2)+x(2)*x(3);
  	//		dy= [x(2),x(1)+x(3),x(2)];
  	//	endfunction
  	//	//Starting point, linear constraints and variable bounds  
  	//	x0=[0.1 , 0.1 , 0.1];
	//	intcon = [2]
  	//	A=[];
  	//	b=[];
  	//	Aeq=[];
  	// 	beq=[];
  	//	lb=[];
  	//	ub=[];
  	//	//Nonlinear constraints  
  	//	function [c,ceq,cg,cgeq]=nlc(x)
  	//		c = [x(1)^2 - x(2)^2 + x(3)^2 - 2 , x(1)^2 + x(2)^2 + x(3)^2 - 10];
  	//  	ceq = [];
  	//		cg=[2*x(1) , -2*x(2) , 2*x(3) ; 2*x(1) , 2*x(2) , 2*x(3)];
  	//		cgeq=[];
  	//	endfunction
  	//	//Options  
  	//	options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "on","GradCon", "on");
  	//    //Calling Ipopt
  	//	[x,fval,exitflag,output] =intfmincon(f, x0,intcon,A,b,Aeq,beq,lb,ub,nlc,options)
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
	//	intcon = [3]
  	//	A=[];
  	//	b=[];
  	//	Aeq=[];
  	// 	beq=[];
  	//	lb=[];
  	//	ub=[0,0,0];
  	//	//Options
  	//	options=list("MaxIter", [1500], "CpuTime", [500]);
  	//    //Calling Ipopt
  	//	[x,fval,exitflag,grad,hessian] =intfmincon(f, x0,intcon,A,b,Aeq,beq,lb,ub,[],options)
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
  	//    function [y,dy]=f(x)
  	//		y=x(1)*x(2)+x(2)*x(3);
  	//		dy= [x(2),x(1)+x(3),x(2)];
  	//	  endfunction
  	//	//Starting point, linear constraints and variable bounds  
  	//	x0=[1,1,1];
	//	intcon = [2]
  	//	A=[];
  	//	b=[];
  	//	Aeq=[];
  	// 	beq=[];
  	//	lb=[0 0.2,-%inf];
  	//	ub=[0.6 %inf,1];
  	//	//Nonlinear constraints  
  	//	function [c,ceq,cg,cgeq]=nlc(x)
  	//		c=[x(1)^2-1,x(1)^2+x(2)^2-1,x(3)^2-1];
  	//		ceq=[x(1)^3-0.5,x(2)^2+x(3)^2-0.75];
  	//		cg = [2*x(1),0,0;2*x(1),2*x(2),0;0,0,2*x(3)];
  	//		cgeq = [3*x(1)^2,0,0;0,2*x(2),2*x(3)];
  	//	endfunction
   	//	//Options  
  	//	options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", "on","GradCon", "on");
  	//    //Calling Ipopt
  	//	[x,fval,exitflag,grad,hessian] =intfmincon(f, x0,intcon,A,b,Aeq,beq,lb,ub,nlc,options)
 	// // Press ENTER to continue
	// Authors
	// Harpreet Singh

	//To check the number of input and output arguments
 	[lhs , rhs] = argn();
	
	//To check the number of arguments given by the user
 	if ( rhs<4 | rhs>11 ) then
  		errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be int [4 5] "), "intfmincon", rhs);
  		error(errmsg);
 	end
 
	//Storing the Input Parameters  
 	fun = varargin(1);
 	x0     = varargin(2);
	intcon = varargin(3); 
	A      = varargin(4);
	b      = varargin(5);
	Aeq    = [];
	beq    = [];
	lb       = [];
	ub       = [];
	nlc      = [];
    
  if (rhs>5) then
    Aeq    = varargin(6);
    beq    = varargin(7);
  end

  if (rhs>7) then
    lb       = varargin(8);
    ub       = varargin(9);
  end

  if (rhs>9) then
    nlc      = varargin(10);
  end

  param = list();
  //To check whether options has been entered by user   
  if ( rhs>  10) then
      param =varargin(11);
  end

  //To check whether the Input arguments
  Checktype("intfmincon", fun, "fun", 1, "function");
  Checktype("intfmincon", x0, "x0", 2, "constant");
  Checktype("intfmincon", intcon, "intcon", 3, "constant");
  Checktype("intfmincon", A, "A", 4, "constant");
  Checktype("intfmincon", b, "b", 5, "constant");
  Checktype("intfmincon", Aeq, "Aeq", 6, "constant");
  Checktype("intfmincon", beq, "beq", 7, "constant");
  Checktype("intfmincon", lb, "lb", 8, "constant");
  Checktype("intfmincon", ub, "ub", 9, "constant");
  Checktype("intfmincon", nlc, "nlc", 10, ["constant","function"]);
  Checktype("intfmincon", param, "options", 11, "list");  


  nbVar = size(x0,"*");
  if(nbVar==0) then
    errmsg = msprintf(gettext("%s: x0 cannot be an empty"), "intfmincon");
    error(errmsg);    
  end

  if(size(lb,"*")==0) then
    lb = repmat(-%inf,nbVar,1);
  end

  if(size(ub,"*")==0) then
    ub = repmat(%inf,nbVar,1);
  end

  //////////////// To Check linear constraints /////////

  //To check for correct size of A(3rd paramter)
  if(size(A,2)~=nbVar & size(A,2)~=0) then
    errmsg = msprintf(gettext("%s: Expected Matrix of size (No of linear inequality constraints X No of Variables) or an Empty Matrix for Linear Inequality Constraint coefficient Matrix A"), "intfmincon");
    error(errmsg);
  end
  nbConInEq=size(A,"r");

  //To check for the correct size of Aeq (5th paramter)
  if(size(Aeq,2)~=nbVar & size(Aeq,2)~=0) then
    errmsg = msprintf(gettext("%s: Expected Matrix of size (No of linear equality constraints X No of Variables) or an Empty Matrix for Linear Equality Constraint coefficient Matrix Aeq"), "intfmincon");
    error(errmsg);
  end
  nbConEq=size(Aeq,"r");

  ///////////////// To check vectors ///////////////// 

  Checkvector("intfmincon", x0, "x0", 2, nbVar);
  x0 = x0(:);
  if(size(intcon,"*")) then
    Checkvector("intfmincon", intcon, "intcon", 3, size(intcon,"*"))
    intcon = intcon(:);
  end
  if(nbConInEq) then
    Checkvector("intfmincon", b, "b", 5, nbConInEq);
    b = b(:);
  end
  if(nbConEq) then
    Checkvector("intfmincon", beq, "beq", 7, nbConEq);
    beq = beq(:);
  end
  Checkvector("intfmincon", lb, "lb", 8, nbVar);
  lb = lb(:);
  
  Checkvector("intfmincon", ub, "ub", 9, nbVar);
  ub = ub(:);

  /////////////// To check integer //////////////////////
  for i=1:size(intcon,1)
      if(intcon(i)>nbVar) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be less than the number of variables"), "intfmincon");
        error(errmsg);
      end

      if (intcon(i)<0) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be greater than 0 "), "intfmincon");
        error(errmsg);
      end

      if(modulo(intcon(i),1)) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be an integer "), "intfmincon");
        error(errmsg);
      end
  end
  
options = list('integertolerance',1d-06,'maxnodes',2147483647,'cputime',10000000,'allowablegap',0,'maxiter',2147483647,'gradobj',"off",'hessian',"off",'gradcon',"off")

  //Pushing param into default value
  
  for i = 1:(size(param))/2
    select convstr(param(2*i-1),'l')
      case 'integertolerance' then
        Checktype("intfmincon_options", param(2*i), param(2*i-1), 2*i, "constant");
        options(2) = param(2*i);
      case 'maxnodes' then
        Checktype("intfmincon_options", param(2*i), param(2*i-1), 2*i, "constant");
        options(4) = options(2*i);
      case 'cputime' then 
        Checktype("intfmincon_options", param(2*i), param(2*i-1), 2*i, "constant");
        options(6) = options(2*i);
      case 'allowablegap' then
        Checktype("intfmincon_options", param(2*i), param(2*i-1), 2*i, "constant");
        options(8) = options(2*i);
      case 'maxiter' then
        Checktype("intfmincon_options", param(2*i), param(2*i-1), 2*i, "constant");
        options(10) = options(2*i);
      case 'gradobj' then
        Checktype("intfmincon_options", param(2*i), param(2*i-1), 2*i, "string");
        if(convstr(param(2*i),'l') == "on") then
          options(12) = "on"
        elseif(convstr(param(2*i),'l') == "off") then
          options(12) = "off"
        else
          error(999, 'Unknown string passed in gradobj.');
        end
      case 'hessian' then
        Checktype("intfmincon_options", param(2*i), param(2*i-1), 2*i, "function");
        options(14) = param(2*i);
      case 'gradcon' then
        Checktype("intfmincon_options", param(2*i), param(2*i-1), 2*i, "string");
        if(convstr(param(2*i),'l') == "on") then
          options(16) = "on"
        elseif(convstr(param(2*i),'l') == "off") then
          options(16) = "off"
        else
          error(999, 'Unknown string passed in gradcon.');
        end
      else    
        error(999, 'Unknown string argument passed.');
    end
  end 

  ///////////////// Functions Check /////////////////

  //To check the match between f (1st Parameter) and x0 (2nd Parameter)
  if(execstr('init=fun(x0)','errcatch')==21) then
    errmsg = msprintf(gettext("%s: Objective function and x0 did not match"), "intfmincon");
    error(errmsg);
  end

  if(options(12) == "on") then
       if(execstr('[grad_y,grad_dy]=fun(x0)','errcatch')==59) then
        errmsg = msprintf(gettext("%s: Gradient of objective function is not provided"), "intfmincon");
        error(errmsg);
      end
	if(grad_dy<>[]) then
      Checkvector("intfmincon_options", grad_dy, "dy", 12, nbVar);
	end
  end

  if(options(14) == "on") then
      if(execstr('[hessian_y,hessian_dy,hessian]=fun(x0)','errcatch')==59) then
        errmsg = msprintf(gettext("%s: Gradient of objective function is not provided"), "intfmincon");
        error(errmsg);
      end
      if ( ~isequal(size(hessian) == [nbVar nbVar]) ) then
        errmsg = msprintf(gettext("%s: Size of hessian should be nbVar X nbVar"), "intfmincon");
        error(errmsg);
      end
  end

  numNlic = 0;
  numNlec = 0;
  numNlc = 0;

  if (type(nlc) == 13 | type(nlc) == 11) then
	[sample_c,sample_ceq] = nlc(x0);
      if(execstr('[sample_c,sample_ceq] = nlc(x0)','errcatch')==21) then
        errmsg = msprintf(gettext("%s: Non-Linear Constraint function and x0 did not match"), "intfmincon");
        error(errmsg);
      end
      numNlic = size(sample_c,"*");
      numNlec = size(sample_ceq,"*");
      numNlc = numNlic + numNlec;
  end

  /////////////// Creating conLb and conUb ////////////////////////

  conLb = [repmat(-%inf,numNlic,1);repmat(0,numNlec,1);repmat(-%inf,nbConInEq,1);beq;]
  conUb = [repmat(0,numNlic,1);repmat(0,numNlec,1);b;beq;]

  //Converting the User defined Objective function into Required form (Error Detectable)
  function [y,check] = _f(x)
    try
      y=fun(x)
      [y,check] =  checkIsreal(y)
    catch
      y=0;
      check=1;
    end
   endfunction
    
   //Defining an inbuilt Objective gradient function 
    function [dy,check] = _gradf(x) 
      if (options(12) =="on") then
        try
          [y,dy]=fun(x)
          [dy,check] =  checkIsreal(dy)
        catch
          dy = 0;
          check=1;
        end
      else
          try 
            dy=numderivative(fun,x)
            [dy,check] =  checkIsreal(dy)
          catch
            dy=0;
          check=1;
          end
      end
    endfunction

    function [y,check] = _addnlc(x)
    x= x(:)
	c = []
	ceq = []
      try
        if((type(nlc) == 13 | type(nlc) == 11) & numNlc~=0) then
           [c,ceq]=nlc(x)
        end
        ylin = [A*x;Aeq*x];
        y = [c(:);ceq(:);ylin(:);];
        [y,check] = checkIsreal(y)
      catch
        y=0;
        check=1;
      end
    endfunction

    //Defining an inbuilt jacobian of constraints function 
    function [dy,check] = _gradnlc(x) 
      if (options(16) =="on") then
        try
          [y1,y2,dy1,dy2]=nlc(x)
		  //Adding derivative of Linear Constraint
		  dylin = [A;Aeq]
		  dy = [dy1;dy2;dylin];
          [dy,check] = checkIsreal(dy)
        catch
          dy = 0;
          check=1;
        end
      else
          try
            dy=numderivative(_addnlc,x)
            [dy,check] = checkIsreal(dy)
          catch
            dy=0;
            check=1;
          end
      end
    endfunction

  //Defining a function to calculate Hessian if the respective user entry is OFF 
  function [hessy,check]=_gradhess(x,obj_factor,lambda)
  x=x(:);
    if (type(options(14)) == "function") then
      try
        [obj,dy,hessy] = fun(x,obj_factor,lambda)
        [hessy,check] =  checkIsreal(hessy)
      catch
        hessy = 0;
        check=1;
      end
    else
      try
        [dy,hessfy]=numderivative(_f,x)
        hessfy = matrix(hessfy,nbVar,nbVar)
        if((type(nlc) == 13 | type(nlc) == 11) & numNlc~=0) then
			function y = nlc2(x)
				[c, ceq] = nlc(x);
				c = matrix(c,-1,1);
				y = [c; ceq];
			endfunction			
        	[dy,hessny]=numderivative(nlc2,x)
        end
        hessianc = lambda(1)*matrix(hessny(1,:),nbVar,nbVar);
        for i = 2:numNlc
            hessianc = hessianc + lambda(i)*matrix(hessny(i,:),nbVar,nbVar);
        end
        hessy = obj_factor*hessfy + hessianc;
        [hessy,check] =  checkIsreal(hessy)
      catch
        hessy=0;
        check=1;
      end
    end
  endfunction

    intconsize = size(intcon,"*")

	[xopt,fopt,exitflag] = inter_fmincon("_f","_gradf","_addnlc","_gradnlc","_gradhess",x0,lb,ub,conLb,conUb,intcon,options,nbConInEq+nbConEq);

    //In the cases of the problem not being solved, return NULL to the output matrices
    if( exitflag~=0 & exitflag~=3 ) then
      gradient = [];
      hessian = [];
    else
      [ gradient, hessian] = numderivative(_f, xopt)
    end

    //To print output message
    select exitflag
    
    case 0 then
        printf("\nOptimal Solution Found.\n");
    case 1 then
        printf("\nInFeasible Solution.\n");
    case 2 then
        printf("\nObjective Function is Continuous Unbounded.\n");
    case 3 then
        printf("\nStop at Tiny Step\n");
    case 4 then
        printf("\nUser Interrupt.\n");
    case 5 then
        printf("\nMINLP Error.\n");
    else
        printf("\nInvalid status returned. Notify the Toolbox authors\n");
        break;
    end
endfunction

function [y, check] = checkIsreal(x)
  if ((~isreal(x))) then
    y = 0
    check=1;
  else
    y = x;
    check=0;
  end 
endfunction
