// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Authors: Animesh Baranawal
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

function [x,fval,maxfval,exitflag] = intfminimax(varargin)
    // Solves minimax constraint problem
    //
    // Calling Sequence
    //  xopt = intfminimax(fun,x0,intcon)
    //  xopt = intfminimax(fun,x0,intcon,A,b)
    //  xopt = intfminimax(fun,x0,intcon,A,b,Aeq,beq)
    //  xopt = intfminimax(fun,x0,intcon,A,b,Aeq,beq,lb,ub)
    //  xopt = intfminimax(fun,x0,intcon,A,b,Aeq,beq,lb,ub,nonlinfun)
    //  xopt = intfminimax(fun,x0,intcon,A,b,Aeq,beq,lb,ub,nonlinfun,options)
    //  [xopt, fval] = intfminimax(.....)
    //  [xopt, fval, maxfval]= intfminimax(.....)
    //  [xopt, fval, maxfval, exitflag]= intfminimax(.....)
    //
    // Parameters
    //  fun: The function to be minimized. fun is a function that accepts a vector x and returns a vector F, the objective functions evaluated at x.
	//   x0 : a vector of double, contains initial guess of variables.
	//   A : a matrix of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b. 
	//   intcon : a vector of integers, represents which variables are constrained to be integers
	//   b : a vector of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b.
	//   Aeq : a matrix of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   beq : a vector of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   lb : a vector of double, contains lower bounds of the variables.
	//   ub : a vector of double,  contains upper bounds of the variables.
    //  nonlinfun: function that computes the nonlinear inequality constraints c⋅x ≤ 0 and nonlinear equality constraints c⋅x = 0.
	//   xopt : a vector of double, the computed solution of the optimization problem.
	//   fopt : a double, the value of the function at x.
    //   maxfval: a 1x1 matrix of doubles, the maximum value in vector fval
	//   exitflag : The exit status. See below for details.
	//   output : The structure consist of statistics about the optimization. See below for details.
	//   lambda : The structure consist of the Lagrange multipliers at the solution of problem. See below for details.
	//
    // Description
    //  intfminimax minimizes the worst-case (largest) value of a set of multivariable functions, starting at an initial estimate. This is generally referred to as the minimax problem.
    //
    //  <latex>
    //  \min_{x} \max_{i} F_{i}(x)\: \textrm{such that} \:\begin{cases}
    //  & c(x) \leq 0 \\
    //  & ceq(x) = 0 \\
    //  & A.x \leq b \\
    //  & Aeq.x = beq \\
    //  & lb \leq x \leq ub
	//  & x_i \in \!\, \mathbb{Z}, i \in \!\, I
    //  \end{cases}
    //  </latex>
    //
    //  Currently, intfminimax calls intfmincon which uses the bonmin algorithm.
    //
    //  max-min problems can also be solved with intfminimax, using the identity
    //
    //  <latex>
    //  \max_{x} \min_{i} F_{i}(x) = -\min_{x} \max_{i} \left( -F_{i}(x) \right)
    //  </latex>
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
    //  The objective function must have header :
    //  <programlisting>
    //      F = fun(x)
    //  </programlisting>
    //  where x is a n x 1 matrix of doubles and F is a m x 1 matrix of doubles where m is the total number of objective functions inside F.
    //  On input, the variable x contains the current point and, on output, the variable F must contain the objective function values.
    //
    //  By default, the gradient options for intfminimax are turned off and and intfmincon does the gradient opproximation of objective function. In case the GradObj option is off and GradConstr option is on, intfminimax approximates Objective function gradient using numderivative toolbox.
    //
    //  If we can provide exact gradients, we should do so since it improves the convergence speed of the optimization algorithm.
    //
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
	// For more details on exitflag see the ipopt documentation, go to http://www.coin-or.org/bonmin/
	//
    // Examples
    //  // A basic case :
    //  // we provide only the objective function and the nonlinear constraint
    //  // function
    //  function f = myfun(x)
    //      f(1)= 2*x(1)^2 + x(2)^2 - 48*x(1) - 40*x(2) + 304;     //Objectives
    //      f(2)= -x(1)^2 - 3*x(2)^2;
    //      f(3)= x(1) + 3*x(2) -18;
    //      f(4)= -x(1) - x(2);
    //      f(5)= x(1) + x(2) - 8;
    //  endfunction
    //  // The initial guess
    //  x0 = [0.1,0.1];
    //  // The expected solution : only 4 digits are guaranteed
    //  xopt = [4 4]
    //  fopt = [0 -64 -2 -8 0]
	//	intcon = [1]
    //  maxfopt = 0
    //  // Run fminimax
    //  [x,fval,maxfval,exitflag] = intfminimax(myfun, x0,intcon)
	// // Press ENTER to continue
	//
    // Examples
    //  // A case where we provide the gradient of the objective
    //  // functions and the Jacobian matrix of the constraints.
    //  // The objective function and its gradient
	//      function [f,G] = myfun(x)
    //          f(1)= 2*x(1)^2 + x(2)^2 - 48*x(1) - 40*x(2) + 304;
    //          f(2)= -x(1)^2 - 3*x(2)^2;
    //          f(3)= x(1) + 3*x(2) -18;
    //          f(4)= -x(1) - x(2);
    //          f(5)= x(1) + x(2) - 8;
    //          G = [ 4*x(1) - 48, -2*x(1), 1, -1, 1;
    //                2*x(2) - 40, -6*x(2), 3, -1, 1; ]'
    //      endfunction
    //      // The nonlinear constraints 
    //      function [c,ceq,DC,DCeq] = confun(x)
    //          // Inequality constraints
    //          c = [1.5 + x(1)*x(2) - x(1) - x(2), -x(1)*x(2) - 10] 
    //          // No nonlinear equality constraints
    //          ceq=[]
    //          DC= [x(2)-1, -x(2);
    //          	x(1)-1, -x(1)]'
    //          DCeq = []'
    //      endfunction
    //      // Test with both gradient of objective and gradient of constraints
    //      minimaxOptions = list("GradObj","on","GradCon","on");
    //      // The initial guess
    //      x0 = [0,10];
    //      intcon = [2]
    //      // Run intfminimax
    //      [x,fval,maxfval,exitflag] = intfminimax(myfun,x0,intcon,[],[],[],[],[],[], confun, minimaxOptions)
    // Authors
    // Harpreet Singh

    // Check number of input and output arguments
    [minmaxLhs,minmaxRhs] = argn()
    Checkrhs("fminimax", minmaxRhs, [2 3 5 7 9 10 11])
    Checklhs("fminimax", minmaxLhs, 1:7)

    // Proper initialisation of objective function
    minmaxObjfun = varargin(1)
    Checktype("fminimax", minmaxObjfun, "minmaxObjfun", 1, "function")

    // Proper initialisation of starting point
    minmaxStartpoint = varargin(2)
    Checktype("fminimax", minmaxStartpoint, "minmaxStartpoint", 2, "constant")


    minmaxNumvar = length(minmaxStartpoint)
    Checkvector("fminimax", minmaxStartpoint, "minmaxStartpoint", 2, minmaxNumvar)
    minmaxStartpoint = minmaxStartpoint(:)

    if(minmaxRhs < 3) then // if A and b are not provided, declare as empty
        intcon = 0;
    else
		intcon = varargin(3);
    end

    // Proper initialisation of A and b
    if(minmaxRhs < 4) then // if A and b are not provided, declare as empty
        minmaxA = []
        minmaxB = []
    else
        minmaxA = varargin(4)
        minmaxB = varargin(5)
    end

    Checktype("fminimax", minmaxA, "A", 4, "constant")
    Checktype("fminimax", minmaxB, "b", 5, "constant")

    // Check if A and b of proper dimensions
    if(minmaxA <> [] & minmaxB == []) then
        errmsg = msprintf(gettext("%s: Incompatible input arguments #%d and #%d: matrix A is empty, but the column vector b is not empty"), "fminimax", 4, 5)
        error(errmsg)
    end

    if(minmaxA == [] & minmaxB <> []) then
        errmsg = msprintf(gettext("%s: Incompatible input arguments #%d and #%d: matrix A is not empty, but the column vector b is empty"), "fminimax", 4, 5)
        error(errmsg)
    end

    minmaxNumrowA = size(minmaxA,"r")
    if(minmaxA <> []) then
        Checkdims("fminimax", minmaxA, "A", 4, [minmaxNumrowA minmaxNumvar])
        Checkvector("fminimax", minmaxB, "b", 5, minmaxNumrowA)
        minmaxB = minmaxB(:)
    end

    // Proper initialisation of Aeq and beq
    if(minmaxRhs < 6) then // if Aeq and beq are not provided, declare as empty
        minmaxAeq = []
        minmaxBeq = []
    else
        minmaxAeq = varargin(6)
        minmaxBeq = varargin(7)
    end

    Checktype("fminimax", minmaxAeq, "Aeq", 6, "constant")
    Checktype("fminimax", minmaxBeq, "beq", 7, "constant")

    // Check if Aeq and beq of proper dimensions
    if(minmaxAeq <> [] & minmaxBeq == []) then
        errmsg = msprintf(gettext("%s: Incompatible input arguments #%d and #%d: matrix Aeq is empty, but the column vector beq is not empty"), "fminimax", 6, 7)
        error(errmsg)
    end

    if(minmaxAeq == [] & minmaxBeq <> []) then
        errmsg = msprintf(gettext("%s: Incompatible input arguments #%d and #%d: matrix Aeq is not empty, but the column vector beq is empty"), "fminimax", 6, 7)
        error(errmsg)
    end

    minmaxNumrowAeq = size(minmaxAeq,"r")
    if(minmaxAeq <> []) then
        Checkdims("fminimax", minmaxAeq, "Aeq", 6, [minmaxNumrowAeq minmaxNumvar])
        Checkvector("fminimax", minmaxBeq, "beq", 7, minmaxNumrowAeq)
        minmaxBeq = minmaxBeq(:)
    end

    // Proper initialisation of minmaxLb and minmaxUb
    if(minmaxRhs < 6) then // if minmaxLb and minmaxUb are not provided, declare as empty
        minmaxLb = []
        minmaxUb = []
    else
        minmaxLb = varargin(6)
        minmaxUb = varargin(7)
    end

    Checktype("fminimax", minmaxLb, "lb", 6, "constant")
    Checktype("fminimax", minmaxUb, "ub", 7, "constant")

    // Check dimensions of minmaxLb and minmaxUb
    if(minmaxLb <> []) then
        Checkvector("fminimax", minmaxLb, "lb", 8, minmaxNumvar)
        minmaxLb = minmaxLb(:)
    end

    if(minmaxUb <> []) then
        Checkvector("fminimax", minmaxUb, "ub", 9, minmaxNumvar)
        minmaxUb = minmaxUb(:)
    end

    // Proper Initialisation of minmaxNonlinfun
    if(minmaxRhs < 10) then // if minmaxNonlinfun is not provided, declare as empty
        minmaxNonlinfun = []
    else
        minmaxNonlinfun = varargin(10)
    end
	if(minmaxNonlinfun<>[]) then
	    Checktype("fminimax", minmaxNonlinfun, "nonlinfun", 10, "function")
	end

    //To check, Whether minimaxOptions is been entered by user
    if ( minmaxRhs<11 ) then
        minmaxUserOptions = list();
    else
        minmaxUserOptions = varargin(11); //Storing the 3rd Input minmaxUserOptionseter in intermediate list named 'minmaxUserOptions'
    end

    //If minimaxOptions is entered then checking its type for 'list'
    if (type(minmaxUserOptions) ~= 15) then
        errmsg = msprintf(gettext("%s: minimaxOptions (10th parameter) should be a list"), "fminimax");
        error(errmsg);
    end

    //If minimaxOptions is entered then checking whether even number of entires are entered
    if (modulo(size(minmaxUserOptions),2)) then
        errmsg = msprintf(gettext("%s: Size of minimaxOptions (list) should be even"), "fminimax");
        error(errmsg);
    end

  /////////////// To check integer //////////////////////
  for i=1:size(intcon,1)
      if(intcon(i)>minmaxNumvar) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be less than the number of variables"), "intfminimax");
        error(errmsg);
      end

      if (intcon(i)<0) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be greater than 0 "), "intfminimax");
        error(errmsg);
      end

      if(modulo(intcon(i),1)) then
        errmsg = msprintf(gettext("%s: The values inside intcon should be an integer "), "intfminimax");
        error(errmsg);
      end
  end
  
    //If minimaxOptions is entered then checking its type for 'list'
    if (type(minmaxUserOptions) ~= 15) then
        errmsg = msprintf(gettext("%s: minimaxOptions (10th parameter) should be a list"), "intfminimax");
        error(errmsg);
    end

    //If minimaxOptions is entered then checking whether even number of entires are entered
    if (modulo(size(minmaxUserOptions),2)) then
        errmsg = msprintf(gettext("%s: Size of minimaxOptions (list) should be even"), "intfminimax");
        error(errmsg);
    end

minmaxoptions = list('integertolerance',1d-06,'maxnodes',2147483647,'cputime',1d10,'allowablegap',0,'maxiter',2147483647,'gradobj',"off",'gradcon',"off")

  //Pushing minmaxUserOptions into default value
  
  for i = 1:(size(minmaxUserOptions))/2
    select convstr(minmaxUserOptions(2*i-1),'l')
      case 'integertolerance' then
        Checktype("intfminimax_options", minmaxUserOptions(2*i), minmaxUserOptions(2*i-1), 2*i, "constant");
        minmaxoptions(2) = minmaxUserOptions(2*i);
      case 'maxnodes' then
        Checktype("intfminimax_options", minmaxUserOptions(2*i), minmaxUserOptions(2*i-1), 2*i, "constant");
        minmaxoptions(4) = minmaxUserOptions(2*i);
      case 'cputime' then 
        Checktype("intfminimax_options", minmaxUserOptions(2*i), minmaxUserOptions(2*i-1), 2*i, "constant");
        minmaxoptions(6) = minmaxUserOptions(2*i);
      case 'allowablegap' then
        Checktype("intfminimax_options", minmaxUserOptions(2*i), minmaxUserOptions(2*i-1), 2*i, "constant");
        minmaxoptions(8) = minmaxUserOptions(2*i);
      case 'maxiter' then
        Checktype("intfminimax_options", minmaxUserOptions(2*i), minmaxUserOptions(2*i-1), 2*i, "constant");
        minmaxoptions(10) = minmaxUserOptions(2*i);
      case 'gradobj' then
        Checktype("intfminimax_options", minmaxUserOptions(2*i), minmaxUserOptions(2*i-1), 2*i, "string");
        if(convstr(minmaxUserOptions(2*i),'l') == "on") then
          minmaxoptions(12) = "on"
        elseif(convstr(minmaxUserOptions(2*i),'l') == "off") then
          minmaxoptions(12) = "off"
        else
          error(999, 'Unknown string passed in gradobj.');
        end
      case 'gradcon' then
        Checktype("intfminimax_options", minmaxUserOptions(2*i), minmaxUserOptions(2*i-1), 2*i, "string");
        if(convstr(minmaxUserOptions(2*i),'l') == "on") then
          minmaxoptions(14) = "on"
        elseif(convstr(minmaxUserOptions(2*i),'l') == "off") then
          minmaxoptions(14) = "off"
        else
          error(999, 'Unknown string passed in gradcon.');
        end
      else   
        error(999, 'Unknown string argument passed.');
    end
  end 

    // Reformulating the problem fminimax to fmincon
    minmaxObjfunval = minmaxObjfun(minmaxStartpoint)
    minmaxStartpoint(minmaxNumvar+1) = max(minmaxObjfunval)

    if(minmaxA <> []) then
        minmaxA = [minmaxA, zeros(minmaxNumrowA,1)]
    end
    if(minmaxAeq <> []) then
        minmaxAeq = [minmaxAeq, zeros(minmaxNumrowAeq,1)]
    end
    if(minmaxLb <> []) then
        minmaxLb(minmaxNumvar+1) = -%inf
    end
    if(minmaxUb <> []) then
        minmaxUb(minmaxNumvar+1) = +%inf
    end

    // function handle defining the additional inequalities
    function temp = minmaxAddIneq(z)
        temp = minmaxObjfun(z) - z(minmaxNumvar+1)
    endfunction

    // function handle defining minmaxNonlinfun derivative using numderivative
    function [dc,dceq] = minmaxNonlinDer(z)
        // function handle extracting c and ceq components from minmaxNonlinfun
        function foo = minmaxC(z)
            [foo,tmp1] = minmaxNonlinfun(z)
            foo = foo'
        endfunction

        function foo = minmaxCEQ(z)
            [tmp1,foo] = minmaxNonlinfun(z)
            foo = foo'
        endfunction

        dc = numderivative(minmaxC,z)
        dceq = numderivative(minmaxCEQ,z)
    endfunction

    // function handle defining new objective function
    function newfunc = newObjfun(z)
        newfunc = z(minmaxNumvar+1)
    endfunction

    // function handle defining new minmaxNonlinfun function
    function [nc,nceq,dnc,dnceq] = newNonlinfun(z)
		dnc = [];
		dnceq = [];
		nc = [];
		nceq= [];
		if (minmaxNonlinfun<>[]) then
	        [nc,nceq] = minmaxNonlinfun(z)
		end
        // add inequalities of the form Fi(x) - y <= 0
        tmp = [minmaxObjfun(z) - z(minmaxNumvar+1)]'
        nc = [nc, tmp]
		if(options(14) =="on") then
			[temp1,temp2,dnc, dnceq] = minmaxNonlinfun(z)
            dnc = [dnc, zeros(size(dnc,'r'),1)]
            dnceq = [dnceq, zeros(size(dnceq,'r'),1)]
        else
            // else use numderivative method to calculate gradient of constraints
			if (minmaxNonlinfun<>[]) then
            	[dnc, dnceq] = minmaxNonlinDer(z)
			end
		end
		
		if(options(12) =="on") then
			[temp,derObjfun] = minmaxObjfun(z);
            mderObjfun = [derObjfun, -1*ones(size(derObjfun,'r'),1)];
            dnc = [dnc; mderObjfun];
        else
            // else use numderivative to calculate gradient of set of obj functions
            derObjfun = numderivative(minmaxAddIneq,z)
            dnc = [dnc; derObjfun]
		end
    endfunction

	if( minmaxoptions(12)=="on"| minmaxoptions(12)=="on") then
		options(14)="on";
	end

	minmaxoptions(12)="off";
    [x,fval,exitflag,gradient,hessian] = ...
    intfmincon(newObjfun,minmaxStartpoint,intcon,minmaxA,minmaxB,minmaxAeq,minmaxBeq,minmaxLb,minmaxUb,newNonlinfun,minmaxoptions)

    x = x(1:minmaxNumvar)
    fval = minmaxObjfun(x)
    maxfval = max(fval)
endfunction
