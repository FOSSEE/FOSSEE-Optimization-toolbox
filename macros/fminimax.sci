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

function [x,fval,maxfval,exitflag,output,lambda] = fminimax(varargin)
    // Solves minimax constraint problem
    //
    // Calling Sequence
    //  xopt = fminimax(fun,x0)
    //  xopt = fminimax(fun,x0,A,b)
    //  xopt = fminimax(fun,x0,A,b,Aeq,beq)
    //  xopt = fminimax(fun,x0,A,b,Aeq,beq,lb,ub)
    //  xopt = fminimax(fun,x0,A,b,Aeq,beq,lb,ub,nonlinfun)
    //  xopt = fminimax(fun,x0,A,b,Aeq,beq,lb,ub,nonlinfun,options)
    //  [xopt, fval] = fminimax(.....)
    //  [xopt, fval, maxfval]= fminimax(.....)
    //  [xopt, fval, maxfval, exitflag]= fminimax(.....)
    //  [xopt, fval, maxfval, exitflag, output]= fminimax(.....)
    //  [xopt, fval, maxfval, exitflag, output, lambda]= fminimax(.....)
    //
    // Parameters
    //  fun: The function to be minimized. fun is a function that accepts a vector x and returns a vector F, the objective functions evaluated at x.
	//   x0 : a vector of double, contains initial guess of variables.
	//   A : a matrix of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b. 
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
    //  fminimax minimizes the worst-case (largest) value of a set of multivariable functions, starting at an initial estimate. This is generally referred to as the minimax problem.
    //
    //  <latex>
    //  \min_{x} \max_{i} F_{i}(x)\: \textrm{such that} \:\begin{cases}
    //  & c(x) \leq 0 \\
    //  & ceq(x) = 0 \\
    //  & A.x \leq b \\
    //  & Aeq.x = beq \\
    //  & minmaxLb \leq x \leq minmaxUb
    //  \end{cases}
    //  </latex>
    //
    //  Currently, fminimax calls fmincon which uses the ip-opt algorithm.
    //
    //  max-min problems can also be solved with fminimax, using the identity
    //
    //  <latex>
    //  \max_{x} \min_{i} F_{i}(x) = -\min_{x} \max_{i} \left( -F_{i}(x) \right)
    //  </latex>
    //
    // The options allows the user to set various parameters of the Optimization problem. 
    // It should be defined as type "list" and contains the following fields.
    // <itemizedlist>
    //   <listitem>Syntax : options= list("MaxIter", [---], "CpuTime", [---], "GradObj", ---, "GradCon", ---);</listitem>
    //   <listitem>MaxIter : a Scalar, containing the Maximum Number of Iteration that the solver should take.</listitem>
    //   <listitem>CpuTime : a Scalar, containing the Maximum amount of CPU Time that the solver should take.</listitem>
    //   <listitem>GradObj : a function, representing the gradient function of the Objective in Vector Form.</listitem>
    //   <listitem>GradCon : a function, representing the gradient of the Non-Linear Constraints (both Equality and Inequality) of the problem. It is declared in such a way that gradient of non-linear inequality constraints are defined first as a separate Matrix (cg of size m2 X n or as an empty), followed by gradient of non-linear equality constraints as a separate Matrix (ceqg of size m2 X n or as an empty) where m2 & m3 are number of non-linear inequality and equality constraints respectively.</listitem>
    //   <listitem>Default Values : options = list("MaxIter", [3000], "CpuTime", [600]);</listitem>
    // </itemizedlist>
    //
    //  The objective function must have header :
    //  <programlisting>
    //      F = fun(x)
    //  </programlisting>
    //  where x is a n x 1 matrix of doubles and F is a m x 1 matrix of doubles where m is the total number of objective functions inside F.
    //  On input, the variable x contains the current point and, on output, the variable F must contain the objective function values.
    //
    //  By default, the gradient options for fminimax are turned off and and fmincon does the gradient opproximation of minmaxObjfun. In case the GradObj option is off and GradConstr option is on, fminimax approximates minmaxObjfun gradient using numderivative toolbox.
    //
    //  If we can provide exact gradients, we should do so since it improves the convergence speed of the optimization algorithm.
    //
    //  Furthermore, we must enable the "GradObj" option with the statement :
    //  <programlisting>
    //      minimaxOptions = list("GradObj",fGrad);
    //  </programlisting>
    //  This will let fminimax know that the exact gradient of the objective function is known, so that it can change the calling sequence to the objective function. Note that, fGrad should be mentioned in the form of N x n where n is the number of variables, N is the number of functions in objective function. 
    //
    //  The constraint function must have header :
    //  <programlisting>
    //      [c, ceq] = confun(x)
    //  </programlisting>
    //  where x is a n x 1 matrix of dominmaxUbles, c is a 1 x nni matrix of doubles and ceq is a 1 x nne matrix of doubles (nni : number of nonlinear inequality constraints, nne : number of nonlinear equality constraints).
    //  On input, the variable x contains the current point and, on output, the variable c must contain the nonlinear inequality constraints and ceq must contain the nonlinear equality constraints.
    //
    //  By default, the gradient options for fminimax are turned off and and fmincon does the gradient opproximation of confun. In case the GradObj option is on and GradCons option is off, fminimax approximates confun gradient using numderivative toolbox.
    //
    //  If we can provide exact gradients, we should do so since it improves the convergence speed of the optimization algorithm.
    //
    //  Furthermore, we must enable the "GradCon" option with the statement :
    //  <programlisting>
    //      minimaxOptions = list("GradCon",confunGrad);
    //  </programlisting>
    //  This will let fminimax know that the exact gradient of the objective function is known, so that it can change the calling sequence to the objective function.
    //
    //  The constraint derivative function must have header :
    //  <programlisting>
    //      [dc,dceq] = confungrad(x)
    //  </programlisting>
    //  where dc is a nni x n matrix of doubles and dceq is a nne x n matrix of doubles.
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
    //  // A basic case :
    //  // we provide only the objective function and the nonlinear constraint
    //  // function
    //  function y = myfun(x)
    //      y(1)= 2*x(1)^2 + x(2)^2 - 48*x(1) - 40*x(2) + 304;     //Objectives
    //      y(2)= -x(1)^2 - 3*x(2)^2;
    //      y(3)= x(1) + 3*x(2) -18;
    //      y(4)= -x(1) - x(2);
    //      y(5)= x(1) + x(2) - 8;
    //  endfunction
    //  // The initial guess
    //  x0 = [0.1,0.1];
    //  // The expected solution : only 4 digits are guaranteed
    //  xopt = [4 4]
    //  fopt = [0 -64 -2 -8 0]
    //  maxfopt = 0
    //  // Run fminimax
    //  [x,fval,maxfval,exitflag,output,lambda] = fminimax(myfun, x0)
	// // Press ENTER to continue
    //
    // Examples
    //  // A case where we provide the gradient of the objective
    //  // functions and the Jacobian matrix of the constraints.
    //  // The objective function and its gradient
    //  function y = myfun(x)
    //      y(1)= 2*x(1)^2 + x(2)^2 - 48*x(1) - 40*x(2) + 304;
    //      y(2)= -x(1)^2 - 3*x(2)^2;
    //      y(3)= x(1) + 3*x(2) -18;
    //      y(4)= -x(1) - x(2);
    //      y(5)= x(1) + x(2) - 8;
    //  endfunction
    //  // Defining gradient of myfun
    //  function G = myfungrad(x)
    //      G = [ 4*x(1) - 48, -2*x(1), 1, -1, 1;
    //            2*x(2) - 40, -6*x(2), 3, -1, 1; ]'
    //  endfunction
    //  // The nonlinear constraints and the Jacobian
    //  // matrix of the constraints
    //  function [c,ceq] = confun(x)
    //      // Inequality constraints
    //      c = [1.5 + x(1)*x(2) - x(1) - x(2), -x(1)*x(2) - 10] 
    //      // No nonlinear equality constraints
    //      ceq=[]
    //  endfunction
    //  // Defining gradient of confungrad
    //  function [DC,DCeq] = cgrad(x)
    //      // DC(:,i) = gradient of the i-th constraint
    //      // DC = [
    //      //   Dc1/Dx1  Dc1/Dx2
    //      //   Dc2/Dx1  Dc2/Dx2
    //      //   ]
    //      DC= [
    //      x(2)-1, -x(2)
    //      x(1)-1, -x(1)
    //      ]'
    //      DCeq = []'
    //  endfunction
    //  // Test with both gradient of objective and gradient of constraints
    //  minimaxOptions = list("GradObj",myfungrad,"GradCon",cgrad);
    //  // The initial guess
    //  x0 = [0,10];
    //  // The expected solution : only 4 digits are guaranteed
    //  xopt = [0.92791 7.93551]
    //  fopt = [6.73443  -189.778  6.73443  -8.86342  0.86342]
    //  maxfopt = 6.73443
    //  // Run fminimax
    //  [x,fval,maxfval,exitflag,output] = fminimax(myfun,x0,[],[],[],[],[],[], confun, minimaxOptions)
    // Authors
    // Animesh Baranawal
    //

    // Check number of input and output arguments
    [minmaxLhs,minmaxRhs] = argn()
    Checkrhs("fminimax", minmaxRhs, [2 4 6 8 9 10])
    Checklhs("fminimax", minmaxLhs, 1:7)

    // Proper initialisation of objective function
    minmaxObjfun = varargin(1)
    Checktype("fminimax", minmaxObjfun, "minmaxObjfun", 1, "function")

    // Proper initialisation of starting point
    minmaxStartpoint = varargin(2)
    Checktype("fminimax", minmaxStartpoint, "minmaxStartpoint", 2, "constant")

    minmaxNumvar = size(minmaxStartpoint,"*")
    Checkvector("fminimax", minmaxStartpoint, "minmaxStartpoint", 2, minmaxNumvar)
    minmaxStartpoint = minmaxStartpoint(:)

    // Proper initialisation of A and b
    if(minmaxRhs < 3) then // if A and b are not provided, declare as empty
        minmaxA = []
        minmaxB = []
    else
        minmaxA = varargin(3)
        minmaxB = varargin(4)
    end

    Checktype("fminimax", minmaxA, "A", 3, "constant")
    Checktype("fminimax", minmaxB, "b", 4, "constant")

    // Check if A and b of proper dimensions
    if(minmaxA <> [] & minmaxB == []) then
        errmsg = msprintf(gettext("%s: Incompatible input arguments #%d and #%d: matrix A is empty, but the column vector b is not empty"), "fminimax", 3, 4)
        error(errmsg)
    end

    if(minmaxA == [] & minmaxB <> []) then
        errmsg = msprintf(gettext("%s: Incompatible input arguments #%d and #%d: matrix A is not empty, but the column vector b is empty"), "fminimax", 3, 4)
        error(errmsg)
    end

    minmaxNumrowA = size(minmaxA,"r")
    if(minmaxA <> []) then
        Checkdims("fminimax", minmaxA, "A", 3, [minmaxNumrowA minmaxNumvar])
        Checkvector("fminimax", minmaxB, "b", 4, minmaxNumrowA)
        minmaxB = minmaxB(:)
    end

    // Proper initialisation of Aeq and beq
    if(minmaxRhs < 5) then // if Aeq and beq are not provided, declare as empty
        minmaxAeq = []
        minmaxBeq = []
    else
        minmaxAeq = varargin(5)
        minmaxBeq = varargin(6)
    end

    Checktype("fminimax", minmaxAeq, "Aeq", 5, "constant")
    Checktype("fminimax", minmaxBeq, "beq", 6, "constant")

    // Check if Aeq and beq of proper dimensions
    if(minmaxAeq <> [] & minmaxBeq == []) then
        errmsg = msprintf(gettext("%s: Incompatible input arguments #%d and #%d: matrix Aeq is empty, but the column vector beq is not empty"), "fminimax", 5, 6)
        error(errmsg)
    end

    if(minmaxAeq == [] & minmaxBeq <> []) then
        errmsg = msprintf(gettext("%s: Incompatible input arguments #%d and #%d: matrix Aeq is not empty, but the column vector beq is empty"), "fminimax", 5, 6)
        error(errmsg)
    end

    minmaxNumrowAeq = size(minmaxAeq,"r")
    if(minmaxAeq <> []) then
        Checkdims("fminimax", minmaxAeq, "Aeq", 5, [minmaxNumrowAeq minmaxNumvar])
        Checkvector("fminimax", minmaxBeq, "beq", 6, minmaxNumrowAeq)
        minmaxBeq = minmaxBeq(:)
    end

    // Proper initialisation of minmaxLb and minmaxUb
    if(minmaxRhs < 7) then // if minmaxLb and minmaxUb are not provided, declare as empty
        minmaxLb = []
        minmaxUb = []
    else
        minmaxLb = varargin(7)
        minmaxUb = varargin(8)
    end

    Checktype("fminimax", minmaxLb, "lb", 7, "constant")
    Checktype("fminimax", minmaxUb, "ub", 8, "constant")

    // Check dimensions of minmaxLb and minmaxUb
    if(minmaxLb <> []) then
        Checkvector("fminimax", minmaxLb, "lb", 7, minmaxNumvar)
        minmaxLb = minmaxLb(:)
    end

    if(minmaxUb <> []) then
        Checkvector("fminimax", minmaxUb, "ub", 8, minmaxNumvar)
        minmaxUb = minmaxUb(:)
    end

    // Proper Initialisation of minmaxNonlinfun
    if(minmaxRhs < 9) then // if minmaxNonlinfun is not provided, declare as empty
        minmaxNonlinfun = []
    else
        minmaxNonlinfun = varargin(9)
    end

    // fmincon library of scilab gives error when 'c' component of minmaxNonlinfun empty
    // add a trivial case of -5 <= 0 to c to bypass this error
    if(minmaxNonlinfun == []) then
        function [c,ceq] = t(z)
            c = []
            ceq = []
        endfunction
        minmaxNonlinfun = t
    end

    Checktype("fminimax", minmaxNonlinfun, "nonlinfun", 9, "function")

    //To check, Whether minimaxOptions is been entered by user
    if ( minmaxRhs<10 ) then
        minmaxUserOptions = list();
    else
        minmaxUserOptions = varargin(10); //Storing the 3rd Input minmaxUserOptionseter in intermediate list named 'minmaxUserOptions'
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

    //Flags to check whether Gradient is "ON"/"OFF" and store values of user options
    flag1=0;
    flag2=0;
    minmaxMaxIter = 3000
    minmaxCPU = 600
    minmaxFGrad=[];
    minmaxCGrad=[];

    //To check the User Entry for Options and storing it
    for i = 1:(size(minmaxUserOptions))/2
        select convstr(minmaxUserOptions(2*i-1),'l')
        case "maxiter" then
            minmaxIter = minmaxUserOptions(2*i);    //Setting the Maximum Iteration as per user entry

        case "cputime" then
            minmaxCPU = minmaxUserOptions(2*i);    //Setting the Maximum CPU Time as per user entry

        case "gradobj" then
            if (type(minmaxUserOptions(2*i))==10) then
                if (convstr(minmaxUserOptions(2*i))=="off") then
                    flag1 = 0;
                else
                    errmsg = msprintf(gettext("%s: Unrecognized String %s entered for the option- %s."), "fminimax",minmaxUserOptions(2*i), minmaxUserOptions(2*i-1));
                    error(errmsg);
                end
            else
                flag1 = 1;
                minmaxFGrad = minmaxUserOptions(2*i);
            end

        case "gradcon" then
            if (type(minmaxUserOptions(2*i))==10) then
                if (convstr(minmaxUserOptions(2*i))=="off") then
                    flag2 = 0;
                else
                    errmsg = msprintf(gettext("%s: Unrecognized String %s entered for the option- %s."), "fminimax",minmaxUserOptions(2*i), minmaxUserOptions(2*i-1));
                    error(errmsg);
                end
            else
                flag2 = 1;
                minmaxCGrad = minmaxUserOptions(2*i);
            end

        else
            errmsg = msprintf(gettext("%s: Unrecognized minmaxUserOptionseter name ''%s''."), "fminimax", minmaxUserOptions(2*i-1));
            error(errmsg)
        end
    end

    // Checking if minmaxFGrad and minmaxCGrad are functions
    if (flag1==1) then
        if (type(minmaxFGrad) ~= 11 & type(minmaxFGrad) ~= 13) then
            errmsg = msprintf(gettext("%s: Expected function for Gradient of Objective"), "fminimax");
            error(errmsg);
        end
    end

    if (flag2==1) then
        if (type(minmaxCGrad) ~= 11 & type(minmaxCGrad) ~= 13) then
            errmsg = msprintf(gettext("%s: Expected function for Gradient of Nonlinfun"), "fminimax");
            error(errmsg);
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

    // function handle defining new objective function
    function newfunc = newObjfun(z)
        newfunc = z(minmaxNumvar+1)
    endfunction

    // function handle defining add_ineq derivative using numderivative
    function func = minmaxObjDer(z)
        func = numderivative(minmaxAddIneq,z)
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

    // function handle defining new minmaxNonlinfun function
    function [nc,nceq] = newNonlinfun(z)
        [nc,nceq] = minmaxNonlinfun(z)
        // add inequalities of the form Fi(x) - y <= 0
        tmp = [minmaxObjfun(z) - z(minmaxNumvar+1)]'
        nc = [nc, tmp]
    endfunction

    // function handle defining new gradient function for non-linear constraints
    // this function passed when the gradient feature is on
    function [dnc,dnceq] = newCGrad(z)

        // if constraint gradient present use it
        if(flag2 == 1) then
            [dnc, dnceq] = minmaxCGrad(z)
            dnc = [dnc, zeros(size(dnc,'r'),1)]
            dnceq = [dnceq, zeros(size(dnceq,'r'),1)]
        else
            // else use numderivative method to calculate gradient of constraints
            [dnc, dnceq] = minmaxNonlinDer(z)
        end

        // if objective gradient is present use it
        if(flag1 == 1) then
            derObjfun = minmaxFGrad(z)
            mderObjfun = [derObjfun, -1*ones(size(derObjfun,'r'),1)]
            dnc = [dnc; mderObjfun]
        else
            // else use numderivative to calculate gradient of set of obj functions
            derObjfun = minmaxObjDer(z)
            dnc = [dnc; derObjfun]
        end
    endfunction

    //    to be passed as minimaxOptions to fmincon
    if(flag1 == 1 | flag2 == 1) then
        minmaxPassOptions = list("MaxIter", minmaxMaxIter, "CpuTime", minmaxCPU, "GradCon", newCGrad)
        [x,fval,exitflag,output,lambda] = ...
        fmincon(newObjfun,minmaxStartpoint,minmaxA,minmaxB,minmaxAeq,minmaxBeq,minmaxLb,minmaxUb,newNonlinfun,minmaxPassOptions)

        x = x(1:minmaxNumvar)
        fval = minmaxObjfun(x)
        maxfval = max(fval)
    else
        minmaxPassOptions = list("MaxIter", minmaxMaxIter, "CpuTime", minmaxCPU)
        [x,fval,exitflag,output,lambda] = ...
        fmincon(newObjfun,minmaxStartpoint,minmaxA,minmaxB,minmaxAeq,minmaxBeq,minmaxLb,minmaxUb,newNonlinfun,minmaxPassOptions)

        x = x(1:minmaxNumvar)
        fval = minmaxObjfun(x)
        maxfval = max(fval)
    end

endfunction
