// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Authors: Prajwala TM,Sheetal Shalini
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

function [x,fval,attainfactor,exitflag,output,lambda] = fgoalattain(varargin)
    // Solves a multiobjective goal attainment problem
    //
    // Calling Sequence
    //  xopt = fgoalattain(fun,x0,goal,weight)
    //  xopt = fgoalattain(fun,x0,goal,weight,A,b)
    //  xopt = fgoalattain(fun,x0,goal,weight,A,b,Aeq,beq)
    //  xopt = fgoalattain(fun,x0,goal,weight,A,b,Aeq,beq,lb,ub)
    //  xopt = fgoalattain(fun,x0,goal,weight,A,b,Aeq,beq,lb,ub,nonlcon)
    //  xopt = fgoalattain(fun,x0,goal,weight,A,b,Aeq,beq,lb,ub,nonlcon,options)
    //  [xopt,fval] = fgoalattain(...)
    //  [xopt,fval,attainfactor] = fgoalattain(...)
    //  [xopt,fval,attainfactor,exitflag] = fgoalattain(...)
    //  [xopt,fval,attainfactor,exitflag,output] = fgoalattain(...)
    //  [xopt,fval,attainfactor,exitflag,output,lambda] = fgoalattain(...)
    //
    // Parameters
    //  fun: a function that accepts a vector x and returns a vector F
	//   x0 : a vector of double, contains initial guess of variables.
	//   A : a matrix of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b. 
	//   b : a vector of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b.
	//   Aeq : a matrix of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   beq : a vector of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   lb : a vector of double, contains lower bounds of the variables.
	//   ub : a vector of double,  contains upper bounds of the variables.
    //   nonlcon: a function, the nonlinear constraints
    //   options : a list, containing the option for user to specify. See below for details.
	//   xopt : a vector of double, the computed solution of the optimization problem.
	//   fopt : a double, the value of the function at x.
    //   attainfactor: The amount of over- or underachievement of the goals,γ at the solution.
	//   exitflag : The exit status. See below for details.
	//   output : The structure consist of statistics about the optimization. See below for details.
	//   lambda : The structure consist of the Lagrange multipliers at the solution of problem. See below for details.
    //
    // Description
    // fgoalattain solves the goal attainment problem, which is one formulation for minimizing a multiobjective optimization problem.
    // Finds the minimum of a problem specified by:
    // Minimise Y such that
    //
    //<latex>
    //\begin{eqnarray}
    //\mbox{min}_{x,\gamma}  & f(x)-weight \ast \gamma \leq goal \\
    //\mbox{subject to} & c(x) \leq 0 \\
    //                  & c_{eq}(x) = 0 \\
    //                  & Ax \leq b \\
    //                  & A_{eq} x = b_{eq} \\
    //                  & lb \leq x \leq ub
    //\end{eqnarray}
    //</latex>
    //
    // The solver makes use of fmincon to find the minimum.
    //
    // The fgoalattain finds out the maximum value of Y for the objectives evaluated at the starting point and
    // adds that as another variable to the vector x
    // This is passed to the fmincon function to get the optimised value of Y
    // Hence, the algorithm used mainly is "ipopt" to obtain the optimum solution
    // The relations between f(x), Y, weights and goals are added as additional non-linear inequality constraints
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
    //  By default, the gradient options for fminimax are turned off and and fmincon does the gradient opproximation of gattainObjfun. In case the GradObj option is off and GradConstr option is on, fminimax approximates gattainObjfun gradient using numderivative toolbox.
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
    //  where x is a n x 1 matrix of doubles, c is a 1 x nni matrix of doubles and ceq is a 1 x nne matrix of doubles (nni : number of nonlinear inequality constraints, nne : number of nonlinear equality constraints).
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
    //  function f1 = gattainObjfun(x)
    //      f1(1)=2*x(1)*x(1)+x(2)*x(2)-48*x(1)-40*x(2)+304
    //      f1(2)=-x(1)*x(1)-3*x(2)*x(2)
    //      f1(3)=x(1)+3*x(2)-18
    //      f1(4)=-x(1)-x(2)
    //      f1(5)=x(1)+x(2)-8
    //  endfunction
    //  x0=[-1,1];
    //  goal=[-5,-3,-2,-1,-4];
    //  weight=abs(goal)
    //  //gval  =[- 0.0000011 -63.999998 -2.0000002 -8 3.485D-08]
    //  //z  = [4 3.99]
    //  //Run fgoalattain
    //  [x,fval,attainfactor,exitflag,output,lambda]=fgoalattain(gattainObjfun,x0,goal,weight)
    // Authors
    // Prajwala TM, Sheetal Shalini , 2015

    // Check number of input and output arguments
    [gattainLhs,gattainRhs] = argn()
    fot_Checkrhs("fgoalattain", gattainRhs, [4 6 8 10 11 12])
    fot_Checklhs("fgoalattain", gattainLhs, 1:6)

    // initialisation of fun
    gattainObjfun = varargin(1)
    fot_Checktype("fgoalattain", gattainObjfun, "fun", 1, "function")

    // initialisation of x0
    gattainStartpoint = varargin(2)
    fot_Checktype("fgoalattain", gattainStartpoint, "x0", 2, "constant")

    gattainNumvar = size(gattainStartpoint,"*")
    gattainStartpoint = gattainStartpoint(:)

    // initialisation of goal
    goal=varargin(3)
    fot_Checktype("fgoalattain",goal,"goal",3,"constant")

    // initialisation of weight
    weight=varargin(4)
    fot_Checktype("fgoalattain",weight,"weight",4,"constant")

    //initialisation of A and b
    if(gattainRhs < 5) then
        gattainA = []
        gattainB = []
    else
        gattainA = varargin(5)
        gattainB = varargin(6)
    end

    fot_Checktype("fgoalattain", gattainA, "A", 5, "constant")
    fot_Checktype("fgoalattain", gattainB, "b", 6, "constant")

    gattainNumrowA = size(gattainA,"r")
    if(gattainA <> []) then
        fot_Checkdims("fgoalattain", gattainA, "A", 5, [gattainNumrowA gattainNumvar])
        fot_Checkvector("fgoalattain", gattainB, "b", 6, gattainNumrowA)
        gattainB = gattainB(:)
    end

    //initialisation of Aeq and beq
    if(gattainRhs < 7) then
        gattainAeq = []
        gattainBeq = []
    else
        gattainAeq = varargin(7)
        gattainBeq = varargin(8)
    end

    fot_Checktype("fgoalattain", gattainAeq, "Aeq", 7, "constant")
    fot_Checktype("fgoalattain", gattainBeq, "beq", 8, "constant")

    gattainNumrowAeq = size(gattainAeq,"r")
    if(gattainAeq <> []) then
        fot_Checkdims("fgoalattain", gattainAeq, "Aeq", 7, [gattainNumrowAeq gattainNumvar])
        fot_Checkvector("fgoalattain", gattainBeq, "beq", 8, gattainNumrowAeq)
        gattainBeq = gattainBeq(:)
    end

    // initialisation of lb and ub
    if(gattainRhs < 9) then
        gattainLb = []
        gattainUb = []
    else
        gattainLb = varargin(9)
        gattainUb = varargin(10)
    end

    fot_Checktype("fgoalattain", gattainLb, "lb", 9, "constant")
    fot_Checktype("fgoalattain", gattainUb, "ub", 10, "constant")

    // Check dimensions of lb and ub
    if(gattainLb <> []) then
        fot_Checkvector("fgoalattain", gattainLb, "lb", 9, gattainNumvar)
        gattainLb = gattainLb(:)
    end

    if(gattainUb <> []) then
        fot_Checkvector("fgoalattain", gattainUb, "ub", 10, gattainNumvar)
        gattainUb = gattainUb(:)
    end

    // Initialisation of nonlcon
    function [c,ceq] = constr(z)
        c = []
        ceq = []
    endfunction

    if(gattainRhs < 11) then
        gattainNonlinfun = constr
    else
        gattainNonlinfun = varargin(11)
    end

    fot_Checktype("fgoalattain", gattainNonlinfun, "nonlcon", 11, "function")

    // initialisation of default options
    if(gattainRhs < 12) then
        gattainUserOptions = list()
    else
        gattainUserOptions = varargin(12)
    end

    //If gattainOptions is entered then checking its type for 'list'
    if (type(gattainUserOptions) ~= 15) then
        errmsg = msprintf(gettext("%s: gattainOptions (10th parameter) should be a list"), "fgoalattain");
        error(errmsg);
    end

    //If minimaxOptions is entered then checking whether even number of entires are entered
    if (modulo(size(gattainUserOptions),2)) then
        errmsg = msprintf(gettext("%s: Size of gattainOptions (list) should be even"), "fgoalattain");
        error(errmsg);
    end

    //Flags to check whether Gradient is "ON"/"OFF" & Hessian is "ON"/"OFF"
    flag1=0;
    flag2=0;
    fgaMaxIter = 3000;
    fgaCPU = 600;
    gattainFGrad=[];
    gattainCGrad=[];

    //To check the User Entry for Options and storing it
    for i = 1:(size(gattainUserOptions))/2
        select convstr(gattainUserOptions(2*i-1),'l')
        case "maxiter" then
            fgaMaxIter = gattainUserOptions(2*i);    //Setting the Maximum Iteration as per user entry

        case "cputime" then
            fgaCPU = gattainUserOptions(2*i);    //Setting the Maximum CPU Time as per user entry

        case "gradobj" then
            if (type(gattainUserOptions(2*i))==10) then
                if (convstr(gattainUserOptions(2*i))=="off") then
                    flag1 = 0;
                else
                    errmsg = msprintf(gettext("%s: Unrecognized String %s entered for the option- %s."), "fgoalattain",gattainUserOptions(2*i), gattainUserOptions(2*i-1));
                    error(errmsg);
                end
            else
                flag1 = 1;
                gattainFGrad = gattainUserOptions(2*i);
            end

        case "gradcon" then
            if (type(gattainUserOptions(2*i))==10) then
                if (convstr(gattainUserOptions(2*i))=="off") then
                    flag2 = 0;
                else
                    errmsg = msprintf(gettext("%s: Unrecognized String %s entered for the option- %s."), "fgoalattain",gattainUserOptions(2*i), gattainUserOptions(2*i-1));
                    error(errmsg);
                end
            else
                flag2 = 1;
                gattainCGrad = gattainUserOptions(2*i);
            end

        else
            errmsg = msprintf(gettext("%s: Unrecognized gattainUserOptionseter name ''%s''."), "fminimax", gattainUserOptions(2*i-1));
            error(errmsg)
        end
    end

    // Checking if gattainFGrad and gattainCGrad are functions
    if (flag1==1) then
        if (type(gattainFGrad) ~= 11 & type(gattainFGrad) ~= 13) then
            errmsg = msprintf(gettext("%s: Expected function for Gradient of Objective"), "fminimax");
            error(errmsg);
        end
    end

    if (flag2==1) then
        if (type(gattainCGrad) ~= 11 & type(gattainCGrad) ~= 13) then
            errmsg = msprintf(gettext("%s: Expected function for Gradient of Nonlinfun"), "fminimax");
            error(errmsg);
        end
    end

    gattainObjfunval = gattainObjfun(gattainStartpoint)
    gattainObjfunval=gattainObjfunval(:)
    goal=goal(:)
    weight=weight(:)
    gaVal=[]

    // appending the gamma value as another variable
    for i=1:size(gattainObjfunval,'r')
        if(weight(i)<>0) then
            gaVal(i)=((gattainObjfunval(i)-goal(i))/weight(i))
        end
    end

    gattainStartpoint(gattainNumvar+1)=max(gaVal)

    if(gattainA <> []) then
        gattainA = [gattainA, zeros(gattainNumrowA,1)]
    end
    if(gattainAeq <> []) then
        gattainAeq = [gattainAeq, zeros(gattainNumrowAeq,1)]
    end
    if(gattainLb <> []) then
        gattainLb(gattainNumvar+1) = -%inf
    end
    if(gattainUb <> []) then
        gattainUb(gattainNumvar+1) = +%inf
    end

    // function handle defining the additional inequalities
    function temp = gattainAddIneq(z)
        gaVar = gattainObjfun(z)
        gattainAddIneqWithWt = []
        gattainAddIneqWitoutWt = []
        for i = 1:size(gaVar,'r')
            if(weight(i) <> 0) then
                gattainAddIneqWithWt = [gattainAddIneqWithWt; ( (gaVar(i)-goal(i))/weight(i) )]
            else
                gattainAddIneqWitoutWt = [gattainAddIneqWitoutWt; gaVar(i)-goal(i)]
            end
        end
        temp = [gattainAddIneqWithWt - z(gattainNumvar+1); gattainAddIneqWitoutWt]
    endfunction

    // function handle defining new objective function
    function newfunc = newObjfun(z)
        newfunc = z(gattainNumvar+1)
    endfunction

    // function handle defining add_ineq derivative using numderivative
    function func = gattainIneqDer(z)
        func = numderivative(gattainAddIneq,z)
    endfunction

    // function handle defining nonlcon derivative using numderivative
    function [dc,dceq] = gattainNonlinDer(z)

        // function handle extracting c and ceq components from nonlcon
        function foo = gattainC(z)
            [foo,tmp1] = gattainNonlinfun(z)
            foo = foo'
        endfunction

        function foo = gattainCEQ(z)
            [tmp1,foo] = gattainNonlinfun(z)
            foo = foo'
        endfunction

        dc = numderivative(gattainC,z)
        dceq = numderivative(gattainCEQ,z)
    endfunction

    // function handle defining new nonlcon function
    function [nc,nceq] = newNonlinfun(z)
        [nc,nceq] = gattainNonlinfun(z)
        tmp = [gattainAddIneq(z)]'
        nc = [nc, tmp]
    endfunction

    function [dnc,dnceq] = newCGrad(z)
        // check if "GradCon" option is turned on
        // if "GradCon" is turned on, use it
        if(flag2 == 1) then
            [dnc,dnceq] = gattainCGrad(z)
            dnc = [dnc, zeros(size(dnc,'r'),1)]
            dnceq = [dnceq, zeros(size(dnceq,'r'),1)]
            // else, calculate it using finite differences
        else
            [dnc,dnceq] = gattainNonlinDer(z)
        end

        // check if "GradObj" option is turned on
        // if "GradObj" is turned on, use it
        if(flag1 == 1) then
            derObjfun = gattainFGrad(z)
            tmp1 = []
            tmp2 = []
            for i = 1:size(gattainObjfun(gattainStartpoint),'r')
                if weight(i) <> 0 then
                    gaVal = [derObjfun(i,:)/weight(i) , -1]
                    tmp1 = [ tmp1; gaVal ]
                else
                    gaVal = [derObjfun(i,:) , 0]
                    tmp2 = [ tmp2; gaVal ]
                end
            end

            dnc = [ dnc; tmp1; tmp2 ]
            // else, calculate it using finite differences
        else
            deraddineq = gattainIneqDer(z)
            dnc = [dnc; deraddineq]
        end
    endfunction

    //to be passed as options to fmincon
    if (flag1 == 1 | flag2 == 1) then
        gattainPassOptions = list("MaxIter", fgaMaxIter, "CpuTime", fgaCPU, "GradCon", newCGrad)
        [x,attainfactor,exitflag,output,lambda] = fmincon(newObjfun,gattainStartpoint,gattainA,gattainB,gattainAeq,gattainBeq,gattainLb,gattainUb,newNonlinfun,gattainPassOptions)
        x= x(1:gattainNumvar)
        fval = gattainObjfun(x)
    else
        gattainPassOptions = list("MaxIter", fgaMaxIter, "CpuTime", fgaCPU)
        [x,attainfactor,exitflag,output,lambda] = fmincon(newObjfun,gattainStartpoint,gattainA,gattainB,gattainAeq,gattainBeq,gattainLb,gattainUb,newNonlinfun,gattainPassOptions)
        x= x(1:gattainNumvar)
        fval = gattainObjfun(x)
    end

endfunction
