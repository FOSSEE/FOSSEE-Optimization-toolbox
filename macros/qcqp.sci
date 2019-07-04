// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Adarsh Shah
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

function [xopt,fopt,lambda,exitflag] = qcqp (varargin)
	// Solves a Linear and Quadratic Constrainted Quadratic problem.
	//
	//   Syntax
	//   xopt = qcqp(x,H,f)
    //   xopt = qcqp(x,H,f,Q,c,r,A,b)
    //   xopt = qcqp(x,H,f,Q,c,r,A,b,Aeq,beq)
	//   xopt = qcqp(x,H,f,Q,c,r,A,b,Aeq,beq,lb,ub)
	//   [xopt,fopt,lambda,exitflag] = qcqp( ... )
	//   
    //   Parameters
    //   x : a matrix of double, represents initial point.
	//   H : a symmetric matrix of double, represents coefficients of quadratic terms in the quadratic problem.
    //   f : a vector of double, represents coefficients of linear terms in the quadratic problem
    //   Q : a [n x n] . q matrix of double, represents coefficients of quadratic terms in the quadratic constraints. x'.Q.x + c'.x <= r
    //   c : a q x n matrix of double, represents coefficients of  linear terms in the quadratic problem. x'.Q.x + c'.x <= r
    //   r : a vector of double, represents the linear constants in the inequality constraints x'.Q.x + c'.x <= r.	
	//   A : a matrix of double, represents the linear coefficients in the inequality constraints A⋅x ≤ b. 
	//   b : a vector of double, represents the linear terms in the inequality constraints A⋅x ≤ b.
	//   Aeq : a matrix of double, represents the linear coefficients in the equality constraints Aeq⋅x = beq.
	//   beq : a vector of double, represents the linear terms in the equality constraints Aeq⋅x = beq.
	//   lb : a vector of double, contains lower bounds of the variables. The default value is 0.
	//   ub : a vector of double, contains upper bounds of the variables. The default value is %inf.
	//   xopt : a vector of double, the computed solution of the optimization problem.
    //   fopt : a double, the value of the function at x.
    //   lambda : a vector of double, final estimation of the Lagrange multipliers
    //   exitflag : a double, this output parameter tells what happened in this subroutine, according to the following conventions. 
    //
	//   Description
    //   Search the minimum of a linear and quadratic constrained quadratic optimization problem specified by :
    //
    // Min x’.H.x + f’.x
    //
    // Subject to
    //
    // A’x <= b
    //
    // Aeq x = beq
    //
    // x’.Q.x + c’x <= r
    //
    // lb <= x <= ub
	//   
	//   The routine calls ALGENCAN for solving the quadratic problem, ALGENCAN is a library written in fortran.
	//
	// For more details on exitflag see the CLP documentation, go to https://www.ime.usp.br/~egbirgin/tango/
    //
    // <itemizedlist>
	//   <listitem> exitflag=0 : Optimal Solution Found </listitem>
	//   <listitem> exitflag=1 : Maximum Number of Output Iterations Exceeded. Output may not be optimal </listitem>
	//   <listitem> exitflag=2 : Maximum Number of Total number of Inner Iterations Exceeded. Output may not be optimal </listitem>
	//   <listitem> exitflag=3 : Maximum Number of Total number of Functional Evaluations Exceeded. Output may not be optimal </listitem>
	//   <listitem> exitflag=4 : The algorithm stopped by ``lack of feasibility progress'', i.e., the current point is infeasible </listitem>
    // </itemizedlist>  
    //
    // Examples
    // //Reference : http://www.minlplib.org/nvs10.html
    // // var i1 := 1, >= 0, <= 200;
    // // var i2 := 1, >= 0, <= 200;
    // // minimize obj: 7*i1^2 + 6*i2^2 - 35*i1 - 80.4*i2;
    // // subject to
    // // e1: (-9*i1^2) - 10*i1*i2 - 8*i2^2 >= -583;
    // // e2: (-6*i1^2) - 8*i1*i2 - 6*i2^2 >= -441;
    // H = [7 0;0 6];
    // f = [-35 ; -80.4];
    // Q = [[9 5;5 8];[3 4;4 3]];
    // c = [0 0;0 0];
    // r = [583;441];
    // x = [1;1];
    // lb = [0;0];
    // ub = [200;200];
    // [xopt,fopt] = qcqp(x,H,f,Q,c,r,[],[],[],[],lb,ub);
	// // Press ENTER to continue 
    //
    //  Examples
    // // var i1 integer := 1, >= 0, <= 200;
    // // var i2 integer := 1, >= 0, <= 200;
    // // var i3 integer := 1, >= 0, <= 200;
    // // minimize obj: 7*i1^2 + 6*i2^2 - 15.8*i1 - 93.2*i2 + 8*i3^2 - 6*i3*i1 + 4*i3*i2
    // //      - 63*i3;
    // // subject to
    // // e1: (-9*i1^2) - 10*i1*i2 - 8*i2^2 - 5*i3^2 - 6*i3*i1 - 10*i3*i2 >= -1000;
    // // e2: (-6*i1^2) - 8*i1*i2 - 6*i2^2 - 4*i3^2 - 2*i3*i1 - 2*i3*i2 >= -550;
    // // e3: (-9*i1^2) - 6*i2^2 - 8*i3^2 + 2*i2*i1 + 2*i3*i2 >= -340;
    // x = [1 ; 1; 1];
    // lb = [0 ; 0; 0];
    // ub = [200 ; 200 ; 200];
    // H = [7 0 -3; 
    //       0 6 2; 
    //      -3 2 8];
    // f = [-15.8 ; -93.2 ; -63];
    //  Q = [
    //  [9 5 3;
    //   5 8 5;
    //   3 5 8];
    //  [3 4 1;
    //   4 6 1;
    //   1 1 4];
    //  [9 -1 0;
    //   -1 6 -1;
    //   0 -1 8];
    //         ];
    //  c = zeros(3,3);
    //  r = [1000;550;340];
    //  [xopt,fopt,lambda,exitflag] = qcqp(x,H,f,Q,c,r,[],[],[],[],lb,ub);
    // // Press Enter to continue
    //
	// Authors
	// Adarsh Shah
    

    [lhs,rhs] = argn();
    if (rhs < 3 || rhs < 6 || rhs < 8 || rhs < 10 || rhs > 12 ) then
        errmsg = msprintf(gettext("%s: Unexpected number of input arguments : %d provided while should be in the set of [3 6 8 10 12]"), "qcqp", rhs);
		error(errmsg);
    end
	x = varargin(1);
    H = varargin(2);
    f = varargin(3);
    
    n = size(x,1);
    
    [t1,t2] = size(H)
    if (isequal(t1,n) == 0 || isequal(t2,n) == 0) then
        errmsg = msprintf(gettext("%s: Invalid Input %d arguement : a n x n square matrix expected"), "qcqp", 2);
		error(errmsg);
    end
    
    [t1,t2] = size(f)
    if (isequal(t1,n) == 0 || isequal(t2,1) == 0) then
        errmsg = msprintf(gettext("%s: Invalid Input %d arguement : a n x 1 matrix expected"), "qcqp", 3);
		error(errmsg);
    end
   
    if(rhs < 6) then
        Q = varargin(4);
        c = varargin(5);
        r = varargin(6);    
        
        [q,t1] = size(c)
        if (isequal(t1,n) == 0) then
            errmsg = msprintf(gettext("%s: Invalid Input %d arguement : a q x n 2 D matrix expected"), "qcqp", 5);
            error(errmsg);
        end
        
        [t1,t2] = size(Q)
        if (isequal(t1,n) == 0 || isequal(t2,q*n)) then
            errmsg = msprintf(gettext("%s: Invalid Input %d arguement : a n x q.n matrix expected"), "qcqp", 4);
            error(errmsg);
        end
        
        [t1,t2] = size(r)
        if (isequal(t1,q) == 0 || isequal(t2,1) == 0) then
            errmsg = msprintf(gettext("%s: Invalid Input %d arguement : a q x 1 matrix expected"), "qcqp", 6);
            error(errmsg);
        end        
    else
        Q = [];
        c = [];
        r = [];
    end
    
    if(rhs < 8) then
        A = varargin(7);
        b = varargin(8);
        
        [m,t1] = size(A);
        if (isequal(t1,n) == 0) then
            errmsg = msprintf(gettext("%s: Invalid Input %d arguement : a m x n 2 D matrix expected"), "qcqp", 7);
            error(errmsg);
        end
        
        [t1,t2] = size(b);
        if (isequal(t1,m) == 0 || isequal(t2,1) == 0) then
            errmsg = msprintf(gettext("%s: Invalid Input %d arguement : a m x 1 matrix expected"), "qcqp", 8);
            error(errmsg);
        end
    else
        A = [];
        b = [];
    end
    
    if(rhs < 10) then
        Aeq = varargin(9);
        beq = varargin(10);
        
        [p,t1] = size(Aeq);
        if (isequal(t1,n) == 0) then
            errmsg = msprintf(gettext("%s: Invalid Input %d arguement : a p x n 2 D matrix expected"), "qcqp", 9);
            error(errmsg);
        end
        
        [t1,t2] = size(b);
        if (isequal(t1,p) == 0 || isequal(t2,1) == 0) then
            errmsg = msprintf(gettext("%s: Invalid Input %d arguement : a p x 1 matrix expected"), "qcqp", 10);
            error(errmsg);
        end
    else
        Aeq = []
        beq = []
    end
       
    if(rhs < 11) then
        lb = varargin(11);
        [t1,t2] = size(lb);
        if (isequal(t1,n) == 0 || isequal(t2,1) == 0) then
            errmsg = msprintf(gettext("%s: Invalid Input %d arguement : a n x 1 matrix expected"), "qcqp", 11);
            error(errmsg);
        end
    else
        lb = zeros(n,1);
    end
    
    if(rhs < 12) then
        ub = varargin(12);
        [t1,t2] = size(ub);
        if (isequal(t1,n) == 0 || isequal(t2,1) == 0) then
            errmsg = msprintf(gettext("%s: Invalid Input %d arguement : a n x 1 matrix expected"), "qcqp", 12);
            error(errmsg);
        end        
    else
        ub = ones(n,1)*%inf;
    end
    
    [xopt,fopt,lambda,exitflag] = sci_qcqp(x,H,f,A,b,Aeq,beq,Q,c,r,lb,ub);

    select exitflag
    case 0 then
        disp('Optimal Solution Found');
    case 1 then
        disp('Maximum Number of Output Iterations Exceeded. Output may not be optimal');
    case 2 then
        disp('Maximum Number of Total number of Inner Iterations Exceeded. Output may not be optimal');
    case 3 then
        disp('Maximum Number of Total number of Functional Evaluations Exceeded. Output may not be optimal');
    case 4 then
        disp('The algorithm stopped by ``lack of feasibility progress'', i.e., the current point is infeasible');           
    end

endfunction
