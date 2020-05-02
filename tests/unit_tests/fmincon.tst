// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// Author: R.Vidyadhar & Vignesh Kannan
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

// <-- JVM NOT MANDATORY -->
// <-- ENGLISH IMPOSED -->


//
// assert_close --
//   Returns 1 if the two real matrices computed and expected are close,
//   i.e. if the relative distance between computed and expected is lesser than epsilon.
// Arguments
//   computed, expected : the two matrices to compare
//   epsilon : a small number
//
function flag = assert_close ( computed, expected, epsilon )
  if expected==0.0 then
    shift = norm(computed-expected);
  else
    shift = norm(computed-expected)/norm(expected);
  end
//  if shift < epsilon then
//    flag = 1;
//  else
//    flag = 0;
//  end
//  if flag <> 1 then pause,end
    flag = assert_checktrue ( shift < epsilon );
endfunction
//
// assert_equal --
//   Returns 1 if the two real matrices computed and expected are equal.
// Arguments
//   computed, expected : the two matrices to compare
//   epsilon : a small number
//
//function flag = assert_equal ( computed , expected )
//  if computed==expected then
//    flag = 1;
//  else
//    flag = 0;
//  end
//  if flag <> 1 then pause,end
//endfunction

//Find x in R^3 such that it minimizes:
//f(x)= x1*x2 + x2*x3
//x0=[0.1 , 0.1 , 0.1]
//constraint-1 (c1): x1^2 - x2^2 + x3^2 <= 2
//constraint-2 (c2): x1^2 + x2^2 + x3^2 <= 10
//Objective function to be minimised
function y=f(x)
y=x(1)*x(2)+x(2)*x(3);
endfunction
//Starting point, linear constraints and variable bounds
x0=[0.1 , 0.1 , 0.1];
A=[];
b=[];
Aeq=[];
beq=[];
lb=[];
ub=[];
//Nonlinear constraints
function [c,ceq]=nlc(x)
c = [x(1)^2 - x(2)^2 + x(3)^2 - 2 , x(1)^2 + x(2)^2 + x(3)^2 - 10];
ceq = [];
endfunction
//Gradient of objective function
function y= fGrad(x)
y= [x(2),x(1)+x(3),x(2)];
endfunction
//Hessian of the Lagrange Function
function y= lHess(x,obj,lambda)
y= obj*[0,1,0;1,0,1;0,1,0] + lambda(1)*[2,0,0;0,-2,0;0,0,2] + lambda(2)*[2,0,0;0,2,0;0,0,2]
endfunction
//Gradient of Non-Linear Constraints
function [cg,ceqg] = cGrad(x)
cg=[2*x(1) , -2*x(2) , 2*x(3) ; 2*x(1) , 2*x(2) , 2*x(3)];
ceqg=[];
endfunction
//Options
options=list("MaxIter", [1500], "CpuTime", [500], "GradObj", fGrad, "Hessian", lHess,"GradCon", cGrad);
//Calling Ipopt
[xopt,fval,exitflag,output] =fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options)

assert_close ( xopt , [ 1.5811388 -2.236068 1.5811388 ]' , 0.0005 );
assert_close ( fval , [-7.0710678 ]' , 0.0005 );
assert_checkequal( exitflag , int32(0) );
printf("Test Successful");
