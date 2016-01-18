// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) Scilab Enterprises - 2015- IIT Bombay - FOSSEE 
// Authors: Prajwala TM, Sheetal Shalini
// This file is distributed under the same license as the Scilab package.
//
// assert_close --
//   Returns 1 if the two real matrices computed and expected are close,
//   i.e. if the relative distance between computed and expected is lesser than epsilon.
//   Arguments
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
//  if flag <> 1 then bugmes();quit;end
  flag = assert_checktrue ( shift < epsilon );
endfunction
//
// assert_equal --
//   Returns 1 if the two real matrices computed and expected are equal.
//   Arguments
//   computed, expected : the two matrices to compare
//   epsilon : a small number
//
//function flag = assert_equal ( computed , expected )
//  if computed==expected then
//    flag = 1;
//  else
//    flag = 0;
//  end
//  if flag <> 1 then bugmes();quit;end
//endfunction
//
// Test fgoalattain 
function [f,G] = objfun(x)
f(1)=x(2)^2-4*x(1)
f(2)=x(2)^3-4*x(1)^2+3
f(3)=x(1)-x(2)
G=[-4,-8*x(1),1;2*x(2),3*x(2)^2,-1]
endfunction

function [c,ceq,dc,dceq]=nonlinfun(x)
c=[x(2)+x(1)-4]
ceq=[x(1)]
dc=[1;1]
dceq=[1]
endfunction

goal=[5.0,-6.3,7.2]
weight=[8,2,3]
x0=[-1,2]
A=[1,2]
b=[3]
Aeq=[-1,4]
beq=[5]
lb=[-1,-1]
ub=[10,10]
options=list("MaxIter", [10000], "CpuTime", [5000], "GradObj", "ON","GradCon", "ON");
[z,gval,attainfactor,exitflag,output,lambda]=fgoalattain(objfun,x0,goal,weight,A,b,Aeq,beq,lb,ub,nonlinfun,options)
assert_close ( x , [ 0.    1.25], 1e-4 )
assert_close ( fval , [1.5625 , 4.953125 , - 1.25 ] , 1e-4 );
assert_close (attainfactor, 5.6265625 ,1e-4);
