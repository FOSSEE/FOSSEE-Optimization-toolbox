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
function f1 = objfun(x)
    f1(1)=2*x(1)*x(1)+x(2)*x(2)-48*x(1)-40*x(2)+304
    f1(2)=-x(1)*x(1)-3*x(2)*x(2)
    f1(3)=x(1)+3*x(2)-18
    f1(4)=-x(1)-x(2)
    f1(5)=x(1)+x(2)-8
endfunction

goal=[-5,-3,-2,-1,-4];
weight=abs(goal)
x0=[-1,1];
[z,gval,attainfactor,exitflag,output,lambda]=fgoalattain(objfun,x0,goal,weight)

assert_close ( x , [4.0   3.99], 1e-4 );
assert_close ( fval ,[- 1.229D-9  - 64. - 2. - 8. 3.842D-11  ] , 1e-4 );
assert_checkequal ( attainfactor , 1.0 );
assert_checkequal ( exitflag , 0 );
assert_checkequal ( output.iterations > 0 , %t );
assert_checkequal ( output.funcCount > 0 , %t );
assert_checkequal ( output.constrviolation <> [] , %t );
assert_checkequal ( output.algorithm , "ipopt" );
assert_checkequal ( output.firstorderopt <> [] , %t );
assert_checkequal ( output.cgiterations , 0 );
assert_checkequal ( output.message(1) , "Optimal solution found.");
assert_checkequal ( lambda.eqlin , [] );
assert_checkequal ( lambda.eqnonlin , [] );
assert_checkequal ( lambda.ineqlin , [] );
assert_checkequal ( lambda.lower , [0 0]' );
assert_checkequal ( lambda.upper , [0 0]' );

