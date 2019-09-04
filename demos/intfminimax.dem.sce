mode(1)
//
// Demo of intfminimax.sci
//

// A basic case :
// we provide only the objective function and the nonlinear constraint
// function
function y = myfun(x)
y(1)= 2*x(1)^2 + x(2)^2 - 48*x(1) - 40*x(2) + 304;     //Objectives
y(2)= -x(1)^2 - 3*x(2)^2;
y(3)= x(1) + 3*x(2) -18;
y(4)= -x(1) - x(2);
y(5)= x(1) + x(2) - 8;
endfunction
// The initial guess
x0 = [0.1,0.1];
// The expected solution : only 4 digits are guaranteed
xopt = [4 4]
fopt = [0 -64 -2 -8 0]
intcon = [1]
maxfopt = 0
// Run fminimax
[x,fval,maxfval,exitflag] = intfminimax(myfun, x0,intcon)
// Press ENTER to continue
halt()   // Press return to continue
 
// A case where we provide the gradient of the objective
// functions and the Jacobian matrix of the constraints.
// The objective function and its gradient
function [y,G] = myfun(x)
y(1)= 2*x(1)^2 + x(2)^2 - 48*x(1) - 40*x(2) + 304;
y(2)= -x(1)^2 - 3*x(2)^2;
y(3)= x(1) + 3*x(2) -18;
y(4)= -x(1) - x(2);
y(5)= x(1) + x(2) - 8;
G = [ 4*x(1) - 48, -2*x(1), 1, -1, 1;
2*x(2) - 40, -6*x(2), 3, -1, 1; ]'
endfunction
// The nonlinear constraints
function [c,ceq,DC,DCeq] = confun(x)
// Inequality constraints
c = [1.5 + x(1)*x(2) - x(1) - x(2), -x(1)*x(2) - 10]
// No nonlinear equality constraints
ceq=[]
DC= [x(2)-1, -x(2);
x(1)-1, -x(1)]'
DCeq = []'
endfunction
// Test with both gradient of objective and gradient of constraints
minimaxOptions = list("GradObj","on","GradCon","on");
// The initial guess
x0 = [0,10];
intcon = [2]
// Run intfminimax
[x,fval,maxfval,exitflag] = intfminimax(myfun,x0,intcon,[],[],[],[],[],[], confun, minimaxOptions)
//========= E N D === O F === D E M O =========//
