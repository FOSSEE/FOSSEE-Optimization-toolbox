
// Check if a user specifies a starting point or not

function y=fun(x)
y=x(1)+x(2);
endfunction

x0 = [];
A = [3,4];
b = [7,9];
options=list("MaxIter", [1500], "CpuTime", [500], "Gradient", "OFF", "Hessian", "OFF");

//  !--error 10000 
// fot_intfmincon: x0 cannot be an empty
// at line     270 of function fot_intfmincon called by :  
// [x0pt,fopt,exitflag,gradient,hessian] = fot_intfmincon(f,x0,intcon, A, b);
// at line      20 of exec file called by :    
// exec fot_intfmincon_x0a.sce


[x0pt,fopt,exitflag,gradient,hessian] = fot_intfmincon(f,x0,intcon, A, b);
