// This example has been taken from the MATLAB quadprog page
// http://in.mathworks.com/help/optim/ug/quadprog.html

H = [1 -1; -1 2]; 
f = [-2; -6];
A = [1 1; -1 2; 2 1];
b = [2; 2; 3];
lb=[0,0];
ub=[%inf, %inf];
options=list('MaxIterations',1);

[xopt,fopt,status,output]=intqpipopt(H,f,[],A,b,[],[],lb,ub,[],options)
