//Example 6:
//Infeasible objective function.
function y=f(x)
y=x(1)^2 - x(1)*x(2)/3 + x(2)^2;
endfunction
x0=[0 , 0];
A=[1,1 ; 1,1/4 ; 1,-1];
b=[2;1;1];
Aeq = [1,1];
beq = 3;
[x,fval,exitflag,output,lambda,grad,hessian] =fmincon(f, x0,A,b,Aeq,beq)
