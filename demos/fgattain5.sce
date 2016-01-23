function f1 = objfun(x)
    f1(1)=2*x(1)*x(1)+x(2)*x(2)-48*x(1)-40*x(2)+304
    f1(2)=-x(1)*x(1)-3*x(2)*x(2)
    f1(3)=x(1)+3*x(2)-18
    f1(4)=-x(1)-x(2)
    f1(5)=x(1)+x(2)-8
endfunction

x0=[-1,1];
goal=[-5,-3,-2,-1,-4];
weight=abs(goal);

A=[1,2]
b=[3]
Aeq=[3,4]
beq=[4]
lb=[6,3];
ub=[7,4];

function [c,ceq]=nonlinfun(x)
c=[x(2)+x(1)-4]
ceq=[4*x(1)-2]
endfunction

[z,gval,attainfactor,exitflag,output,lambda]=fgoalattain(objfun,x0,goal,weight,A,b,Aeq,beq,lb,ub,nonlinfun)
