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
