function f3=objfun3(x)
f3(1)=x(1)+x(2)-9*x(1)*x(2)
f3(2)=x(3)-87*x(4)+x(1)
f3(3)=x(5)-6*x(4)*x(3)
endfunction

x0=[1,3,2,5,4];
A=[-3,7,6,1,2];
b=[9];
Aeq=[1,2,3,4,5];
beq=[3];
goal=[-6,8,-2];
weight=[7,0,3];

[z,gval,attainfactor,exitflag,output,lambda]=fgoalattain(objfun3,x0,goal,weight,A,b,Aeq,beq)
