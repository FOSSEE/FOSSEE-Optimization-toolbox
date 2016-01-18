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

[z,gval,attainfactor,exitflag,output,lambda]=fgoalattain(objfun,x0,goal,weight)

