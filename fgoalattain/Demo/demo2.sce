function f2=objfun1(x)
f2(1) = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
f2(2)= x(2)-x(1)*5+x(2)*x(2)
endfunction
x0 = [-1,2];
A = [1,2];
b = [3];
goal=[5,-6];
weight=abs(goal);
[z,gval,attainfactor,exitflag,output,lambda]=fgoalattain(objfun1,x0,goal,weight,A,b)

