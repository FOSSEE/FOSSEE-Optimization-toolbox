// Ref : http://www.minlplib.org/st_cqpjk2.html
// var x1 >= 0, <= 1;
// var x2 >= 0, <= 1;
// var x3 >= 0, <= 1;

// minimize obj: 9*x1*x1 - 15*x1 + 9*x2*x2 - 12*x2 + 9*x3*x3 - 9*x3;

// subject to

// e1:    x1 + x2 + x3 <= 10000000000;

H = eye(3,3)*18;
f = [-15;-12;-9];
A = [1 1 1];
b = [10000000000];
lb = [0;0;0];
ub = [1;1;1];

[xopt,fopt] = quadprogCLP(H,f,A,b,[],[],lb,ub)