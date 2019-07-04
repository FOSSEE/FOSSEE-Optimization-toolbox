// Ref : http://www.minlplib.org/st_e23.html

// var x1 >= 0, <= 5;
// var x2 >= 0, <= 5;

// minimize obj: x1*x2 - x1 - x2;

// subject to

// e1:  - 6*x1 + 8*x2 <= 3;

// e2:    3*x1 - x2 <= 3;

H = [0 1;
     1 0];
f = [-1;-1];
A = [-6 8;
     3 -2];
b = [3;3];
lb = [0;0];
ub = [5;5];

[xopt,fopt] = quadprogCLP(H,f,A,b,[],[],lb,ub)