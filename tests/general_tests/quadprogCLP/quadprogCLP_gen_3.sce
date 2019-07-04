// Ref : http://www.minlplib.org/st_cqpjk1.html
// var x1 >= 0;
// var x2;
// var x3 >= -10000, <= 10000;
// var x4 >= -10000, <= 10000;

// minimize obj: 2*x1*x1 - 1.33333*x1 + 4*x2*x2 - 2.66667*x2 + 6*x3*x3 - 4*x3 + 
//     0.5*x4*x4 - 10*x4;

// subject to

// e1:  - x1 - x2 - x3 - x4 <= -1;

// e2:    x1 + x2 + x3 + x4 <= 1;


H = [4 0 0 0;
     0 8 0 0;
     0 0 12 0;
     0 0 0 1];
f = [-1.33333;-2.66667;-4;-10];
A = [-1 -1 -1 -1;
      1  1  1  1];
b = [-1 ; 1];
lb = ones(4, 1)*-10000;
ub = ones(4, 1)*10000;

[xopt,fopt] = quadprogCLP(H,f,A,b,[],[],lb,ub)