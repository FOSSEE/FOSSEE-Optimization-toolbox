mode(1)
//
// Demo of quadprogCLP.sci
//

//Ref : example 14 :
//https://www.me.utexas.edu/~jensen/ORMM/supplements/methods/nlpmethod/S2_quadratic.pdf
// min. -8*x1*x1 -16*x2*x2 + x1 + 4*x2
// such that
//    x1 + x2 <= 5,
//    x1 <= 3,
//    x1 >= 0,
//    x2 >= 0
H = [2 0
0 8];
f = [-8; -16];
A = [1 1;1 0];
b = [5;3];
lb = [0; 0];
ub = [%inf; %inf];
[xopt,fopt,exitflag,output,lambda] = quadprogCLP(H,f,A,b,[],[],lb,ub)
// Press ENTER to continue
halt()   // Press return to continue
 
//Find x in R^6 such that:
Aeq= [1,-1,1,0,3,1;
-1,0,-3,-4,5,6;
2,5,3,0,1,0];
beq=[1; 2; 3];
A= [0,1,0,1,2,-1;
-1,0,2,1,1,0];
b = [-1; 2.5];
lb=[-1000; -10000; 0; -1000; -1000; -1000];
ub=[10000; 100; 1.5; 100; 100; 1000];
//and minimize 0.5*x'*H*x + f'*x with
f=[1; 2; 3; 4; 5; 6]; H=eye(6,6);
[xopt,fopt,exitflag,output,lambda]=quadprogCLP(H,f,A,b,Aeq,beq,lb,ub)
halt()   // Press return to continue
 
//Solving Linear Programming Problem
// min -x0 - x1
// subject to
//    x0+2x1 <= 3
//  2x0+x1 <= 3
f = [-1;-1];
A = [1,2;2,1];
b = [3;3];
[xopt,fopt] = quadprogCLP([],f,A,b);
//========= E N D === O F === D E M O =========//
