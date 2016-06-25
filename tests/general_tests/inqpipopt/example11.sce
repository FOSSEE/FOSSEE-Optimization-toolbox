//Find x in R^6 such that:
// An example with equality, inequality constraints, variable bounds and intial guess

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
x0 = repmat(0,6,1);
intcon=[1];
options=list('MaxNodes',0);

path = get_absolute_file_path('example11.sce');

[xopt,fopt,status,output]=intqpipopt(H,f,intcon,A,b,Aeq,beq,lb,ub,x0,options,path+'bonmin.opt')
