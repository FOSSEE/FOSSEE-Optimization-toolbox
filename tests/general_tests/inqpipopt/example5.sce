//intcon out of range

H = [1 -1; -1 2]; 
f = [-2; -6];
intcon = [200];
A = [1 1; -1 2; 2 1];
b = [2; 2; 3];
lb=[0,0];
ub=[%inf, %inf];

path = get_absolute_file_path('example5.sce');

[xopt,fopt,status,output]=intqpipopt(H,f,intcon,A,b,[],[],lb,ub,[0 0],list(),[path + '/bonmin.opt'])

clear path;

// "bonmin.opt" contains bonmin file options
// Current options enables printing output to the terminal which can be used for debugging and other purposes.
// For more details visit their website http://coin-or.org/Bonmin