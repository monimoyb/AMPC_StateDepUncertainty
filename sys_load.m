%% Defining System and Constraint Matrices with this function 
%  Monimoy Bujarbaruah

function [A,B,C,D,b,X,U,nx,nu,Ld,x_0,Xlb,Xub,Ulb,Uub] = sys_load()
 
%%%% Considering two states 
A = [1.2, 1.3; 0, 1.5]; 
B = [0;1];   
nx = size(A,2); 
nu = size(B,2); 

%%%% Considering constraints of the form -a<=x(i)<=a and -ulb<=u<=uub
%          and expressing in Cx+Du <=b format 
Xlb =  -[1; 1];  
Xub =   [1.4; 3];
Ulb =  -4;
Uub =   1; 
C = [1 0; 0 1; -1 0; 0 -1; 0 0; 0 0]; D = [0; 0; 0; 0; 1; -1]; b = [Xub;-Xlb;Uub;-Ulb]; 
X = Polyhedron('lb',Xlb,'ub',Xub);
U = Polyhedron('lb',Ulb,'ub',Uub);   
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Uncertainty Lipschitz condition
Ld = .05;                                   % Lipschitz constant 

%% Initial condition for starting control process
x_0 = [-1;2]; 

end