%% Robust Adaptive MPC with State-Dependent Uncertainty
%     This will print out the results of the paper https://ieeexplore.ieee.org/document/9143777
%     By: Monimoy Bujarbaruah
%% Flow: Get uncertainty bounds from any x_t
%  elld4xpoint -> polOutEll -> robSucc (use for next rob succ) -> ellOutXpol -> elld4Xset -> polOutEll
%  Repeat this along horizon 

%% Start the main code here 
clear all
close all
clc
warning off
yalmip 'clear'
%% All system parameters
adap_flag = 1;                                                                   % should adapt or not?
N = 3;                                                                               % MPC horizon length.
[A,B,C,D,b,X,U,nx,nu,Ld,x_0,Xlb,Xub,Ulb,Uub] = sys_load(); 
UA = U.A; Ub = U.b; 

%%% Matrices for cost function and terminal set 
Q =  10*eye(nx);
R =   2*eye(nu);

%%% Start condition here. This on, will go into data collected
x_start = [1;1];                                     % Picking a starting condition
d_start = true_uncert(x_start,Ld);           % Getting Ld value ;

%%% Arrays will be built from here on 
x_prev = x_start;                                  % start from x_start onward 
d_prev = d_start; 

%% Keep gathering data with rand inputs to get a non-empty terminal set 
term_emptyFlag = 1;                         % Assume terminal set is empty!
x_nxt = x_prev; 
count = 0; 
[pX, qX, flX]  =   ellOutXpol(X.A,X.b); 

while term_emptyFlag == 1
 uP = rand(nu,1);
 x_nxt = A*x_nxt + B*uP + true_uncert(x_nxt,Ld); 
 x_prev  = [x_prev, x_nxt]; 
 d_prev  = [d_prev, true_uncert(x_nxt,Ld)];    

%%% Calculating Terminal Set and terminal cost weight 
[pdX,qdX,flag] =    elld4Xset(pX,qX,x_prev(:,1:end),d_prev(:,1:end), Ld);

 if flag ~= 1

    [WTermA, WTermb]  = polOutdEll(pdX,qdX); 
    WTerm = Polyhedron('A',WTermA,'b',WTermb); 

    [Xn,Pinf] = term_set(A,B,C,D,b,Q,R,U,WTerm,nx,nu); 
 
    if isEmptySet(Xn) == 0
        term_emptyFlag = 0; 
    end
    
 end
 
 count = count+1;
 yalmip 'clear'

end

Xn0 = Xn; 
dim_t = size(C,1)*N + size(Xn.A,1);                           
%% Forming matrices appearing in the optimization problem 
[capA, capE, capB, capC, capD, Aw_batch, Bu_batch, A_batch] = obtain_mat(A,B,C,D,Xn,nx,nu,N,dim_t);
matF = capC*capB + capD; 
matG = capC*capE; 
matH = -capC*capA; 
mat_c = [kron(ones(N,1),b); Xn.b];        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% START CONTROL PROCESS FROM HERE ON 
simsteps = 10;                                                          % Number of steps in cl-loop
options = sdpsettings('solver','gurobi','verbose',0);
x_cl = zeros(nx,simsteps+1);
u_cl = zeros(nu,simsteps);  
x_cl(:,1) = x_0;
polprev = [];
for j = 1:N
    polprev = [polprev, Polyhedron('lb',-inf*ones(nx,1),'ub',inf*ones(nx,1))]; 
end

figure;

for i=1:simsteps

    if (mod(i,5)==0)                                              % Recompute only after certain intervals
       
       %% Calculating Terminal Set and terminal cost weight 
       [pdX,qdX,~] =   elld4Xset(pX,qX,x_prev(:,1:end),d_prev(:,1:end), Ld);
       [WTermA, WTermb]  = polOutdEll(pdX,qdX); 

       %%% Including intersection of previous polytope too
       WTerm = intersect(WTerm, Polyhedron('A',WTermA,'b',WTermb));

       [Xn,Pinf] = term_set(A,B,C,D,b,Q,R,U,WTerm,nx,nu); 
       dim_t = size(C,1)*N + size(Xn.A,1);                           

        %% Forming matrices appearing in the optimization problem 
        [capA, capE, capB, capC, capD, Aw_batch, Bu_batch, A_batch] = obtain_mat(A,B,C,D,Xn,nx,nu,N,dim_t);

        matF = capC*capB + capD; 
        matG = capC*capE; 
        matH = -capC*capA; 
        mat_c = [kron(ones(N,1),b); Xn.b]; 

    end
  
    constraints = []; 
  
    %% Form All the Successor States and Uncertainties here
    [pd0,qd0,~] = elld4xpoint(x_cl(:,i),x_prev(:,1:end),d_prev(:,1:end),Ld); 

    Xt = Polyhedron('ub',x_cl(:,i),'lb',x_cl(:,i));
    XtA = Xt.A; Xtb = Xt.b;  
    Hs=[]; hs =[]; 

    for k = 1:N-1
      % polytope bound of uncertainty 
      [d0A, d0b]  = polOutdEll(pd0,qd0); 
      polD = Polyhedron('A',d0A,'b',d0b);
      polIn = intersect(polD,polprev(k)); 

      d0A = polIn.A; d0b = polIn.b; 

      Hs = blkdiag(Hs, d0A); hs = [hs; d0b]; 
      % computing succesor states
      [XtA, Xtb] = robSucc(A,B,UA,Ub,XtA,Xtb,d0A,d0b,X.A,X.b); 
      % outer ellipse to successor states
      [px, qx, ~] = ellOutXpol(XtA,Xtb); 
      % obtaining uncertainty with successor outer ellipse 
      [pd0,qd0,~] = elld4Xset(px,qx,x_prev(:,1:end),d_prev(:,1:end),Ld);

      if k ~=1
        polprev(k-1) = polD; 
      end

    end

    % polytope bound of uncertainty 
    [d0A, d0b]  = polOutdEll(pd0,qd0); 
    polD = Polyhedron('A',d0A,'b',d0b);
    polIn = intersect(polD,polprev(N)); 

    d0A = polIn.A; d0b = polIn.b; 
    Hs = blkdiag(Hs, d0A); hs = [hs; d0b]; 

    polprev(N-1) = polD; 
    polprev(N) = WTerm; 

    dim_a = size(Hs,1);   % 'a' dimension from https://www.sciencedirect.com/science/article/pii/S0005109806000021 notations  
   
   %% Creating Open Loop Optimization Variables for MPC 
    sl = sdpvar(size(matF,1),1);  % the slacks are small (~10^-6). This is for numerical reasons
    Wsl = 10000;                      % high penalty on nonzero slacks

    % M and v are input decision variables
    M = sdpvar(nu*N,nx*N); 
    for j=1:nu*N
        for k = 2*j-1:nx*N
                M(j,k) =0;
        end
    end
    v = sdpvar(nu*N,1); 
    % dual variables
    Z = sdpvar(dim_a,dim_t);     
  
    %% Solving the Optimization Problem 
    x_pred = x_cl(:,i); 
    cost_state = x_pred'*Q*x_pred; 
    
    for k=1:N
        x_pred  = A*x_pred + B*v(1+(k-1)*nu:k*nu,1);
        if k~=N
            cost_state = cost_state + x_pred'*Q*x_pred;
        end
    end
    cost_state = cost_state + x_pred'*Pinf*x_pred + sl'*Wsl*sl;                
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    constraints = [constraints; matF*v + Z'*hs <= mat_c + matH*x_cl(:,i) + sl];
    constraints = [constraints; Z>=0; sl>=0];
    constraints = [constraints, matF*M + matG == Z'*Hs];

    obj_ol = v'*kron(R,N)*v + cost_state; 
    diagn=solvesdp(constraints, obj_ol, options);
    
    %% Obtaining Closed Loop Parameters 
    v_hor = double(v);
    u_cl(:,i) = v_hor(1:nu,:);                % Closed loop control 

    %% Actual System Simulation in Closed Loop 
    w = true_uncert(x_cl(:,i),Ld);                      
    x_cl(:,i+1) = A*x_cl(:,i) + B*u_cl(:,i) +  w;  

  %% Apprending to measured pairs for envelope improvement 
    if adap_flag == 1
        x_prev = [x_prev, x_cl(:,i)];
        d_prev = [d_prev, w];
    end
    
  %%% Pick a specific point and show how adding data is helping (2D ONLY)
    x_probe = [[-1;2],[1;1], [-1;1], [-2;-1]]; 
    
    for kk = 1:size(x_probe,2)
        [pd_probe,qd_probe,flg_probe] = elld4xpoint(x_probe(:,kk),x_prev(:,1:end),d_prev(:,1:end), Ld); 
        [evec_prb, eval_prb] = eig(qd_probe); 
 
        ang_probe = dot(evec_prb(:,1),[1;0])/(norm(evec_prb(:,1),2)*norm([1;0],2)); 
        ang_probe = acos(ang_probe);
        subplot(2,2,kk)
        ellipse(1/sqrt(min(eig(qd_probe))),1/sqrt(max(eig(qd_probe))),ang_probe...
                    ,pd_probe(1),pd_probe(2)); 
        uncr = true_uncert(x_probe(:,kk),Ld);
        hold on; plot(uncr(1),uncr(2),'*','MarkerSize',10,'Color','k');
        
        grid on; hold on;
    end
    
    % some error displays
    feasmpc = diagn.problem
    iter_count = i
    flg_probe
    
end

%% Plotting
figure; 
plot(Xn0, 'color','k','alpha',0,'linestyle',':','linewidth',2); 
hold on; plot(Xn,'color','r','alpha',0);
legend({'Nonempty $\mathcal{X}_N$ at $t=0$',...
    'Improved $\mathcal{X}_N$ at $t>0$'},'Interpreter','latex','fontsize',20);

%% Inputs
figure; 
plot(u_cl,'linewidth',2); grid on; 
hold on; plot(Uub(1)*ones(1,simsteps+1),'linewidth',2,'color','k');
hold on; plot(Ulb(1)*ones(1,simsteps+1),'linewidth',2,'color','k');
legend('u','Constraints');

%% States
figure; 
subplot(2,1,1)
plot(x_cl(1,:),'linewidth',2); grid on;
hold on; plot(Xub(1)*ones(1,simsteps+1),'linewidth',2,'color','k');
hold on; plot(Xlb(1)*ones(1,simsteps+1),'linewidth',2,'color','k');
legend('x_1','Constraints');

subplot(2,1,2)
plot(x_cl(2,:),'linewidth',2); grid on;
hold on; plot(Xub(2)*ones(1,simsteps+1),'linewidth',2,'color','k');
hold on; plot(Xlb(2)*ones(1,simsteps+1),'linewidth',2,'color','k');
legend('x_2','Constraints');