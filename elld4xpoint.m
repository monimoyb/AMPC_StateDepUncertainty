%% Function to get ellipsoidal approximation of uncertainty for a starting point x_0
%     Monimoy Bujarbaruah. 

function [pd,qd,flag]  = elld4xpoint(x_0,xprev,dprev,Ld)
    % inputs: initial state x_0, L_d and recorded data

    % Data size grows with time. Function complexity increases.

    % output is uncertainty ellipse parameters c and R^-1
    % uncertainty ellipse: (y-pd)^T qd^-1 (y-pd) <= 1

    % Form the SDP to be solved for getting d bound (ACC SysID paper: Appendix)
    
    n = size(xprev,1);
    tau = sdpvar(size(xprev,2),1); 
    qdv  = sdpvar(n,n); pdv = sdpvar(n,1);

    Px = -sum(tau)*eye(n);
    Qx = dprev*tau;
    Rx = -1;

    for j = 1:size(xprev,2)
    Rx = Rx + Ld^2*tau(j)*(x_0'*x_0) ...
                -  2*Ld^2*x_0'*xprev(:,j)*tau(j)...
                + Ld^2*tau(j)*xprev(:,j)'*xprev(:,j)...
                  - tau(j)*dprev(:,j)'*dprev(:,j);
    end


    options3 = sdpsettings('solver','mosek','verbose',0);
    cost = trace((qdv));
    constraints = [tau >= 0; qdv >= 0];

    constraints = [constraints; [Px, Qx, -eye(n);...
                                 Qx', Rx, pdv';...
                                 -eye(n), pdv, -(qdv)]<=0]; 

    exitflag = solvesdp(constraints, cost, options3);
    flag = exitflag.problem; 

    pd = double(pdv); qd = double(qdv); qd = inv(qd);       % uncertainty matrix parameters 

end