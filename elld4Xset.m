%% Function to get ellipsoidal approximation of uncertainty for a set of X_t
%  By: Monimoy Bujarbaruah. 

function [pd,qd,flag]  = elld4Xset(px,qx,xprev,dprev,Ld)
    % inputs are parameters of the X_t ellipse, L_d and recorded data

    % X_t ellipse is: (x-px)^Tqx(x-px) <= 1
    % Data size grows with time. Function complexity increases.

    % output is uncertainty ellipse parameters pd and qd
    % uncertainty ellipse: (y-pd)^Tqd(y-pd) <= 1

    % Form the SDP to be solved for getting d bound (Appendix of the paper)
    n = size(xprev,1);
    rho = sdpvar(1,1); tau = sdpvar(size(xprev,2),1); 
    qdv  = sdpvar(n,n); pdv = sdpvar(n,1);

    Px = -rho*qx + Ld^2*sum(tau)*eye(n);
    Qx = rho*qx'*px - Ld^2*xprev*tau;
    Rx = -sum(tau)*eye(n);
    Sx = dprev*tau;
    Tx = -1+ rho*(1-px'*qx*px);

    for j = 1:size(xprev,2)
    Tx = Tx + Ld^2*tau(j)*xprev(:,j)'*xprev(:,j) ...
            - tau(j)*dprev(:,j)'*dprev(:,j);
    end


    options2 = sdpsettings('solver','mosek','verbose',0);
    cost = trace((qdv));
    constraints = [rho>=0; tau >= 0; qdv >=0];

    constraints = [constraints; [Px, zeros(n,n), Qx, zeros(n,n);...
                                 zeros(n,n), Rx, Sx, -eye(n);...
                                 Qx', Sx', Tx, pdv';...
                                 zeros(n,n), -eye(n), pdv, -(qdv)]<=0]; 

    exitflag = solvesdp(constraints, cost, options2);
    flag = exitflag.problem; 

    pd = double(pdv); qd = double(qdv); qd = inv(qd);       % uncertainty matrix parameters 

end