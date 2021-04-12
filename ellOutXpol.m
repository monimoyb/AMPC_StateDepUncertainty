%% Function to get ellipsoidal outer-approximation to a Polytope
% Follow El Ghaoui, Calafiore Optimization Models book Section 11.4.2.2
% By: Monimoy Bujarbaruah.

function [px, qx, flag] = ellOutXpol(H,h)
    % Outputs are center and shape of the outer ellipse
    % Ellipse is: (x-px)^Tqx(x-px) <=1
    
    % Input polytope representation is Hx <= h
    % Convert to vertex notations and obtain vertices
    
    n = size(H,2);
    X = Polyhedron('A',H,'b',h);
    Vert = X.V;                              % extracted vertices

    % Form the SDP to be solved
    Av = sdpvar(n,n);
    bv = sdpvar(n,1);                      % decision variables

    cost = -log(det(Av));
    constraints = Av>=0;

    for j = 1:size(Vert,1)
         constraints = [constraints, norm(Av*Vert(j,:)' - bv,2)^2 <=1];
    end

    options1 = sdpsettings('solver','mosek','verbose',0);

    exitflag = solvesdp(constraints,cost,options1);
    flag = exitflag.problem;                                                % solver flag

    A = double(Av); b = double(bv);

    qx = A^2; px = A\b;                                                   % Book chapter 11.4.2.2

end
