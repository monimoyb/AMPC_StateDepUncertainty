%% Function to get robust successor states 
%  By: Monimoy Bujarbaruah. 
%  Following Borrelli MPC book notations: Section 10.3.2

function [HX_RRS, hX_RRS] = robSucc(A,B,UA,Ub,XA,Xb,PdA,Pdb,XfullA, Xfullb)
    % inputs: system matrices. Input range. Uncertainty outer polytope.
    % Starting states range (X)
    % outputs: robust successor states

    W =  Polyhedron('A',PdA,'b',Pdb);
    X  = Polyhedron('A',XA,'b',Xb); 
    U  = Polyhedron('A',UA,'b',Ub); 

    model = LTISystem('A',A,'B',B); 
    RobReach = model.reachableSet('X', X, 'U', U, 'direction', 'forward') + W;

    Xfull = Polyhedron('A',XfullA,'b',Xfullb);
    RobReach = intersect(RobReach,Xfull);

    HX_RRS = RobReach.A;
    hX_RRS = RobReach.b;
    
end
