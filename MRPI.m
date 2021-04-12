%% Maximal Disturbance Invariant Set Computation 
% Monimoy Bujarbaruah

function P = MRPI(model, X, U,  W)
% computation of a robust invariant set for given LTImodel 
%
maxIterations = 100;
X0 = X; % initial set constraint

for j = 1:maxIterations
    % subtract noise    
    S = X0 - W;
    % backward reachable set
    R = model.reachableSet('X', S, 'U', U,...
        'direction', 'backward');
    % intersect with the state constraints
    P = R.intersect(X0);

    if P==X
        break
    else
        X0 = P;
    end
end
