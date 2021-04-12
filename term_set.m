%% Calculating the Terminal RPI Set X_N for Robust AMPC
%  Monimoy Bujarbaruah. 

function [Xn,Pinf] = term_set(A,B,C,D,b,Q,R,U,W,nx,nu)

  [Finf,Pinf] = dlqr(A,B,Q,R);
  Acl = A-B*Finf;                   % Closed loop system for terminal set 

  model = LTISystem('A',Acl); 
  S = Polyhedron('A',C-D*Finf,'b',b); 
  S = S.intersect(Polyhedron('H',[-U.H(:,1:nu)*Finf U.H(:,nu+1)])); 

  Xn = MRPI(model, S, Polyhedron.emptySet(nx), W);

end

