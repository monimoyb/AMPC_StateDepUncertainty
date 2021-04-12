%% Function to get polytopic outer approximation to uncertainty set
%  By: Monimoy Bujarbaruah. 

function [dpolA, dpolb] = polOutdEll(pd,qd)
% input: uncertainty ellipse
% output: its polytopic outer approximation

  [V,D] = eig(qd);                                  % sorted evals and evecs 
  D = D*ones(size(D,1),1);                          % get a column vector     

  % find semi major axis edge points
  unit_vecSemiMaj = V(:,1)/norm(V(:,1),2);
  extrem1 =  pd + 1/sqrt(D(1))*unit_vecSemiMaj;
  extrem2 =  pd - 1/sqrt(D(1))*unit_vecSemiMaj;
  
  % from there build vertices along all other directions
  nsearch_dir = size(V,2)-1; 
  vertices = zeros(4*nsearch_dir,length(D));

  
  for j = 1:nsearch_dir
      unit_vec = V(:,j+1)/norm(V(:,j+1),2);
      vertices(4*(j-1)+1,:) =  extrem1 + 1/sqrt(D(j+1))*unit_vec;
      vertices(4*(j-1)+2,:) =  extrem1 - 1/sqrt(D(j+1))*unit_vec;      
      vertices(4*(j-1)+3,:) =  extrem2 + 1/sqrt(D(j+1))*unit_vec;
      vertices(4*(j-1)+4,:) =  extrem2 - 1/sqrt(D(j+1))*unit_vec;  
  end
  
  pold = Polyhedron('V',vertices);
  dpolA = pold.A; dpolb = pold.b; 

end
