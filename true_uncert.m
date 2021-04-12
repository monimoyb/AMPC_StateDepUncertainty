%% Uncertainty function for true system (L_d Lipschitz)
%  Monimoy Bujarbaruah

function val = true_uncert(x,Ld)
%%% d(x) value at any given x
      val = Ld*[atan(x(1));x(2)]; 
end
