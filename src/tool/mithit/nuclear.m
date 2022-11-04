%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                %

function M = nuclear(mu,Q,I)

%                                                                                                %
% This function computes the reduced nuclear matrix elements                                     %
% in atomic units.                                                                               %
%                                                                                                %
% Written by Per Jonsson and Martin Andersson, July 2006                                         %
% Modified by Wenxian Li, March 2019                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

me = 5.48579909070D-04;      % electron mass in amu
a0cm = 5.2917721067e-9;      % Bohr radius in cm
mp = 1.007276466879e0 ;      % proton mass in amu
cvac = 1.37035999139e2;
mu_au = mu/2*me/(mp*cvac);
Q_au  = Q*1e-24/a0cm^2;
if (I > 0)
  M(1) = mu_au*sqrt((I+1)*I)/I;
  if (I > 0.5)
    M(2) = Q_au/2*sqrt((2*I+3)*(I+1))/(I*(2*I-1));
  else
    M(2) = 0;
  end
else
  M(1) = 0;
  M(2) = 0;
end

M(1) = M(1);
M(2) = M(2);
