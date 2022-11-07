%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                %

function w3j = wig3j(j1,j2,j3,m1,m2,m3)

% Computes the Wigner 3j-symbol (j1 j2 j3)                                                       %
%                               (m1 m2 m3)                                                       %
%                                                                                                %
% Expression from R. Cowan, "The Theory of Atomic Structure and Spectra",                        %
% formula 5.1.                                                                                   %
%                                                                                                %
% Written by Per Jonsson, july 2005                                                              %
%                                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if j1 >= abs(m1) & j2 >= abs(m2) & j3 >= abs(m3) & ...             % ji >= abs(mi) >= 0
   mod(4*j1,2) == 0 & mod(4*j2,2) == 0 & mod(4*j3,2) == 0 & ...    % ji integral or half integral
   mod(2*(j1-m1),2) == 0 & mod(2*(j2-m2),2) == 0 & ...             % mi integral when ji integral
   mod(2*(j3-m3),2) == 0 & mod(2*(j1+m1),2) == 0 & ...             % and mi half integral when ji
   mod(2*(j2+m2),2) == 0 & mod(2*(j3+m3),2) == 0 & ...             % half integral
   mod(2*(m1+m2+m3),2) == 0 & mod(2*(j1-j2-m3),2) == 0 & ...       % m1+m2+m3 and j1-j2-m3 integral
   mod(2*(j1+j2+j3),2) == 0 & ...                                  % j1+j2+j3 integral
   j1+j2 >= j3 & j2+j3 >= j1 & j3+j1 >= j2 & ...                   % triangle relations
   m1+m2+m3 == 0

   kmin = max([0 j2-j3-m1 j1-j3+m2]);
   kmax = min([j1+j2-j3 j1-m1 j2+m2]);

   w3j = 0;
   for k = kmin:kmax
     numerator = 1;
     denominator = factorial(k)*factorial(j1+j2-j3-k)*factorial(j1-m1-k)* ...
         factorial(j2+m2-k)*factorial(j3-j2+m1+k)*factorial(j3-j1-m2+k);
     w3j = w3j + (-1)^k*numerator/denominator;
   end
   numerator = factorial(j1+j2-j3)*factorial(j1-j2+j3)*factorial(-j1+j2+j3)* ...
       factorial(j1-m1)*factorial(j1+m1)*factorial(j2-m2)* ...
       factorial(j2+m2)*factorial(j3-m3)*factorial(j3+m3);
   denominator = factorial(j1+j2+j3+1);
   w3j = w3j*(-1)^(j1-j2-m3)*sqrt(numerator/denominator);
else
   w3j = 0;
end


