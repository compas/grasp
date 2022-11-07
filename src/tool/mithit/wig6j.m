%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                %

function w6j = wig6j(j1,j2,j3,l1,l2,l3)

% Computes the Wigner 6j-symbol {j1 j2 j3}                                                       %
%                               {l1 l2 l3}                                                       %
%                                                                                                %
% Expression from R. Cowan, "The Theory of Atomic Structure and Spectra",                        %
% formula 5.23.                                                                                  %
%                                                                                                %
% Written by Per Jonsson, july 2005                                                              %
%                                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if j1 >= 0 & j2 >= 0 & j3 >= 0 & l1 >= 0 & l2 >= 0 & l3 >= 0 & ... % quantities >= 0
   mod(4*j1,2) == 0 & mod(4*j2,2) == 0 & mod(4*j3,2) == 0 & ...    % ji,li integral
   mod(4*l1,2) == 0 & mod(4*l2,2) == 0 & mod(4*l3,2) == 0 & ...    % or half integral
   mod(2*(j1+j2+j3),2) == 0 & ...                                  % j1+j2+j3 integral
   j1 + j2 >= j3 & j2 + j3 >= j1 & j3 + j1 >= j2 & ...             % triangle relations
   mod(2*(j1+l2+l3),2) == 0 & ...                                  % j1+l2+l3 integral
   j1 + l2 >= l3 & l2 + l3 >= j1 & l3 + j1 >= l2 & ...             % triangle relations
   mod(2*(l1+j2+l3),2) == 0 & ...                                  % l1+j2+l3 integral
   l1 + j2 >= l3 & j2 + l3 >= l1 & l3 + l1 >= j2 & ...             % triangle relations
   mod(2*(l1+l2+j3),2) == 0 & ...                                  % l1+l2+j3 integral
   l1 + l2 >= j3 & l2 + j3 >= l1 & j3 + l1 >= l2                   % triangle relations

   a = j1; b = j2; c = j3;
   d1 = sqrt(factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/factorial(a+b+c+1));
   a = j1; b = l2; c = l3;
   d2 = sqrt(factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/factorial(a+b+c+1));
   a = l1; b = j2; c = l3;
   d3 = sqrt(factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/factorial(a+b+c+1));
   a = l1; b = l2; c = j3;
   d4 = sqrt(factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/factorial(a+b+c+1));

   kmin = max([j1+j2+j3 j1+l2+l3 l1+j2+l3 l1+l2+j3]);
   kmax = min([j1+j2+l1+l2 j2+j3+l2+l3 j3+j1+l3+l1]);

   w6j = 0;
   for k = kmin:kmax
     numerator = factorial(k+1);
     denominator = factorial(k-j1-j2-j3)*factorial(k-j1-l2-l3)*factorial(k-l1-j2-l3)* ...
         factorial(k-l1-l2-j3)*factorial(j1+j2+l1+l2-k)*factorial(j2+j3+l2+l3-k)*factorial(j3+j1+l3+l1-k);
     w6j = w6j + (-1)^k*numerator/denominator;
   end
   w6j = w6j*d1*d2*d3*d4;
else
   w6j = 0;
end

