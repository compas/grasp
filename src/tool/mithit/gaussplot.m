%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                %
%                                                                                                %

function [wl,intens] = gaussplot(A,xc,FWHM,wmin,wmax)

%                                                                                                %
%This m-file is the Gauss function called to model Gauss profile spectra.                        %
%                                                                                                %
% Written by Wenxian Li, March 2019                                                              %
%                                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nw=10000;             %number of grid points in the spectra plots
dd=FWHM/(2*sqrt(2*log(2)));
wl = linspace(wmin-FWHM*4,wmax+FWHM*4,Nw);  % wavelength range
intens = A/(sqrt(2*pi)*dd)*exp(-(wl-xc).^2/(2*dd^2));
plot(wl,intens);
xlim([wmin-FWHM*4,wmax+FWHM*4])
hold on;

end
% end function
