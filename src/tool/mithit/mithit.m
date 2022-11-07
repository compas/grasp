%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%											        %%
%%  mithit										        %%
%%                                                                                              %%
%%  Written by: Wenxian Li, March 2019                                                          %%
%%                                                                                              %%
%% This m-file calls subroutines to compute and plot Zeeman splittings of fine and hyperfine    %%
%% levels as a function of magnetic fields, and then compute transition rates between magnetic  %%
%% fine- and hyperfine structure sublevels or the rates of hyperfine induced transitions in the %%
%% field free limit. At last, the program can generate synthetic spectra profiles through       %%
%% convolution with Gaussian function. The mixing between the levels for Bmax are saved in the  %%
%% output file <name>.(c)zm. The transition rates between magnetic sublevels are saved in the   %%
%% output files <name1>.<name2>.tr and <name1>.<name2>.tr.mtrans. The spectra data are written  %%
%% to the output file <name1>.<name2>.s.                                                        %%
%% 						                                                %%
%% This m-file calls for: nuclear.m, plothfszeeman.m, mixingC.m, trans.m and gaussplot.m        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% relcalc        calculations based on a relativistic CI calculation (Y/N)                       %
% I              nuclear spin   								 %
% mu             nuclear magnetic dipole moment						         %
% Q              nuclear electric quadrupole moment 					         %
% UNITB          unit of magnetic field, Tesla (0) or Gauss (1)				         %
% Bmax           upper limit of the magnetic field					         %
% UNITE          energy unit in plot, a.u. (0), cm-1 (1) or MHz (2)			         %
% ITtype         transition type, MIT-fs(0), HIT(1) or MIT-hfs(2)			         %
% JE_i           Breit-Pauli eigenvectors of initial states 				         %
% FE_i           Breit-Pauli-HFS eigenvectors of initial states				         %
% JE_f           Breit-Pauli eigenvectors of final states				         %
% FE_f           Breit-Pauli-HFS eigenvectors of final states				         %
% EM_i           mixing coeficient of initial states					         %
% EM_f           mixing coeficient of final states					         %
% Ttype          transition type, E1/M2(0) M1/E2(1) 					         %
% Parity_i       parity of initial states						         %
% Parity_f       parity of final states							         %
% wl             wavelength vector of synthetic spectra in unit of cm^{-1}		         %
% intens         intensity vector of synthetic spectra					         %
% AM_Mall        transition information between magnetic sublevles selected by user	         %
% AMplots        transition information between magnetic sublevles selected by user to plot      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
format long
warning off;

% Open and read data from file <name>.cgjhfs and <name>.czm
name_i = input("Name of the Initial state:   ",'s');
name_f = input("\nName of the Final state:   ",'s');

relcalc = input("\nAre the calculations based on a relativistic CI calculation? (Y/N) ",'s');

if (relcalc ~= 'y' & relcalc ~= 'Y' & relcalc ~= 'n' & relcalc ~= 'N' )
  relcalc = input('You have to answere (Y/N)','s');
end

ITtype = input ("\nMIT-fs(0), HIT(1) or MIT-hfs(2):  ");

%Input nuclear data
I=0;
mu=0;
Q=0;

if ITtype~=0
  I  = input("\nNuclear spin I: ");
  mu = input("\nNuclear magnetic dipole moment mu:  ");
  Q  = input("\nNuclear electric quadrupole moment Q:  ");
end
%End input nuclear data

%Read unit for magnetic field strength
if ITtype~=1
  % Read unit for magnetic field strength
  UNITB = input("\nB-field in Tesla (0) or Gauss (1):  ");
  Bmax = input("\nGive the upper limit for the B-field:  ");
  N_B = 300;
  Bfield = linspace(0,Bmax,N_B);
  % End read unit for magnetic field strength

  %  Calculate Bau from B.
  if UNITB == 0
    % Tesla.
    Bau = Bfield;
  else
    % Gauss.
    Bau = Bfield/1e4;
  end
  %  end calculate Bau from B.
  UNITE = input("\nEnergies in a.u. (0), cm-1 (1) or MHz (2) ? ");
else
  UNITB=0;
  Bmax=0;
  UNITE=1;
end

N_plots=0;
%Compute energies and mixing coefficients of the magnetic sublevels
display(' ')
display('Start Computation of Energies and Mixing Coefficients')
display('of the Magnetic Sublevels of Initial States')
plothfszeeman(name_i,relcalc,I,mu,Q,UNITB,Bmax,UNITE,ITtype,N_plots);

if ~(strcmp(name_i,name_f))
  display(' ')
  display('Start Computation of Energies and Mixing Coefficients')
  display('of the Magnetic Sublevels of Final States')
  plothfszeeman(name_f,relcalc,I,mu,Q,UNITB,Bmax,UNITE,ITtype,N_plots);
end
%End compute

%Compute transition rates between magnetic sublevels or hyperfine sublevels
Transition = input("\nWould you like to compute the transition rates? (Y/N)  ", 's');
if (Transition == 'y' | Transition == 'Y')

  %Read energies and mixing coefficients of the magnetic sublevels
  [N_eigvec_i,JE_i,FE_i,B_i,unitB_i,EM_i,Parity_i] = mixingC(name_i,relcalc,ITtype,I);
  [N_eigvec_f,JE_f,FE_f,B_f,unitB_f,EM_f,Parity_f] = mixingC(name_f,relcalc,ITtype,I);

  %End read data
  %Define transition type
  if Parity_i==Parity_f
  Ttype = 1;
  else
  Ttype = 0;
  end
  %End define

  [AM_Mall,n]=trans(name_i,name_f,relcalc,ITtype,JE_i,FE_i,JE_f,FE_f,Bmax,EM_i,EM_f,Parity_i,Parity_f,UNITB,I,mu,Q,Ttype);

  %plot spectra
  %synthetic spectra
  Spectra = input("\nWould you like a plot of synthetic spectra? (Y/N)  ", 's');
  if (Spectra == 'y' | Spectra == 'Y')
    MITind_i_s = input("\nGive an index vector of the initial levels(initial level):  ");
    MITind_f_s = input("\nGive an index vector of the final levels(upper level):  ");
    FWHM = input("\nGive the FWHM (in cm-1):  ");

    if ITtype==0
      fileS = strcat(name_i,'.',name_f,'.fs.mit.s'); % MIT-isotopes with nuclear spin I=0
    elseif ITtype==1
      fileS = strcat(name_i,'.',name_f,'.hfs.hit.s'); % HIT-with B=0
    elseif ITtype==2
      fileS = strcat(name_i,'.',name_f,'.hfs.mit.s'); % MIT-isotopes with nuclear spin I~=0
    end

    fpS = fopen(fileS,'w');

    AMplots=[];
    ysum=0;
    Ncol=size(AM_Mall,2);
    for i=1:size(MITind_f_s,2)
      for j=1:size(MITind_i_s,2)
        AMplot=AM_Mall(find(AM_Mall(:,1)==MITind_f_s(i) & AM_Mall(:,7)==MITind_i_s(j)),[Ncol-1:Ncol,3]);
        AMplots=[AMplots;AMplot];
        j=j+1;
      end
      i=i+1;
    end

    wmin=min(AMplots(:,2));
    wmax=max(AMplots(:,2));

    figure(N_plots+1);
    subplot(2,1,1);
    for j=1:size(AMplots,1)
      [wl,intens]= gaussplot(AMplots(j,1),AMplots(j,2),FWHM,wmin,wmax);
      ysum=ysum+intens;
      j=j+1;
  %    pause
    end
    ylabel('Intensity')
    hold off;

    subplot(2,1,2);
    plot(wl,ysum);
    xlabel('Wavelength (cm^{-1})')
    ylabel('Intensity')
    xlim([wmin-FWHM*4,wmax+FWHM*4]);
    hold off;

    for k=1:size(wl,2)
    fprintf(fpS,'%16.9f\t\t%16.9f\n',wl(k),ysum(k));
    end
    fclose(fpS);
  end
  %end synthetic spectra

end
% End Transition

display(' ')
display('MITHIT finished')
display(' ')
