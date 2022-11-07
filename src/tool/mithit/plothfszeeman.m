%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                          %%

function plothfszeeman(name,relcalc,I,mu,Q,UNITB,Bmax,UNITE,ITtype,N_plots)

%%  Written by: Martin Andersson and Per J\"onsson, july 2006               %%
%%  Modified by: Wenxian Li, March 2018                                     %%
%%                                                                          %%
%%                                                                          %%
%%  This m-file opens the file <name>.gjhfs or <name>.cgjhfs and reads the  %%
%%  matrix of the electronic part of the Zeeman  interaction, the magnetic  %%
%%  dipole  hyperfine interaction  and the electric  quadropole  hyperfine  %%
%%  interaction.  The interaction matrix for a number  of magnetic fields,  %%
%%  B, are constructed. For each B value the magnetic sublevels are obtai-  %%
%%  ned by diagonalizing the interaction matrix. The magnetic sublevels of  %%
%%  user defined  combinations of  Breit-Pauli  fine-structure  levels, of  %%
%%  Breit-Pauli  hyperfine structure levels are then  plotted as functions  %%
%%  of B in the interval [0, Bmax]. The mixing between the levels for Bmax  %%
%%  are saved in the output file <name>.zm or <name>.czm                    %%
%%                                                                          %%
%%  This file calls for: nuclear.m, wig3j.m and wig6j.m                     %%
%%                                                                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                            %
%  B           B-field vector in unit (Tesla or Gauss) defined by user.      %
%  Bau         B-field in a.u.                                               %
%  Bmax        Upper limit of the magnetic field                             %
%  BPLEV       vector keeping trak of which hyperfine levels belongs to      %
%                which Breit-Pauli level                                     %
%  EBP         vector of energies for the Breit-Pauli eigenvectors.          %
%  EHFS        vector of energies for the Breit-Pauli-HFS eigenvectors.      %
%  Fmax        max F-value. If no nuclear spin max J-value                   %
%  Fmin        min F-value. If no nuclear spin min J-value                   %
%  Fval        vector of F values for the Breit-Pauli-HFS eigenvectors. If   %
%                no nuclear spin, same as J                                  %
%  HDCB(:,:)   hartree-dirac-coulomb-breit-pauli energy matrix               %
%  hfslev      vector used to print mixing coeficients in the right          %
%                possition in the output file                                %
%  Hhfs        hyperfine interaction matrix                                  %
%  HM          Magnetic part of the interaction matrix.                      %
%  Htot        Total interaction matrix.                                     %
%  I           nuclear spin                                                  %
%  J           vector of J values for the Breit-Pauli eigenvectors.          %
%  Lev         vector keeping trak of which Breit Pauli level for a certain  %
%                J-value                                                     %
%  M(1:2)      reduced nuclear matrix elements in a.u., computed in          %
%                nuclear.m                                                   %
%  MF          MF quantum number. If no nuclear spin same as MJ              %
%  mix         vector used to print the mixing coeficient to the ouput file  %
%  mu          nuclear magnetic dipole moment                                %
%  N_B         number of grid points in the plots                            %
%  N_eigvec    number of Breit-pauli eigenvectors.                           %
%  N_Feigvec   number of Breit-pauli-HFS eigenvectors.                       %
%  N_plots     number of plots which will be created.                        %
%  NT(:,:)     electronic part of the zeeman interaction matrix              %
%  Parity      vector of parity of the Breit-Pauli eigenvectors.             %
%  Q           nuclear electric quadrupole moment                            %
%  T(:,:,1)    electronic part of the magnetic dipole hyperfine interaction  %
%                matrix                                                      %
%  T(:,:,2)    electronic part of the electric quadropole hyperfine          %
%                interaction matrix                                          %
%  UNITB       unit of magnetic field. 0=Tesla and 1=Gauss                   %
%  UNITE       energy unit in plot. 0=a.u., 1=cm-1 and 2=MHz                 %
%                                                                            %
%                                                                            %
%  VARIABLEout, VARIABLEsort, VARIABLEindex, VARIABLEdummy all used to sort  %
%  variables for different purposes                                          %
%                                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

format long
warning('off', 'all');

% Open and read data from file <name>.gjhfs or <name>.cgjhfs
if (relcalc == 'y' | relcalc == 'Y')
  file = strcat(name,'.cgjhfs');
elseif (relcalc == 'n' | relcalc == 'N')
  file = strcat(name,'.gjhfs');
end

fp = fopen(file,'r');
[string1,c] = fscanf(fp,'%s',4);
[N_eigvec,c] = fscanf(fp,'%u',1);
[string2,c] = fscanf(fp,'%s',4);

% Read J values Parity and Energies

Lev = [];
J = [];
Parity = [];
EBP = [];
for i =1:N_eigvec
  Lev = [Lev fscanf(fp,'%f',1)];
  J = [J fscanf(fp,'%f',1)];
  Parity = [Parity fscanf(fp,'%s',1)];
  EBP = [EBP fscanf(fp,'%f',1)];
end

% Checking is all levels have the same parity

for i = 1:N_eigvec
  if (Parity(1) ~= Parity(i))
    error('All levels have to have the same parity')
  end
end

% Read the Zeeman and hyperfine interaction matrices.
% Matlab puts each line as a column so the matrices are
% mirrored in the diagonal.

[string3,c] = fscanf(fp,'%s',3);
[NT,c] = fscanf(fp,'%e',[N_eigvec,N_eigvec]);
[string3,c] = fscanf(fp,'%s',6);
[T1,c] = fscanf(fp,'%e',[N_eigvec,N_eigvec]);
[string3,c] = fscanf(fp,'%s',6);
[T2,c] = fscanf(fp,'%e',[N_eigvec,N_eigvec]);
T = zeros(N_eigvec,N_eigvec,2);
T(:,:,1) = T1;
T(:,:,2) = T2;

fclose(fp);
if (relcalc == 'y' | relcalc == 'Y')
  file = strcat(name,'.czm');
else
  file = strcat(name,'.zm');
end
fp = fopen(file,'w');

% compute reduced nuclear matrix elements
M = nuclear(mu,Q,I);

% determine the range of allowed F values
Fmax = max(I+J);       % J is a vector
Fmin = min(abs(I-J));  %

% determine vector with F values

Fval = [];
for F_ = Fmax:-1:Fmin
  % for current F determine the range of allowed J values
  Jmin = abs(F_-I);
  Jmax = F_+I;
  % loop over the J dependent states and check if they
  % can be coupled to the current F
  for i = 1:N_eigvec
    if Jmin <= J(i) & J(i) <= Jmax
      Fval = [Fval F_];
    end
  end
end

N_Feigvec = length(Fval);

% Read unit for magnetic field strength
N_B = 300;
B = linspace(0,Bmax,N_B);

%  Calculate Bau from B.
if UNITB == 0
   % Tesla to a.u.
   Bau = B/(2.35051808e5);
else
   % Gauss to a.u.
   Bau = B/(2.35051808e9);
end

%  End calculate Bau from B.

% Read energy unit in plots

fprintf(fp,'%s','  B  = ');
fprintf(fp,'%15.7f', Bmax);
if UNITB == 0
   fprintf(fp,'%s',' Tesla');
else
   fprintf(fp,'%s',' Gauss');
end
fprintf(fp,'%s','   N_EIGVEC = ');
fprintf(fp,'%4.0f',N_Feigvec);
fprintf(fp,'\n');
fprintf(fp,'\n');


EHFS = [];
BPLEV = [];
for F = Fmax:-1:Fmin
  clear H;
  % BPLEVindex and EHFSdiag are used for creating BPLEV which keeps trak on which
  % Breit-Pauli level belongs to which hyperfine eigenvalue
  BPLEVindex = [];
  EHFSdiag = [];
  Jmin = abs(F-I);
  Jmax = F+I;
  nrow = 0;
  for i = 1:N_eigvec
    ncol = 0;
    if (Jmin <= J(i) & J(i) <= Jmax)
      BPLEVindex = [BPLEVindex i];
      nrow = nrow + 1;
      for j = 1:i
        if (Jmin <= J(j) & J(j) <= Jmax)
          ncol = ncol + 1;
          % add hyperfine interaction contributions
          s = I+J(i)+F;
          w6j = wig6j(I,J(j),F,J(i),I,1);
          Hdipol = (-1)^s*w6j*sqrt((2*J(j)+1)*(2*I+1))* T(i,j,1)*M(1);
          w6j = wig6j(I,J(j),F,J(i),I,2);
          Hquadro = (-1)^s*w6j*sqrt((2*J(j)+1)*(2*I+1))*T(i,j,2)*M(2);
          H(nrow,ncol) = Hdipol + Hquadro;
          if i == j
            H(nrow,ncol) = H(nrow,ncol)  + EBP(j);
            EHFSdiag = [EHFSdiag H(nrow,ncol)];
          end
          H(ncol,nrow) = H(nrow,ncol);
        end
      end
    end
  end

  [EHFSdiag sort3] = sort(EHFSdiag);
  BPLEVindex = BPLEVindex(sort3);
  BPLEV = [BPLEV BPLEVindex];

% diagonalize the interaction matrix
[c,EF] = eig(H);             % Obtain matrix of eigenvectors and eigenvalus
[EF,ind] = sort(diag(EF));
c = c(:,ind');               % Rearrange the eigenvectors so that they
                             % match the sorted eigenvalues
EHFS = [EHFS EF'];
end

[EHFSSort,index] = sort(EHFS);
BPLEVSort = BPLEV(index);

FvalSort = Fval(index);

% In case of I=0, display energies and J quantun numbers for the
% Breit-Pauli levels.
% Otherwise, diplay hyperfine energies, J and F quantum numbers and which
% Breit-Pauli level the hyperfine level is derived from

if(I == 0)
  header = 'level    E_fs (a.u.)     J ';
  Jstring = rats(J',5);
else
  header = 'level    E_hfs (a.u.)   FS-LEV    J      F ';
  HFSJval = zeros(1,N_Feigvec);
  for i =1:N_Feigvec
    HFSJval(i) = J(BPLEV(i));
  end
  Jstring = rats(HFSJval',6);
  Fstring = rats(Fval',6);
end

blank = blanks(length(EHFS))';
disp(blanks(2)')
disp(header);
disp('-----------------------------------------------')
blank =  ' ';
levelstring = int2str((1:length(EHFS))');
EFstring = num2str(EHFS','%.9f');
BPLEVstring = num2str(BPLEV','%3u');
if(I == 0)
  disp([levelstring,blank(ones(length(EHFS),6)),EFstring,blank(ones(length(EHFS),2)),...
         Jstring])
else
  disp([levelstring blank(ones(length(EHFS),6)) EFstring blank(ones(length(EHFS),4)) ...
         BPLEVstring  blank(ones(length(EHFS),2)) Jstring Fstring])
end

% End display information about fine/hyperfine structure levels

% Make some sorting of the energy levels for the output file.

BPout = [];
Fout = [];
EHFSout = [];
for i = 1:max(BPLEVSort)
  BPdummy = [];
  Fdummy = [];
  EHFSdummy = [];
  for j=1:N_Feigvec
    if (i==BPLEVSort(j))
      BPdummy = [BPdummy i];
      Fdummy = [Fdummy FvalSort(j)];
      EHFSdummy = [EHFSdummy EHFSSort(j)];
    end
  end
  [SLASK INDX] = sort(-Fdummy);    % Sort with highest value first
  BPout = [BPout BPdummy(INDX)];
  Fout = [Fout Fdummy(INDX)];
  EHFSout = [EHFSout EHFSdummy(INDX)];
end

% Define plots
if ITtype~=1
  N_plots0 = N_plots+1;        % Counter for the number of plots.
  y = 1;
  Y = 1;
  n = 0;
  N = 0;
  moreplots = 1;
  moreplots = input("\nWould you like a plot of Zeeman splitting with B field? (Y/N) ");
  while moreplots == 1
     plotind = input("\nGive an index vector of the levels for which the\nzeeman patterns should be plotted ");
     plotind = unique(plotind);   % Omit levels occuring more than ones in the index vector.
     if any(plotind>length(EHFS)) | any(plotind<1)
        warning('Index is less than one or larger than the number of Breit-Pauli levels.')
     else
        N_plots = N_plots + 1;
        % If plotmatrix(n,m) = 1 then level m is ploted in figure n.
        % If plotmatrix(n,m) = 0 then level m is not ploted in figure n.
        plotmatrix(N_plots,:)=[plotind zeros(1,length(EHFS)-length(plotind))];
        if (I == 0)
          plottext(N_plots) = input("\nEnergies and the J-value printed in the plot ? (Y/N) ");
        else
          plottext(N_plots) = input("\nEnergies and the F-value printed in the plot ? (Y/N) ");
        end
        moreplots = input("\nMore plots ? (Y/N) ");
        % Calculate average energy for each plot.
        Eave(N_plots) = (2*Fval(plotind)+1)*EHFS(plotind)'/sum(2*Fval(plotind)+1);
     end
  end

  % Clear plots
  for n = N_plots0:N_plots
     figure(n)
     clf
  end
  % End clear plots
end %ITtype

% Print information about the different fine/hyperfine structure levels
% in the <name>.(c)zm file

if(I == 0)
  fprintf(fp,'Eigenv.  J-val    FS-energy \n\n');
  for i = 1:N_Feigvec
    fprintf(fp,'  %2u    ',i);
    fprintf(fp,'   %3.1f  ',J(i));
    fprintf(fp,' %.9f\n',EBP(i));
  end
else
  fprintf(fp,'Eigenv. J-val F-val   HFS-energy \n\n');
  for i = 1:N_Feigvec
    fprintf(fp,'  %2u    ',i);
    fprintf(fp,' %3.1f  ',J(BPout(i)));
    fprintf(fp,' %3.1f  ',Fout(i));
    fprintf(fp,' %.9f\n',EHFSout(i));
  end
end
fprintf(fp,'\n');


% Loop over MF values
% hfslev, hfslevINDEX and Ediag are used to print the mixing coeficients in the right
% possition in the output file for a certain MF-value

if ITtype ~= 1
  MFmax=Fmax;
  MFmin=-Fmax;
else
  MFmax=Fmin;
  MFmin=Fmin;
end

for MF = MFmax:-1:MFmin

  mdim = sum(Fval >= abs(MF));   % Find the size of the interaction matrix

  HDCB = zeros(mdim);
  Hhfs = zeros(mdim);
  HM = zeros(mdim);
  nrow = 0;
  hfslev = [];
  hfslevINDX = 0;
  Ediag = [];
  for i=1:N_eigvec
    FBMIN = abs(J(i)-I);
    FBMAX = J(i)+I;
    for iF = FBMAX:-1:FBMIN
      hfslevINDX = hfslevINDX + 1;
      if (iF >= abs(MF))
        hfslev = [hfslev hfslevINDX];
        ncol = 0;
        nrow = nrow + 1;
        for j=1:i
          FAMIN = abs(J(j)-I);
          FAMAX = J(j)+I;
          for jF = FAMAX:-1:FAMIN
            if (jF >= abs(MF))
              ncol = ncol + 1;
              % Hyperfine interaction matrix
              if (iF == jF)
                s = I+J(i)+iF;
                w6j = wig6j(I,J(j),iF,J(i),I,1);
                Hdipol = (-1)^s*w6j*sqrt((2*J(j)+1)*(2*I+1))* T(i,j,1)*M(1);
                w6j = wig6j(I,J(j),iF,J(i),I,2);
                Hquadro = (-1)^s*w6j*sqrt((2*J(j)+1)*(2*I+1))*T(i,j,2)*M(2);
                Hhfs(nrow,ncol) = Hdipol + Hquadro;
              % add diagonal J dependent energies
                if (i == j)
                  HDCB(nrow,ncol) = EBP(j);
                  Ediag = [Ediag HDCB(nrow,ncol)+Hhfs(nrow,ncol)];
                end
                Hhfs(ncol,nrow) = Hhfs(nrow,ncol);
              end
            % construct the magnetic part of the interaction matrix
            % assuming B = 1 (see eq 44 and 45)
              HM(nrow,ncol) = (-1)^(iF-MF)*wig3j(iF,1,jF,-MF,0,MF) ...
                     *(-1)^(I+J(j)+iF+1)* ...
                     sqrt(2*iF+1)*sqrt(2*jF+1) * sqrt(2*J(j)+1) * ...
                     wig6j(I,J(i),iF,1,jF,J(j))*NT(i,j);
              HM(ncol,nrow) = HM(nrow,ncol);
            end
          end
        end
      end
    end
  end
  [Ediag sortINDEX] = sort(Ediag);
  hfslev = hfslev(sortINDEX);

  % Loop over B and construct the total interaction matrix
  if ITtype~=1
    eigvalB = zeros(mdim,N_B);
    for nB = 1:N_B
      HTOT = HDCB + Hhfs + Bau(nB)*HM;
      eigvalB(:,nB) = sort(eig(HTOT)); % Sort the eigenvalues in ascending order.
    end
  else
    HTOT = HDCB + Hhfs;
  end

   % Print eigenvectors and eigenvalues to <name>.(c)zm for maximal value
   % of the B field
  [Evec,Eval] = eig(HTOT);        % Obtain eigenvectors and eigenvalus
  [Evals,ind] = sort(eig(HTOT));  % Sort the eigenvalues in ascending order
  [ndim1,ndim2] = size(Evec);
  Evecs = Evec(:,ind');           % Rearrange the eigenvectors so that they
                                   % match the sorted eigenvalues
  if ITtype ~=1
    fprintf(fp,'\n');
    fprintf(fp,'%s','  2*M = ');
    fprintf(fp,'%4.0f', 2*MF);
    fprintf(fp,'%s','  NUMBER = ');
    fprintf(fp,'%3.0f',ndim1);
    fprintf(fp,'\n');
    fprintf(fp,'\n');
  end

  for i = 1:ndim1

    % mix is a vector which is used to print the mixing coeficient in <name>.(c)zm
    mix = zeros(size(Fval));

    if max(Evecs(:,i)) < abs(min(Evecs(:,i)))
      Evecs(:,i) = -Evecs(:,i);
    end

    for j = 1:ndim1
      mix(hfslev(j)) = Evecs(sortINDEX(j),i);
    end

    if ITtype==1
      fprintf(fp,'\n');
      fprintf(fp,'%s','  2*F = ');
      fprintf(fp,'%4.0f', 2*Fout(hfslev(i)));
      fprintf(fp,'%s','  NUMBER = ');
      fprintf(fp,'%3.0f',1);
      fprintf(fp,'\n');
      fprintf(fp,'\n');
    end

    fprintf(fp,'%6d',hfslev(i));
    fprintf(fp,'%16.9f',Evals(i));
    fprintf(fp,'\n');


    for j = 1:N_Feigvec
      fprintf(fp,'%16.6E ',mix(j));
    end
    fprintf(fp,'\n');
  end

  % End print to file


  % Start plotting

  % Sort Breit-Pauli or hyperfine energies in acsending order
  % and keep the original ordering in vector ind.
  [EM,ind] = sort(EHFS(1:mdim));
  if ITtype ~=1
    if(N_plots > 0)
    % Set energy units in plot
      eigvalBplot = eigvalB;
      EHFSplot = EHFS;
      Eaveplot = Eave;
      if UNITE==1
        eigvalBplot = 2*109737.31*eigvalB;
        EHFSplot = 2*109737.31*EHFS;
        Eaveplot = Eave*2*109737.31;
      elseif UNITE == 2
        eigvalBplot = eigvalB/1.519829e-10;
        EHFSplot = EHFS/1.519829e-10;
        Eaveplot = Eave/1.519829e-10;
      end
    % End set units in plots
    % Plot energies
      for i = 1:mdim
        for n = N_plots0:N_plots
          if any(plotmatrix(n,:)==ind(i))
            figure(n)
            % Plot zeeman energy relative average energy.
            plot(B,eigvalBplot(i,:)-Eaveplot(n),'k')
            hold on
          end
        end
      end
    end
    % End plot energies.
  end %ITtype
end %MF

if ITtype~=1
if (N_plots > 0)
  for n = N_plots0:N_plots
    figure(n)
    % Text in plot
    if plottext(n) == 1
      clear EHFSFindex
      EHFSFindex = plotmatrix(n,find(plotmatrix(n,:)));
      blank = blanks(length(EHFSFindex))';
      EHFSstring = num2str(EHFS(EHFSFindex)','%.9f');
      Fstring = rats(Fval(EHFSFindex)');
      if (I == 0)
        EFtext = [blanks(5),'E_{FS}(a.u.)',blanks(13),'J'];
      else
        EFtext = [blanks(5),'E_{HFS}(a.u.)',blanks(11),'F'];
      end
      figtext = strvcat(EFtext,[EHFSstring,blank,Fstring]);
%      text(1/20,3/20,figtext,'sc')
      xc=xlim();
      yc=ylim();
      text(1/20*(xc(2)-xc(1))+xc(1),3/20*(yc(2)-yc(1))+yc(1),figtext)
    end
     % End text in plot

    if UNITB == 0
      xlabel('B (Tesla)')
    else
      xlabel('B (Gauss)')
    end

    if UNITE == 0
      ylabel('Energy (Hartree)')
    elseif UNITE == 1
      ylabel('Energy (cm^{-1})')
    else
      ylabel('Energy (MHz)')
    end
  end
  hold off
end
end %ITtype

display(' ')
display(['Finished ',name])
display(' ')

fclose(fp);
end
