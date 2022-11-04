%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                %
%                                                                                                %

function [N_eigvec,JE,FE,B,unitB,EM,Parity] = mixingC(name,relcalc,ITtype,I)

%                                                                                                %
% This function is to read eigenvectors and mixing coefficients from <name>.(c)gjhfs             %
% and <name>.(c)zm files.                                                                        %
%                                                                                                %
% Written by Wenxian Li, March 2019                                                              %
%                                                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Read Breit-Pauli energy levels
rightansw = 0;
while rightansw == 0
  if (relcalc == 'y' | relcalc == 'Y')
    file = strcat(name,'.cgjhfs');
    rightansw = 1;
  elseif (relcalc == 'n' | relcalc == 'N')
    file = strcat(name,'.gjhfs');
    rightansw = 1;
  end
  if (rightansw == 0)
    relcalc = input('You have to answere (Y/N)','s');
  end
end
fp = fopen(file,'r');
[string1,c] = fscanf(fp,'%s',4);
[N_eigvec,c] = fscanf(fp,'%u',1);
[string2,c] = fscanf(fp,'%s',4);
[JEtemp,c] = fscanf(fp,'%f%f%s%f',[4,N_eigvec]);
Parity=JEtemp(3,1);
JE(1:N_eigvec,1)=1:N_eigvec;
JE(:,2) = JEtemp(2,:)';
JE(:,3) = JEtemp(4,:)';
JE(:,4) = JEtemp(1,:)';

[string3,c] = fscanf(fp,'%s',3);
[NT,c] = fscanf(fp,'%e',[N_eigvec,N_eigvec]);
[string3,c] = fscanf(fp,'%s',6);
[T1,c] = fscanf(fp,'%e',[N_eigvec,N_eigvec]);
[string3,c] = fscanf(fp,'%s',6);
[T2,c] = fscanf(fp,'%e',[N_eigvec,N_eigvec]);
T1=T1';
T2=T2';

N_eigvecJ=N_eigvec;
fclose(fp);
% End read Breit-Pauli energy levels

% Read hyperfine energy levels
rightansw = 0;
while rightansw == 0
  if (relcalc == 'y' | relcalc == 'Y')
    file = strcat(name,'.czm');
    rightansw = 1;
  elseif (relcalc == 'n' | relcalc == 'N')
    file = strcat(name,'.zm');
    rightansw = 1;
  end
  if (rightansw == 0)
    relcalc = input('You have to answere (Y/N)','s');
  end
end

fp = fopen(file,'r');
[string3,c]=fscanf(fp,'%s',2);
B=fscanf(fp,'%f');
unitB=fscanf(fp,'%s',1);
[string,c] = fscanf(fp,'%s',2);
[N_eigvec,c] = fscanf(fp,'%u',1);

if ITtype~=0
  [string,c]=fscanf(fp,'%s',4);
  [FEtemp,c] = fscanf(fp,'%f%f%f%f',[4,N_eigvec]);
  FE(:,1) = FEtemp(1,:)';
  FE(:,2) = FEtemp(2,:)';
  FE(:,3) = FEtemp(3,:)';
  FE(:,4) = FEtemp(4,:)';
  i=1;
  j=1;
  while i<=size(FE,1)
    n=abs(I+FE(i,2))-abs(I-FE(i,2))+1;
    FE(i:i+n-1,5)=JE(j,1);
    FE(i:i+n-1,6)=JE(j,4);
    i=i+n;
    j=j+1;
  end
else
  FE(:,1:2)=JE(:,1:2);
  FE(:,3)=JE(:,2);
  FE(:,4)=JE(:,3);
  FE(:,5)=JE(:,1);
  FE(:,6)=JE(:,4);
end
% End read hyperfine energy levels

%Read and reformat the <name>.(c)zm file and save the mixing
%coefficients to EM
i=1;
EM=zeros;
while ~feof(fp)
  tline=fgetl(fp);
  if ~isempty(strfind(tline, '='))
    [string string M]=strread(tline,'%s %s %f');
    for j=1:M(2,1)
      temp_EM=fscanf(fp,'%f',2);
      EM(i,1:2)=temp_EM';
      EM(i,3)=M(1,1);
      EM(i,4:N_eigvec+3)=fscanf(fp,'%f',[N_eigvec,1]);
      i=i+1;
    end
  end
end

% Sort the rows of EM according to the index vector
EM=sortrows(EM,1);
%End read and reformat the <name>.(c)zm file

end
%end function
