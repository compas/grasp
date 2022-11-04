%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                %
%                                                                                                %

 function [AM_Mall,n] = trans(name_i,name_f,relcalc,ITtype,JE_i,FE_i,JE_f,FE_f,Bmax, ...
                              EM_i,EM_f,Parity_i,Parity_f,UNITB,I,mu,Q,Ttype)

%                                                                                                %
% This function opens the file <name1>.<name2>.(c)t and reads the reduced transition matrix.     %
% With the saved mixing coefficients between the levels the transition rates between magnetic    %
% fine- and hyperfine structure sublevels in the presence of external magnetic field and the     %
% rates of hyperfine induced transitions in the field free limit are computed and saved in the   %
% output file <name1>.<name2>.<tr>.mtrans.                                                       %
%                                                                                                %
% Written by Wenxian Li, March 2019                                                              %
%                                                                                                %
% This m-file calls for: mtrans.m                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open and read transition data from file <name_i>.<name_j>.ct
rightansw = 0;
while rightansw == 0
  if (relcalc == 'y' | relcalc == 'Y')
    file = strcat(name_i,'.',name_f,'.ct');
    rightansw = 1;
  elseif (relcalc == 'n' | relcalc == 'N')
    file = strcat(name_i,'.',name_f,'.t');
    rightansw = 1;
  end
  if (rightansw == 0)
    relcalc = input('You have to answere (Y/N)','s');
  end
end

i=1;
tr=zeros;

fp = fopen(file,'r');
a=fscanf(fp,'%s',9);
while ~feof(fp)
  tline=fgetl(fp);
  if ~isempty(strfind(tline, 'f')) &  ~isempty(strfind(tline, 'C')) & isempty(strfind(tline, 'gf'))
    [f2 lev_2 J_2 P_2 f1 lev_1 J_1 P_1 dE C A_C gf_C S_C M_C]=strread(tline,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s');
    if cell2mat(P_2)==Parity_f & cell2mat(P_1)==Parity_i
      tline=fgetl(fp);
      [B A_B gf_B S_B M_B]=strread(tline,'%s %s %s %s %s');
      tr(i,1)=JE_f(intersect(find(JE_f(:,2)==str2num(cell2mat(J_2))),find(JE_f(:,4)==str2num(cell2mat(lev_2)))),1);
      tr(i,2)=2*str2num(cell2mat(J_2));
      tr(i,3)=JE_f(intersect(find(JE_f(:,2)==str2num(cell2mat(J_2))),find(JE_f(:,4)==str2num(cell2mat(lev_2)))),3);
      tr(i,4)=JE_i(intersect(find(JE_i(:,2)==str2num(cell2mat(J_1))),find(JE_i(:,4)==str2num(cell2mat(lev_1)))),1);
      tr(i,5)=2*str2num(cell2mat(J_1));
      tr(i,6)=JE_i(intersect(find(JE_i(:,2)==str2num(cell2mat(J_1))),find(JE_i(:,4)==str2num(cell2mat(lev_1)))),3);
      tr(i,7)=str2num(cell2mat(A_B));
      tr(i,8)=str2num(cell2mat(M_B));
      tr(i,9)=str2num(cell2mat(A_C));
      tr(i,10)=str2num(cell2mat(M_C));
      i=i+1;
    end
  elseif ~isempty(strfind(tline, 'f')) & ~isempty(strfind(tline, 'M')) & isempty(strfind(tline, 'gf'))
    [f2 lev_2 J_2 P_2 f1 lev_1 J_1 P_1 dE M A_M gf_M S_M M_M]=strread(tline,'%s %s %s %s %s %s %s %s %s %s %s %s %s %s');
    if cell2mat(P_2)==Parity_f & cell2mat(P_1)==Parity_i
      tr(i,1)=JE_f(intersect(find(JE_f(:,2)==str2num(cell2mat(J_2))),find(JE_f(:,4)==str2num(cell2mat(lev_2)))),1);
      tr(i,2)=2*str2num(cell2mat(J_2));
      tr(i,3)=JE_f(intersect(find(JE_f(:,2)==str2num(cell2mat(J_2))),find(JE_f(:,4)==str2num(cell2mat(lev_2)))),3);
      tr(i,4)=JE_i(intersect(find(JE_i(:,2)==str2num(cell2mat(J_1))),find(JE_i(:,4)==str2num(cell2mat(lev_1)))),1);
      tr(i,5)=2*str2num(cell2mat(J_1));
      tr(i,6)=JE_i(intersect(find(JE_i(:,2)==str2num(cell2mat(J_1))),find(JE_i(:,4)==str2num(cell2mat(lev_1)))),3);
      tr(i,7)=str2num(cell2mat(A_M));
      tr(i,8)=str2num(cell2mat(M_M));
      i=i+1;
    end
  end
end
fclose(fp);
%end read transition data

%transition type
if Parity_f==Parity_i
  ptype=1;
else
  ptype=0;
end
%end transition type

%calculate and print the magnetic-field(MIT) or hyperfine(HIT) induced transition rates
%Define output file name
if ITtype==0
  fileM = strcat(name_i,'.',name_f,'.fs.mit.mtrans'); % MIT-isotopes with nuclear spin I=0
elseif ITtype==1
  fileM = strcat(name_i,'.',name_f,'.hfs.hit.trans'); % HIT-with B=0
elseif ITtype==2
  fileM = strcat(name_i,'.',name_f,'.hfs.mit.mtrans'); % MIT-isotopes with nuclear spin I~=0
end
fpM = fopen(fileM,'w');
%End define

% Read and print magnetic field (unit), nuclear paremeters and energy levels
if ITtype~=0
  if ITtype==2
    fprintf(fpM,'%s\n','Magnetic field and nuclear data');
    fprintf(fpM,'%s','B                              ');
    fprintf(fpM,'%15.7f', Bmax);
    if UNITB == 0
      fprintf(fpM,'%s\n',' T');
    else
      fprintf(fpM,'%s\n',' ');
    end
  elseif ITtype==1
    fprintf(fpM,'%s\n','Nuclear data');
  end
  fprintf(fpM,'%s%f%s\n','Nuclear spin                         ',I,' au');
  fprintf(fpM,'%s%f%s\n','Nuclear magnetic dipole moment       ',mu,' n.m.');
  fprintf(fpM,'%s%f%s\n','Nuclear electric quadrupole moment   ',Q,' barns');
  fprintf(fpM,'\n');
  fprintf(fpM,'\n');
elseif ITtype==0
  fprintf(fpM,'%s\n','Magnetic field');
  fprintf(fpM,'%s','  B  = ');
  fprintf(fpM,'%15.7f', Bmax);
  if UNITB == 0
      fprintf(fpM,'%s',' Tesla');
  else
      fprintf(fpM,'%s',' Gauss');
  end
  fprintf(fpM,'\n');
  fprintf(fpM,'\n');
end

if ITtype==0
  fprintf(fpM,'%s\n','Fine structure energies in a.u.');
  fprintf(fpM,'%s\n',name_i);
  fprintf(fpM,'%s\t%s\t%s\n','level','J','E_BP (a.u.)');
  for i=1:size(JE_i,1)
    fprintf(fpM,'%d\t%4.1f\t%f\n',JE_i(i,1:3));
  end
  fprintf(fpM,'\n');
  fprintf(fpM,'%s\n',name_f);
  fprintf(fpM,'%s\t%s\t%s\n','level','J','E_BP (a.u.)');
  for i=1:size(FE_f,1)
    fprintf(fpM,'%d\t%4.1f\t%f\n',JE_f(i,1:3));
  end
else
  fprintf(fpM,'%s\n','Hyperfine structure energies in a.u.');
  fprintf(fpM,'%s\n',name_i);
  fprintf(fpM,'%s\t%s\t%s\t%s\t%s\n','level','J','F','E_hfs (a.u.)','FS-LEV');
  for i=1:size(FE_i,1)
    fprintf(fpM,'%d\t%4.1f\t%4.1f\t%f\t\t%d\n',FE_i(i,1:5));
  end
  fprintf(fpM,'\n');
  fprintf(fpM,'%s\n',name_f);
  fprintf(fpM,'%s\t%s\t%s\t%s\t%s\n','level','J','F','E_hfs (a.u.)','FS-LEV');
  for i=1:size(FE_f,1)
    fprintf(fpM,'%d\t%4.1f\t%4.1f\t%f\t\t%d\n',FE_f(i,1:5));
  end
end
% End read and print magnetic field (unit), nuclear paremeters and energy levels

fprintf(fpM,'\n');
fprintf(fpM,'\n');
fprintf(fpM,'%s\n','Transition rates and wavelength in Kays');
fprintf(fpM,'\t\t\t%s\t\t\t\t\t\t\t%s\n','Upper','Lower');
if ITtype==0
  fprintf(fpM,'%s\t%s\t%s\t%s\t%s\t','level',' J',' M_J','    E_hfs (a.u.)','FS-LEV');
  fprintf(fpM,'%s\t%s\t%s\t%s\t%s\t','level',' J',' M_J','    E_hfs (a.u.)','FS-LEV');
  fprintf(fpM,'%s\t%s\n',' A (s-1)','   E (Kays)');
elseif ITtype==1
  fprintf(fpM,'%s\t%s\t%s\t%s\t%s\t','level',' J',' F','    E_hfs (a.u.)','FS-LEV');
  fprintf(fpM,'%s\t%s\t%s\t%s\t%s\t','level',' J',' F','    E_hfs (a.u.)','FS-LEV');
  fprintf(fpM,'%s\t%s\n',' A (s-1)','   E (Kays)');
elseif ITtype==2
  fprintf(fpM,'%s\t%s\t%s\t%s\t%s\t%s\t','level',' J',' F',' M_F','    E_hfs (a.u.)','FS-LEV');
  fprintf(fpM,'%s\t%s\t%s\t%s\t%s\t%s\t','level',' J',' F',' M_F','    E_hfs (a.u.)','FS-LEV');
  fprintf(fpM,'%s\t%s\n',' A (s-1)','   E (Kays)');
end

% Display energies and quantun numbers for the transition levels.
%In case of I=0, display energies and J quantun numbers for the Breit-Pauli levels.
% Otherwise, diplay hyperfine energies, J and F quantum numbers and which
% Breit-Pauli level the hyperfine level is derived from
J_i = FE_i(:,2);
J_f = FE_f(:,2);
Jstring_i = rats(J_i);
levelstring_i = int2str((1:size(FE_i,1))');
Jstring_f = rats(J_f);
levelstring_f = int2str((1:size(FE_f,1))');

blank = blanks(length(J_i))';
disp(blanks(2)')
if ITtype==0
  disp('level    E_BP (a.u.)         J'); %
  disp('-------------------------------')
  blank =  ' ';
  EBP_i = FE_i(:,4);
  EBP_f = FE_f(:,4);
  EBPstring_i = num2str(EBP_i,'%.9f');
  EBPstring_f = num2str(EBP_f,'%.9f');
  display('Initial levels:')
  disp([levelstring_i,blank(ones(length(EBP_i),6)),EBPstring_i,blank(ones(length(EBP_i),2)),Jstring_i])
  display('Final levels:')
  disp([levelstring_f,blank(ones(length(EBP_f),6)),EBPstring_f,blank(ones(length(EBP_f),2)),Jstring_f])
else
  disp('level    E_hfs (a.u.)     FS-LEV        J           F');
  disp('------------------------------------------------------------')
  blank =  ' ';
  F_i=FE_i(:,3);
  EHFS_i=FE_i(:,4);
  FS_i=FE_i(:,5);
  Fstring_i=rats(F_i);
  EHFSstring_i = num2str(EHFS_i,'%.9f');
  FSstring_i=rats(FS_i);
  F_f=FE_f(:,3);
  EHFS_f=FE_f(:,4);
  FS_f=FE_f(:,5);
  Fstring_f=rats(F_f);
  EHFSstring_f = num2str(EHFS_f,'%.9f');
  FSstring_f=rats(FS_f);
  display('Initial levels:')
  disp([levelstring_i,blank(ones(length(EHFS_i),4)),EHFSstring_i,FSstring_i,Jstring_i,Fstring_i])
  display('Final levels:')
  disp([levelstring_f,blank(ones(length(EHFS_f),4)),EHFSstring_f,FSstring_f,Jstring_f,Fstring_f])
end
%end display energies and J quantun numbers for the transition levels.

%Input the initial level of the computed transitions
MITind_ilist = input("\nGive an index vector of the initial levels(lower level):  ");
MITind_flist = input("\nGive an index vector of the final levels(upper level):  ");

AM_Mall=[];
for ii = 1:size(MITind_ilist,2)
  MITind_i=MITind_ilist(ii);

  while ~(any(FE_i(:,1)==MITind_i))
    MITind_i = input("\nThe level you give is not belonging to initial state, redo it):  ");
  end

  n=1;
  MITJ_i=FE_i(MITind_i,2);
  MITF_i=FE_i(MITind_i,3);

  % Loop upper states and print to output files
  for ff=1:size(MITind_flist,2)
    MITind_f=MITind_flist(ff);
    AM=[];
    AM_M=[];
    MITF_f=FE_f(MITind_f,3);
    % Extract perturbers index
    if ITtype==0
      MITind_mix=[find(FE_f(:,2)==MITF_f+1);find(FE_f(:,2)==MITF_f-1);find(FE_f(:,2)==MITF_f)]; % include d0
      MITind_mix(:,2)=MITind_mix(:,1);
    elseif ITtype==1
      MITind_mix=[find(FE_f(:,3)==MITF_f)];%include d0
      MITind_mix(:,2)=FE_f(MITind_mix,5);
    elseif ITtype==2
      MITind_mix=[find(FE_f(:,3)==MITF_f+1);find(FE_f(:,3)==MITF_f-1);find(FE_f(:,3)==MITF_f)];%include d0
      MITind_mix(:,2)=FE_f(MITind_mix,5);
    end
    % End extract perturbers index

    % Calculate transition rates
    j=1;
    for m=1:size(MITind_mix,1)
        MITJ_mix(m,1)=FE_f(MITind_mix(m,1),2);
        MITJ_mix(m,2)=FE_f(MITind_mix(m,1),3);
    end

    if ITtype~=1
      MITF_fmin=-MITF_f;
      MITF_fmax=MITF_f;
      MITF_imin=-MITF_i;
      MITF_imax=MITF_i;
    elseif ITtype==1
      MITF_fmin=MITF_f;
      MITF_fmax=MITF_f;
      MITF_imin=MITF_i;
      MITF_imax=MITF_i;
    end

    for MITM_f=MITF_fmax:-1:MITF_fmin %for
      sum_mL_EM=zeros((2*MITF_f+1)*(2*MITF_i+1),1);
      sum_mL_MM=zeros((2*MITF_f+1)*(2*MITF_i+1),1);
      for i=1:size(MITind_mix,1)
        MITind_F_i=FE_i(find(FE_i(:,1)==MITind_i),5);
        TM_L_temp=tr(intersect(find(tr(:,1)==MITind_mix(i,2)),find(tr(:,4)==MITind_F_i)),8);
        TM_V_temp=tr(intersect(find(tr(:,1)==MITind_mix(i,2)),find(tr(:,4)==MITind_F_i)),10);
        MC=EM_f(intersect(find(EM_f(:,1)==MITind_f),find(EM_f(:,3)==2*MITM_f)),MITind_mix(i,1)+3);
        for k=1:size(TM_L_temp,1)
          TM_L=TM_L_temp(k,:);  % Transition matrix element-B gauge
          TM_V=TM_V_temp(k,:);  % Transition matrix element-C gauge
          % B gauge is used for calculation
          iMi=1;
          if ~isempty(TM_L) & ~isempty(TM_V) & TM_V~=0
            EMtype=1;
            for MITM_i=MITF_imax:-1:MITF_imin
              MITind_M=iMi+(j-1)*(2*MITF_i+1);
              sum_mL_EM(MITind_M)=sum_mL_EM(MITind_M)+MC*mtrans(ITtype,I,MITJ_i,MITJ_mix(i,1),MITF_i,MITJ_mix(i,2), ...
                                                                MITM_i,MITM_f,TM_L,EMtype,MITM_i-MITM_f,ptype);
              sum_mL_MM(MITind_M)=sum_mL_MM(MITind_M);
              iMi=iMi+1;
            end
          elseif ~isempty(TM_L) & ~isempty(TM_V) & TM_V==0
            EMtype=2;
            for MITM_i=MITF_imax:-1:MITF_imin
              MITind_M=iMi+(j-1)*(2*MITF_i+1);
              sum_mL_EM(MITind_M)=sum_mL_EM(MITind_M);
              sum_mL_MM(MITind_M)=sum_mL_MM(MITind_M)+MC*mtrans(ITtype,I,MITJ_i,MITJ_mix(i,1),MITF_i,MITJ_mix(i,2), ...
                                                                MITM_i,MITM_f,TM_L,EMtype,MITM_i-MITM_f,ptype);
              iMi=iMi+1;
            end
          else
            for MITM_i=MITF_imax:-1:MITF_imin
              MITind_M=iMi+(j-1)*(2*MITF_i+1);
              sum_mL_EM(MITind_M)=sum_mL_EM(MITind_M);
              sum_mL_MM(MITind_M)=sum_mL_MM(MITind_M);
              iMi=iMi+1;
            end
          end
        end
      end

      E_f=EM_f(intersect(find(EM_f(:,1)==MITind_f),find(EM_f(:,3)==2*MITM_f)),2); % Upper level energy
      AM(j,1)=MITM_f;
      AM(j,2)=0;
      if Ttype==0 %if Ttype
        iMi=1;
        for MITM_i=MITF_imax:-1:MITF_imin
          E_i=EM_i(intersect(find(EM_i(:,1)==MITind_i),find(EM_i(:,3)==2*MITM_i)),2); % Lower level energy
          MITind_M=iMi+(j-1)*(2*MITF_i+1);
          AM_M(MITind_M,1:3)= FE_f(MITind_f,1:3);
          AM_M(MITind_M,4)=MITM_f;
          AM_M(MITind_M,5)=E_f;
          AM_M(MITind_M,6)=FE_f(MITind_f,5);
          AM_M(MITind_M,7:9)=FE_i(MITind_i,1:3);
          AM_M(MITind_M,10)=MITM_i;
          AM_M(MITind_M,11)=E_i;
          AM_M(MITind_M,12)=FE_i(MITind_i,5);
          if ITtype==1
            AM_M(MITind_M,13)=2.1420E10/(2*MITF_f+1)*abs((E_f-E_i))^3*sum_mL_EM(MITind_M)^2+ ...
                              7.5926E-1/(2*MITF_f+1)*abs((E_f-E_i))^5*sum_mL_MM(MITind_M)^2;
          else
            AM_M(MITind_M,13)=2.1420E10*abs((E_f-E_i))^3*sum_mL_EM(MITind_M)^2+ ...
                              7.5926E-1*abs((E_f-E_i))^5*sum_mL_MM(MITind_M)^2;
          end
          AM_M(MITind_M,14)=abs(AM_M(MITind_M,5)-AM_M(MITind_M,11))*219474.63068;
          iMi=iMi+1;
          E_i=FE_i(find(FE_i(:,1)==MITind_i),4); % Lower level energy   %used for AM
          AM(j,2)=AM(j,2)+AM_M(MITind_M,13);
        end
      elseif Ttype==1
         iMi=1;
         for MITM_i=MITF_imax:-1:MITF_imin
           E_i=EM_i(intersect(find(EM_i(:,1)==MITind_i),find(EM_i(:,3)==2*MITM_i)),2); % Lower level energy
           MITind_M=iMi+(j-1)*(2*MITF_i+1);
           AM_M(MITind_M,1:3)= FE_f(MITind_f,1:3);
           AM_M(MITind_M,4)=MITM_f;
           AM_M(MITind_M,5)=E_f;
           AM_M(MITind_M,6)=FE_f(MITind_f,5);
           AM_M(MITind_M,7:9)=FE_i(MITind_i,1:3);
           AM_M(MITind_M,10)=MITM_i;
           AM_M(MITind_M,11)=E_i;
           AM_M(MITind_M,12)=FE_i(MITind_i,5);
           if ITtype==1
             AM_M(MITind_M,13)=5.7032E4/(2*MITF_f+1)*abs((E_f-E_i))^5*sum_mL_EM(MITind_M)^2+ ...
                               2.8516E5/(2*MITF_f+1)*abs((E_f-E_i))^3*sum_mL_MM(MITind_M)^2;
           else
             AM_M(MITind_M,13)=5.7032E4*abs((E_f-E_i))^5*sum_mL_EM(MITind_M)^2+ ...
                               2.8516E5*abs((E_f-E_i))^3*sum_mL_MM(MITind_M)^2;
           end
           AM_M(MITind_M,14)=abs(AM_M(MITind_M,5)-AM_M(MITind_M,11))*219474.63068;
           iMi=iMi+1;
           E_i=FE_i(find(FE_i(:,1)==MITind_i),4); % Lower level energy   %used for AM
           AM(j,2)=AM(j,2)+AM_M(MITind_M,13);
         end
      end % end if Ttype
      j=j+1;
    end %for
    % End calculate transition rates


    % Print MIT or HIT results
    if ~isempty(AM)
      iMi=1;
      for k=1:size(AM,1)
        %print to <name>.trans
        if ITtype==0
          iMi=1;
          for MITM_i=MITJ_i:-1:-MITJ_i
            MITind_M=iMi+(k-1)*(2*MITJ_i+1);
            Ekays=abs((AM_M(MITind_M,5)-AM_M(MITind_M,4)))*219474.63068;
            fprintf(fpM,'%2.0f\t%3.1f\t%3.1f\t%16.9f\t%2.0f\t%2.0f\t%3.1f\t%3.1f\t%16.9f\t%2.0f\t%5.4E\t%16.9f',...
                    AM_M(MITind_M,1:2),AM_M(MITind_M,4:8),AM_M(MITind_M,10:14));
            iMi=iMi+1;
            fprintf(fpM,'\n');
          end
        elseif ITtype==1
          Ekays=abs((FE_f(MITind_f,4)-FE_i(MITind_i,4)))*219474.63068;
          fprintf(fpM,'%2.0f\t%3.1f\t%3.1f\t%16.9f\t%2.0f\t%2.0f\t%3.1f\t%3.1f\t%16.9f\t%2.0f\t%5.4E\t%16.9f',...
             AM_M(MITind_M,1:3),AM_M(MITind_M,5:9),AM_M(MITind_M,11:14));
          fprintf(fpM,'\n');
        elseif ITtype == 2
          iMi=1;
          for MITM_i=MITF_imax:-1:MITF_imin
            MITind_M=iMi+(k-1)*(2*MITF_i+1);
            Ekays=abs((AM_M(MITind_M,5)-AM_M(MITind_M,4)))*219474.63068;
            fprintf(fpM,'%2.0f\t%3.1f\t%3.1f\t%3.1f\t%16.9f\t%2.0f\t%2.0f\t%3.1f\t%3.1f\t%3.1f\t%16.9f\t%2.0f\t%5.4E\t%16.9f',...
                             AM_M(MITind_M,1:14));
            iMi=iMi+1;
            fprintf(fpM,'\n');
          end
        end % end ITtype
      end % end for size(AM)
    end % End print MIT or HIT results

    AM_Mall=[AM_Mall;AM_M];

  end % Loop and sum the contributions from different perturbers

end

fclose(fpM);

end
%end function
