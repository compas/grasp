%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                                %
%                                                                                                %

 function transm = mtrans(ITtype,I,j1,j2,f1,f2,m1,m2,TM,EMtype,qn,ptype)

%                                                                                                %
% This function computes transition matrix between M-sublevels.                                  %
%                                                                                                %
% Written by Wenxian Li, March 2019                                                              %
%                                                                                                %
% This m-file calls for: wig3j.m and wig6j.m                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

transm=0;

if ITtype==0     %mit-fs
  if ((EMtype==1 & ptype==0) | (EMtype==2 & ptype==1))
    transm=transm+(-1)^(j1-m1)*wig3j(j1,1,j2,-m1,qn,m2)*TM;
  elseif ((EMtype==2 & ptype==0) | (EMtype==1 & ptype==1))
    transm=transm+(-1)^(j1-m1)*wig3j(j1,2,j2,-m1,qn,m2)*TM;
  end
elseif ITtype==1  %hit
  if ((EMtype==1 & ptype==0) | (EMtype==2 & ptype==1))
    transm=transm+(-1)^(I+j2+f1+1)*sqrt((2*f1+1)*(2*f2+1))*wig6j(j1,f1,I,f2,j2,1)*TM;
  elseif ((EMtype==2 & ptype==0) | (EMtype==1 & ptype==1))
    transm=transm+(-1)^(I+j2+f1+2)*sqrt((2*f1+1)*(2*f2+1))*wig6j(j1,f1,I,f2,j2,2)*TM;
  end
elseif ITtype==2  %mit-hfs
  if ((EMtype==1 & ptype==0) | (EMtype==2 & ptype==1))
    transm=transm+(-1)^(f1-m1+j2+I+f1+1)*sqrt((2*f1+1)*(2*f2+1))*wig3j(f1,1,f2,-m1,qn,m2)*...
                    wig6j(j1,f1,I,f2,j2,1)*TM;
  elseif ((EMtype==2 & ptype==0) | (EMtype==1 & ptype==1))
    transm=transm+(-1)^(f1-m1+j2+I+f1+2)*sqrt((2*f1+1)*(2*f2+1))*wig3j(f1,2,f2,-m1,qn,m2)*...
                    wig6j(j1,f1,I,f2,j2,2)*TM;
  end
end
