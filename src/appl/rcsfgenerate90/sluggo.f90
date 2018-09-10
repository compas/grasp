 
!     last edited September 23, 1995
      subroutine sluggo(i, j, varmax, varupp, varned, ansats, org, lock, low, &
         start, stopp) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: i 
      integer , intent(in) :: j 
      integer , intent(in) :: varmax 
      integer , intent(out) :: start 
      integer , intent(out) :: stopp 
      logical , intent(in) :: lock 
      integer , intent(inout) :: varupp(15,0:10) 
      integer , intent(inout) :: varned(15,0:10) 
      integer , intent(in) :: ansats(15,0:10,0:1) 
      integer , intent(in) :: org(15,0:10) 
      integer , intent(in) :: low(15,0:10) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: minmax, iold, jold 
!-----------------------------------------------
      if (i == 1) then 
         varupp(1,0) = 0 
         varned(1,0) = 0 
      else 
         if (j == 0) then 
            iold = i - 1 
            jold = min(10,iold - 1) 
         else 
            iold = i 
            jold = j - 1 
         endif 
         if (jold == 0) then 
            varupp(i,j) = varupp(iold,jold) + max(0,ansats(iold,jold,0)-org(&
               iold,jold)) 
            varned(i,j) = varned(iold,jold) + max(0,org(iold,jold)-ansats(iold,&
               jold,0)) 
         else 
            varupp(i,j) = varupp(iold,jold) + max(0,ansats(iold,jold,0)+ansats(&
               iold,jold,1)-org(iold,jold)) 
            varned(i,j) = varned(iold,jold) + max(0,org(iold,jold)-ansats(iold,&
               jold,0)-ansats(iold,jold,1)) 
         endif 
      endif 
      if (lock) then 
         start = org(i,j) 
         stopp = org(i,j) 
         return  
      endif 
      if (j >= 5) then 
         minmax = 4 
      else 
         minmax = 4*j + 2 
      endif 
      start = min(minmax,org(i,j)+(varmax-varupp(i,j))) 
      stopp = max(low(i,j),org(i,j)-(varmax-varned(i,j))) 
      return  
      end subroutine sluggo 
