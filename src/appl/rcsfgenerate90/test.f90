! last edited October 31, 1996
      subroutine test(p1, p2, pop1, pop2, nmax) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06  
!...Switches:                     
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nmax 
      logical , intent(out) :: p1 
      logical , intent(out) :: p2 
      integer , intent(in) :: pop1(15,0:10,0:1) 
      integer , intent(in) :: pop2(15,0:10,0:1) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: n, l, k 
!-----------------------------------------------
 
      p1 = .TRUE. 
      p2 = .TRUE. 
      do n = 1, nmax 
         do l = 0, min(10,n - 1) 
            if (pop1(n,l,1) + pop1(n,l,0) < pop2(n,l,1) + pop2(n,l,0)) then 
               p1 = .FALSE. 
               return  
            else if (pop1(n,l,1) + pop1(n,l,0) > pop2(n,l,1) + pop2(n,l,0)) &
                  then 
               p2 = .FALSE. 
               return  
            else if (pop1(n,l,1) < pop2(n,l,1)) then 
               p1 = .FALSE. 
               return  
            else if (pop1(n,l,1) > pop2(n,l,1)) then 
               p2 = .FALSE. 
               return  
            else if (pop1(n,l,0) < pop2(n,l,0)) then 
               p1 = .FALSE. 
               return  
            else if (pop1(n,l,0) > pop2(n,l,0)) then 
               p2 = .FALSE. 
               return  
            endif 
         end do 
      end do 
      return  
      end subroutine test 
