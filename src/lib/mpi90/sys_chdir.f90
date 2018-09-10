!***********************************************************************
      subroutine sys_chdir (dir, lendir, ierr)
!
!   This routine changes current working directory to dir.
!   lendir is the length of character string dir; 
!   ierr will be zero if successful, otherwise non-zero;
!
!   Xinghong He  98-08-21
!
!************************************************************************     
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      implicit none
      character(len=*), intent(in):: dir
      integer, intent(in):: lendir
      integer, intent(out):: ierr


      integer chdir

      ierr = chdir (dir(1:lendir))

      return
      end subroutine sys_chdir
