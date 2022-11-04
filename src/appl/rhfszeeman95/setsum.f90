!************************************************************************
!*                                                                      *
       SUBROUTINE SETSUM(NAME,NCI,NOFFD)
!*                                                                      *
!*   Open the .gjhfs file on stream 111 and .ch on  29                  *
!*                                                                      *
!*   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
!*                                                                      *
!*   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
!*                                                                      *
!*   Updated by Per Jonsson                               28 Oct 1999   *
!*   Updated by Per and Martin                                          *
!*                                                                      *
!*   Translated by Wenxian Li F77 to F90 12/28/18                       *
!************************************************************************
!
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s 
!-----------------------------------------------
      USE openfl_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s 
!-----------------------------------------------
      INTEGER, INTENT(IN) :: NCI 
      INTEGER, INTENT(IN) :: NOFFD 
      CHARACTER, INTENT(IN) :: NAME*24 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s 
!-----------------------------------------------
      INTEGER   :: K, IERR 
      CHARACTER :: FILNAM1*256, FILNAM2*256, FORM*11, STATUS*3
!-----------------------------------------------
!
!   File <name>.gjhfs and <name>.h or <name>.cgjhfs and 
!   <name>.ch is FORMATTED
!
      K = INDEX(NAME,' ')
      IF (NCI == 0) THEN
         FILNAM1 = NAME(1:K-1)//'.ch'
         IF (NOFFD == 0) THEN
            FILNAM2 = NAME(1:K-1)//'.cgjhfs'
         ENDIF
      ELSE
         FILNAM1 = NAME(1:K-1)//'.h'
         IF (NOFFD == 0) THEN
            FILNAM2 = NAME(1:K-1)//'.gjhfs'
         ENDIF
      ENDIF
      FORM = 'FORMATTED'
      STATUS = 'NEW'
!
      CALL OPENFL (29, FILNAM1, FORM, STATUS, IERR)
      IF (IERR /= 0) THEN
         PRINT *, 'Error when opening',FILNAM1
         STOP
      ENDIF
!
      IF (NOFFD == 0) THEN
         CALL OPENFL (111,FILNAM2,FORM,STATUS,IERR)
         IF (IERR /= 0) THEN
            PRINT *, 'Error when opening',FILNAM2
            STOP
         ENDIF
      ENDIF
!
      RETURN
      END SUBROUTINE SETSUM
