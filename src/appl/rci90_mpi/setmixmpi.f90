!***********************************************************************
!                                                                      *
      SUBROUTINE SETMIXmpi(NAME, IDBLK)
!                                                                      *
!   Opens the  .mix  file on stream 25; writes a header to this file;  *
!   calls LODMIX to interactively determine the eigenpairs required.   *
!                                                                      *
!   Call(s) to: [LIB92]: OPENFL.                                       *
!               [RCI92]: LODMIX.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
!   Modified by Xinghong He               Last revision: 23 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   C O M M O N    B l o c k s
!-----------------------------------------------
      USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      USE lodmixmpi_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER(LEN=*) , INTENT(IN) :: NAME
      CHARACTER(LEN=8), DIMENSION(*) :: IDBLK
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER*11, PARAMETER :: FORM = 'UNFORMATTED'
      CHARACTER*7, PARAMETER  :: STATUS = 'UNKNOWN'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K
!-----------------------------------------------
!
      IF (myid .EQ. 0) THEN
         k = INDEX (name,' ')
         CALL openfl (25, name(1:k-1)//'.cm', FORM, STATUS, ierr)
         IF (ierr .NE. 0) THEN
            CALL stopmpi ('setmix: Error when opening .cm file', myid)
         ENDIF

         WRITE (25) 'G92MIX'
      ENDIF

      CALL LODMIXmpi (IDBLK)

      RETURN
      END SUBROUTINE SETMIXmpi
