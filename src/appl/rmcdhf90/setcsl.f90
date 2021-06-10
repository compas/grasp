!***********************************************************************
!                                                                      *
      SUBROUTINE SETCSL(NAME, NCORE, IDBLK)
!                                                                      *
!  A container which calls lib92/lodcsh to get                         *
!     ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()                *
!  It then calls lodcsl to load the entire file (all blocks),          *
!  obtaining   IQA, JQSA, JCUPA                                        *
!  Note that unlike setcsl of mcp[mpi] and rci[mpi]vu, the following   *
!  quantities are assumed to be known:                                 *
!     nblock, ncfblk(), idblk(), ncftot.                               *
!  In this rscfmpi, they have been read into memory from mcp files via *
!  a call to setmcp.                                                   *
!                                                                      *
!  Xinghong He 98-08-06                                                *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param,  ONLY: DOUBLE
      USE parameter_def,    ONLY: NNNW
      USE memory_man
      USE hblock_C
      USE def_C
      USE orb_C, NCFTOT=>NCF
      USE stat_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE openfl_I
      USE lodcsh_I
      USE lodcsh2gg_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NCORE
      CHARACTER  :: NAME*(*)
      CHARACTER  :: IDBLK(*)*8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER ::  IERR, IOS
      CHARACTER :: STR*15
!-----------------------------------------------
!
! opens, reads the header of the file

      CALL OPENFL (21, NAME, 'FORMATTED', 'OLD', IERR)
      IF (IERR == 1) THEN
         WRITE (6, *) 'Error when opening ', NAME(1:LEN_TRIM(NAME))
         ERROR STOP
      ENDIF

      READ (21, '(1A15)', IOSTAT=IOS) STR
      IF (IOS/=0 .OR. STR/='Core subshells:') THEN
         WRITE (6, *) 'Not a Configuration Symmetry List File;'
         CLOSE(21)
         ERROR STOP
      ENDIF

      !..Load header of <name> file
      CALL LODCSH (21, NCORE)

! Allocate memories for all blocks and then load the entire file

      CALL ALLOC (iqa, NNNW, NCFTOT, 'IQA', 'SETCSL')
!GG      CALL ALLOC (jqsa, NNNW,3,NCFTOT, 'JQSA', 'SETCSL')
!GG      CALL ALLOC (jcupa, NNNW, NCFTOT, 'JCUPA', 'SETCSL')

      CALL LODCSH2GG (21, NCORE, -119)
      ! -119 means "load all blocks"

      RETURN
      END SUBROUTINE SETCSL
