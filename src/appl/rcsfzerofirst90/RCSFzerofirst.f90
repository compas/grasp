!***********************************************************************
!                                                                      *
      PROGRAM RCSFzerofirst
!-----------------------------------------------
!                                                                      *
!   From a set of CSLs this program identifies the ones that           *
!   interact with a given multireference                               *
!                                                                      *
!   This program is a slight modification of the GENMCP program        *
!                                                                      *
!   Written by  G. Gaigalas                     Vilnius, May 2016      *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE BLK_C,            only: NBLOCK,NCFBLK
      USE rang_Int_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE set_CSF_ZFlist_I
      USE lodcsl_Zero_I
      USE lodcsl_Part_I
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL :: NEXT_BLOCK
      INTEGER :: CSF_Number, ncount1
!-----------------------------------------------
      call starttime (ncount1, 'RCSFzerofirst')
!
      print *, ""
      print *, "RCSFzerofirst: Takes a list of CSFs and partitions each symmetry"
      print *, "               block into a zero- and first-order CSF space from"
      print *, "               a zero-order list."
      print *, "               (C)   Copyright by G. Gaigalas and Ch. F. Fischer"
      print *, "               (Fortran 95 version)                NIST  (2017)."
      print *, "               Input files:     list with CSFs to be partitioned"
      print *, "                                list with CSFs defining "
      print *, "                                            the zero-order space"
      print *, "               Output file:     rcsf.out"
      print *, ""
      NBLOCK = 0
      CALL SET_CSF_ZFlist
!cychen: output the zero-order space for re-use in mcp, rci
      open(301,file='icut',form='formatted',status='unknown')
      WRITE (6, *) "  Block    Zero-order Space   Complete Space"
      DO
         CALL LODCSL_Zero (NEXT_BLOCK)
         CALL LODCSL_Part (CSF_Number)
         WRITE (6,'(3X,I2,6X,I14,3X,I17)')                            &
                                     NBLOCK,NCFBLK(NBLOCK),CSF_Number-1
         write(301,*)NCFBLK(NBLOCK)
         deallocate (Found)
         deallocate (C_shell)
         deallocate (C_quant)
         deallocate (C_coupl)
         IF(.NOT. NEXT_BLOCK) EXIT
         WRITE(22,'(A2)') ' *'
      END DO
      close(301)
      call stoptime (ncount1, 'RCSFzerofirst')
      STOP
      END PROGRAM RCSFzerofirst
