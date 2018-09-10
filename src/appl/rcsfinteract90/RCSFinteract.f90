!***********************************************************************
!                                                                      *
      PROGRAM RCSFinteract
!-----------------------------------------------
!                                                                      *
!   From a set of CSLs this program identifies the ones that           *
!   interact with a given multireference                               *
!                                                                      *
!   This program is a slight modification of the GENMCP program        *
!                                                                      *
!   Written by  G. Gaigalas                   NIST, December 2015      *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE memory_man
      USE default_C 
      USE BLK_C,            only:  NBLOCK,NCFBLK
      USE orb_C
      USE STAT_C
      USE rang_Int_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE Interact_MR_I
      USE set_CSF_list_I 
      USE lodcsl_MR_I
      USE lodcsl_CSF_I
      USE factt_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      LOGICAL ::  NEXT_BLOCK, NEXT_CSF
      INTEGER :: ICOLBREI,NCORE_MR,NPEEL_MR,NCORE,NPEEL,CSF_Number,NCFD
      INTEGER :: I_Count, ncount1
!-----------------------------------------------
      call starttime (ncount1, 'RCSFinteract')
!
      print *, ""
      print *, "RCSFinteract: Determines all the CSFs (rcsf.inp) that interact"
      print *, "              with the CSFs in the multireference (rcsfmr.inp)"
      print *, "              (C)  Copyright by G. Gaigalas and Ch. F. Fischer"
      print *, "              (Fortran 95 version)               NIST  (2017)."
      print *, "              Input files: rcsfmr.inp, rcsf.inp"
      print *, "              Output file: rcsf.out"
      print *, ""
      print *, 'Reduction based on Dirac-Coulomb (1) or'
      print *, 'Dirac-Coulomb-Breit (2) Hamiltonian?'
      READ(*,*) ICOLBREI
!
      NBLOCK = 0 
      CALL FACTT 
      CALL SET_CSF_list (NCORE,NPEEL)
      WRITE (6, *) " Block    MR NCSF   Before NCSF   After NCSF"
      DO
         I_Count = 0
         CALL LODCSL_MR (NCORE,NPEEL,NCFD,NEXT_BLOCK)
         CSF_Number = 0
         DO 
            CSF_Number = CSF_Number + 1
            CALL LODCSL_CSF (NCFD,CSF_Number,NCORE,NPEEL,NEXT_CSF)
            IF(.NOT. NEXT_CSF) EXIT
            IF(Found_CSF == 1) CYCLE
            CALL Interact_MR (ICOLBREI,I_Count)
         END DO
         WRITE (6,'(3X,I2,6X,I7,3X,I10,3X,I10)')                       &
               NBLOCK,NCFBLK(NBLOCK),CSF_Number-1,I_count+NCFBLK(NBLOCK)
!
         deallocate (Found)
         deallocate (C_shell)
         deallocate (C_quant)
         deallocate (C_coupl)
         deallocate (iqa) 
         deallocate (jqsa) 
         deallocate (jcupa)
!
         IF(.NOT. NEXT_BLOCK) EXIT
         WRITE(22,'(A2)') ' *' 
      END DO
      CLOSE(24) 
      call stoptime (ncount1, 'RCSFinteract')
      STOP  
      END PROGRAM RCSFinteract
