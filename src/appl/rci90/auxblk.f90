!************************************************************************     
      SUBROUTINE AUXBLK(J2MAX, ATWINV) 
!************************************************************************     
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE parameter_def,   ONLY: NNNW
      USE decide_C
      USE def_C
      USE grid_C
      USE ncdist_C
      USE tatb_C
      USE coeils_C
      USE vpilst_C,indoei=>indvpi, valoei=>valvpi
      USE bilst_C
      USE keilst_C
      USE vinlst_C
      USE orb_C, ONLY: ncf, nw, iqa
      USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ncharg_I 
      USE vacpol_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(IN) :: J2MAX 
      REAL(DOUBLE), INTENT(OUT) :: ATWINV 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, K 
      LOGICAL :: LDBPG 
!-----------------------------------------------
!
      FRSTCO = .TRUE. 
      NCOEI = 0 
 
      IF (LTRANS) THEN 
!        ...Check the maximum numbers of orbtitals allowed in brint.f
         I = NNNW 
         SELECT CASE (J2MAX)  
         CASE (11)  
            I = 114 
         CASE (12)  
            I = 112 
         CASE (13)  
            I = 110 
         CASE (14)  
            I = 108 
         CASE (15)  
            I = 106 
         CASE (16)  
            I = 105 
         CASE (17)  
            I = 103 
         CASE (18)  
            I = 101 
         CASE (19)  
            I = 100 
         CASE (21:)  
            I = 90 
         END SELECT 
 
         IF (I < NW) THEN 
            WRITE (ISTDE, *) 'In setham. The number of orbitals is too' 
            WRITE (ISTDE, *) 'large for the brint routine' 
            STOP  
         ENDIF 
 
         FIRST = .TRUE. 
         NTPI = 0 
      ENDIF 
!
!   Initialisations for the vacuum polarisation corrections
!
      IF (LVP) THEN 
         CALL NCHARG 
         CALL VACPOL 
         ZDIST(2:N) = TB(2:N)*RP(2:N) 
         FRSTVP = .TRUE. 
         NVPI = 0 
      ENDIF 
!
!   Initialisations for nuclear translational energy corrections
!
      IF (EMN > 0.D0) THEN 
         ATWINV = 1.D0/EMN 
         IF (LNMS) THEN 
            FRSTKI = .TRUE. 
            NKEI = 0 
         ENDIF 
         IF (LSMS) THEN 
            FRSTVI = .TRUE. 
            NVINTI = 0 
         ENDIF 
      ELSE 
         ! atwinv will not be used
         LNMS = .FALSE. 
         LSMS = .FALSE. 
      ENDIF 
 
      RETURN  
      END SUBROUTINE AUXBLK 
