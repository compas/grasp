!***********************************************************************
!                                                                      *
      SUBROUTINE MCTIN(IOPAR, JKP, NAME)
!                                                                      *
!   This routine loads  coefficients with parity and  rank specified   *
!   by KP(JKP) into the arrays ISLDR and XSLDR.  IOPAR is the parity   *
!   (+/- 1) and is determined from the sign of  KP(JKP).               *
!                                                                      *
!                                         Last revision: 28 Dec 1992   *
!   Updated by Jacek Bieron               Last revision: 10 Mar 1994   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:26:50   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE parameter_def,   ONLY: KEYORB, NNNW
      USE debug_C
      USE decide_C
      USE def_C
      USE foparm_C
      USE OFFD_C
      USE orb_C
      USE OSC_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE alcnsa_I
      USE alcnta_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IOPAR
      INTEGER, INTENT(IN) :: JKP
      CHARACTER  :: NAME(2)*24
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: NFILE = 93
      INTEGER, PARAMETER :: NFILE1 = 237
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NFILE2, M, K, IBLKI, IBLKF, I, LABL, NCSF, J
      REAL(DOUBLE), DIMENSION(NNNW) :: TSHELL
      LOGICAL :: AVAIL
!-----------------------------------------------
!
!   Read the data back as required by OSCL conventions
!
      NFILE2 = NFILE1 + JKP

      M = 0
      K = 0
!
      READ (NFILE2) IBLKI, IBLKF, NW, NKP
      READ (NFILE2) NINT
!
      DO I = 1, NINT
!
         READ (NFILE2) LABL, NCSF
!
         M = M + 1
!bieron
         IF (M >= NSDIM) CALL ALCNSA (JJA, JJB, HB1, HB2, HC1, &
            HC2, HM1, HM2, LAB, NPTR, NSDIM, 2)
         LAB(M) = LABL
!
!   Read configuration pairs and coefficients for this integral
!
    4    CONTINUE
         IF (NCSF + K > NTDIM) THEN
            CALL ALCNTA (ISLDR, ISLDR1, XSLDR, NTDIM, 2)
            GO TO 4
         ENDIF
         NPTR(M) = K
         READ (NFILE2) (ISLDR(J + K),ISLDR1(J + K),XSLDR(J + K),J=1,NCSF)
!          write(*,*) (ISLDR(J+K),XSLDR(J+K),J = 1,NCSF)
         K = K + NCSF
!
      END DO
!
!   Close (and hence release) the scratch file
!
!     IF (AVAIL) THEN
!
!        CLOSE (unit=NFILE,status="DELETE")
!      ENDIF
!
      NPTR(M+1) = K
      NINTEG = M
!
      RETURN
!
  301 FORMAT(/,/,/,1X,I8,' MCT coefficients generated for rank ',I2,&
         ' and parity ',I2,/,/)
      RETURN
!
      END SUBROUTINE MCTIN
