!***********************************************************************
!                                                                      *
      SUBROUTINE SMSREAD(VINT,VINT2)
!
!   Call(s) to: [LIB92]: ALCBUF                                        *
!                                                                      *
!   WRITTEN  by C. Naz\'e  Oct. 2011                                   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE parameter_def,    ONLY: KEYORB, NNNW
      USE prnt_C
      USE ris_C
      USE orb_C
      USE eigv_C
      USE BUFFER_C
      USE debug_C,          ONLY: CUTOFF
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE alcbuf_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE), DIMENSION(NNNW,NNNW), INTENT(IN) :: VINT, VINT2
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: KEY = KEYORB
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(DOUBLE) :: CONTRI, CONTRIK1, COEFFSMS
      INTEGER :: IC, IR, I, J, IIA, IIB, IIC, IID, LOC, LAB, IOS
!-----------------------------------------------
!
      CALL ALCBUF (1)

      REWIND(51)
  16  READ (51,IOSTAT = IOS) IC,IR
!C      WRITE(*,*) IC,IR
      IF (IOS .EQ. 0) THEN
!cc        DO J = 1,NVEC
          READ(51) COEFFSMS,LAB
          IID = MOD (LAB, KEY)
          LAB = LAB/KEY
          IIB  = MOD (LAB, KEY)
          LAB = LAB/KEY
          IIC = MOD (LAB, KEY)
          IIA = LAB/KEY
        DO J = 1,NVEC
          LOC = (J-1)*NCF
            CONTRIK1 = - EVEC(IC+LOC)*EVEC(IR+LOC)                     &
               * COEFFSMS                                              &
               * VINT (IIA,IIC)*VINT(IIB,IID)
            CONTRI = - EVEC(IC+LOC)*EVEC(IR+LOC)                       &
               * COEFFSMS                                              &
               * ( VINT2(IIA,IIC)*VINT(IIB,IID)                        &
               + VINT2(IIB,IID)*VINT(IIA,IIC))/2.0D00
          IF (IR.NE.IC) THEN
            CONTRI = 2.0D00 * CONTRI
            CONTRIK1 = 2.0D00 * CONTRIK1
          ENDIF
          SMSC1(J) = SMSC1(J) + CONTRIK1
          SMSC2(J) = SMSC2(J) + CONTRI
        ENDDO

      GOTO 16
      ENDIF
      CALL ALCBUF (3)
      RETURN
      END SUBROUTINE SMSREAD
