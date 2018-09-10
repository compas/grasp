!***********************************************************************
!                                                                      *
      SUBROUTINE Interact_CSF(JA,JB,ICOLBREI,int_CSF)
!                                                                      *
!                                                                      *
!                                                                      *
!                                                                      *
!   Written by  G. Gaigalas                   NIST, December 2015      *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE CONS_C
      USE m_C
      USE orb_C
      USE dumx_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
!      USE cord_I
      USE itjpo_I
      USE ispar_I
      USE ichkq2_I
      USE setqna_I
      USE onescalarINT_I
      USE el1INT_I
      USE el2INT_I
      USE el3INT_I
      USE el4INT_I
      USE el5INT_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
!      EXTERNAL CORD
      INTEGER, INTENT(IN) :: JA,JB,ICOLBREI
      INTEGER, INTENT(OUT) :: int_CSF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I,IB1,IB2,IDQ,IDQG,II,I1,IJ,IJW,IM,IPCA,IT1,IT2, &
                 JA1,JB1,JA2,JB2,J,JW,JT1,JT2,JT3,NDQ, &
                 KLAST,KW,KWA,K1,KW1,KW2,NPEELM, IA, IB
      INTEGER :: INTERACT
!-----------------------------------------------
!
!   The Hamiltonian is an even scalar operator
      int_CSF = 0
      IF (ITJPO (JA) .NE. ITJPO (JB)) RETURN
      IF (ISPAR (JA) .NE. ISPAR (JB)) RETURN
      IF (JA .EQ. JB)  THEN
         PRINT*, "Interact_CSF: The JA and JB is the same.",JA,JB
         STOP
      END IF
      IF (ICHKQ2(JA,JB) .EQ. 0) RETURN
      CALL SETQNA (JA,JB)
!
!   1.0 Analyse peel shell interactions
!
!   1.1 Analyse electron distribution in peel. (The full procedure is
!       needed only if the number of peel orbitals NPEEL .GE. 2)
      IF (NW .LT. 1) THEN
         int_CSF = 1
         PRINT *, 'Interact_CSF: No subshellsi NW.'
         STOP
      ENDIF
      IF (NPEEL .EQ. 0) THEN
         RETURN
         PRINT *, 'Interact_CSF: No subshells NPEEL.'
         STOP
      END IF
      CALL ONESCALARINT(JA,JB,IA,IB,INTERACT)
!      print*,"ONESCALARINT",INTERACT
      IF(INTERACT .NE. 0) THEN
         int_CSF = 1
         RETURN
      END IF
      IF (NPEEL .EQ. 1) THEN
!
!   1.2.4 IDQ .EQ. 0 - diagonal case. Include all pairs with
!         JA1 = JA2, JB1 = JB2.
         DO KW1 = 1,NPEEL
            K1 = JLIST(KW1)
            JB1 = KW1
            JB2 = KW1
            DO KW2 = 1,KW1
               JA1 = KW2
               IF (JA1 .EQ. JB1) THEN
                 IF (NQ1(K1) .LE. 1) CYCLE
               END IF
               JA2 = JA1
               IF(JA.NE.JB) THEN
                 IF(IDQG.NE.2) THEN
                   WRITE(99,996)
                   WRITE(99,*)"JA,JB,JA1,JB1",JA,JB,JA1,JB1
                   WRITE(99,*)"IDQG IDQ",IDQG,IDQ
                   RETURN
                 ELSE
!
!                   TARP KONFIGURACIJU
!                           KAI      N'1=N1+1
!                                    N'2=N2-1
!                                    N'3=N3
                    CALL EL3INT(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI,INTERACT)
!                    print*, "CIA EL3INT"
                    IF(INTERACT .NE. 0) THEN
                       int_CSF = 1
                       RETURN
                    END IF 
                 END IF
               ELSE
!
!                TARP TU PACIU BUSENU
                 CALL EL1INT(JA,JB,JA1,JB1,0,ICOLBREI,INTERACT)
!                 print*, "EL1INT"
                 IF(INTERACT .NE. 0) THEN
                    int_CSF = 1
                    RETURN
                 END IF 
               END IF
            END DO
         END DO
         RETURN
      END IF
!
!   Find differences in occupations, NDQ, for each peel orbital in
!   turn and use to set up labels of active orbitals maintaining the
!   convention JA1 .LE. JB1, JA2 .LE. JB2.
      IDQ = 0
      JA1 = 0
      JB1 = 0
      JA2 = 0
      JB2 = 0
      IDQG = 0
      DO JW = 1,NPEEL
         J = JLIST(JW)
         NDQ = NQ1(J) - NQ2(J)
         IF (IABS (NDQ) .GT. 2) RETURN
         IF (NDQ .LT. 0) THEN
           IF (NDQ+1 .GT. 0) THEN
             CYCLE
           ELSE IF (NDQ+1 .EQ. 0) THEN
             IF (JA2 .GT. 0) THEN
               JB2 = JW
             ELSE
               JA2 = JW
             END IF
             IDQ = IDQ+1
             IDQG=1+IDQG
           ELSE IF (NDQ+1 .LT. 0) THEN
             JA2 = JW
             IDQ = IDQ+2
             IDQG=20+IDQG
           END IF
         ELSE
           IF (NDQ-1 .GT. 0) THEN
             JA1 = JW
             IDQ = IDQ+2
             IDQG=20+IDQG
             CYCLE
           ELSE IF (NDQ-1 .EQ. 0) THEN
             IF (JA1 .GT. 0) THEN
               JB1 = JW
             ELSE
               JA1 = JW
             END IF
             IDQ = IDQ+1
             IDQG=1+IDQG
             CYCLE
           ELSE
             CYCLE
           END IF
         END IF
      END DO
!
!   1.2 Calculate coefficients for all possible sets of active shells.
!
!   There are 4 cases, depending on the value of IDQ, the sum of the
!   absolute differences NDQ:
!
!   1.2.1 IDQ .GT. 4: matrix element null
      IF (IDQ .GT. 4) RETURN
      IF (IDQ .EQ. 4) THEN
!
!   1.2.2 IDQ .EQ. 4: matrix element uniquely defined
        IF (JB1 .EQ. 0) THEN
          JB1 = JA1
        END IF
        IF (JB2 .EQ. 0) THEN
          JB2 = JA2
        END IF
        CONTINUE
        IF(IDQG.NE.40) THEN
          IF(IDQG.NE.22) THEN
!
!         TARP KONFIGURACIJU
!                        KAI      N'1=N1+-1
!                                 N'2=N2+-1
!                        KAI      N'3=N3+-1
!                                 N'4=N4+-1
            CALL EL5INT(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI,INTERACT)
!            print*, "EL5INT"
            IF(INTERACT .NE. 0) THEN
               int_CSF = 1
               RETURN
            END IF
          ELSE
            CALL EL4INT(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI,INTERACT)
!            print*, "EL4INT"
            IF(INTERACT .NE. 0) THEN
               int_CSF = 1
               RETURN
            END IF
          END IF
        ELSE
!
!       TARP KONFIGURACIJU
!                        KAI      N'1=N1+1
!                                 N'2=N2-1
          CALL EL2INT(JA,JB,JA1,JA2,ICOLBREI,INTERACT)
!          print*, "EL2INT",JA,JB
          IF(INTERACT .NE. 0) THEN
             int_CSF = 1
             RETURN
          END IF
        END IF
        RETURN
      END IF
!
      KLAST = NPEEL
      IF (IDQ .NE. 2) THEN
        IF (IDQ .NE. 0) THEN
!
!         3.0 Diagnostic print - NW .LT. 1
          WRITE (*,300)
          STOP
        END IF
      ELSE
        KLAST = 1
        GOTO 16
      END IF
!
!
!   1.2.3 IDQ .EQ. 2: One orbital fixed each side include all
!                     possible spectators.
!
!   Also IDQ .EQ. 0 for a matrix element off-diagonal in coupling
!   only. Must sum over all pairs of orbitals excluding core-core
!   terms
   16 DO KWA = 1,KLAST
         IF (IDQ .NE. 2) THEN
           JA1 = KWA
           JA2 = KWA
         END IF
         JT1 = JA1
         JT2 = JA2
         IT1 = JLIST(JA1)
         IT2 = JLIST(JA2)
         DO KW = KWA,NPEEL
            K1 = JLIST(KW)
            IF (NQ1(K1)*NQ2(K1) .EQ. 0) CYCLE
            JB1 = KW
            JB2 = KW
            JA1 = JT1
            JA2 = JT2
!
!   Interchange JA1 and JB1 and/or JA2 and JB2 if necessary
            IF (JA1-JB1 .GT. 0) THEN
              JT3 = JB1
              JB1 = JA1
              JA1 = JT3
            ELSE IF (JA1-JB1 .EQ. 0) THEN
              IB1 = JLIST(JB1)
              IF (NQ1(IB1) .LE. 1) CYCLE
            END IF
            IF (JA2-JB2 .GT. 0) THEN
              JT3 = JB2
              JB2 = JA2
              JA2 = JT3
            ELSE IF (JA2-JB2 .EQ. 0) THEN
              IB2 = JLIST(JB2)
              IF (NQ2(IB2) .LE. 1) CYCLE
            END IF
            IF(IDQ.NE.0) THEN
!
!     TARP KONFIGURACIJU
!                        KAI      N'1=N1+1
!                                 N'2=N2-1
!                                 N'3=N3
              CALL EL3INT(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI,INTERACT)
!              print*, "EL3INT 2"
              IF(INTERACT .NE. 0) THEN
                int_CSF = 1
                RETURN
              END IF 
            ELSE
!
!     TARP TU PACIU KONFIGURACIJU
              CALL EL1INT(JA,JB,JA1,JB1,1,ICOLBREI,INTERACT)
!            print*, "EL1INT   2"
              IF(INTERACT .NE. 0) THEN
                int_CSF = 1
                RETURN
              END IF 
            END IF
         END DO
         IF ((IDQ .EQ. 0) .AND. (NCORE .EQ. 0)) CYCLE
         IF ((NCORE .EQ. 0) .OR. (NAK(IT1) .NE. NAK(IT2))) RETURN
!
!   This section calculates the terms arising from active electrons
!   which are in closed shells
         NPEELM = NPEEL-1
         DO I = 1,NPEEL
            JLIS(I) = JLIST(I)
         END DO
         DO I = 1,NPEELM
            JC1S(I) = JJC1(I)
            JC2S(I) = JJC2(I)
         END DO
         DO KW = 1,NCORE
            IJW = KLIST(KW)
            DO I = 1,NPEEL
               IJ = JLIST(I)
               IF (IJW .LT. IJ) GOTO 29
            END DO
            I = NPEEL+1
            GOTO 31
   29       IM = NPEEL-I+1
            DO II = 1,IM
               JLIST(NPEEL+2-II) = JLIST(NPEEL+1-II)
               IF (NPEEL .EQ. II) GOTO 31
               JJC1(NPEEL+1-II) = JJC1(NPEEL-II)
               JJC2(NPEEL+1-II) = JJC2(NPEEL-II)
            END DO
   31       CONTINUE
            IF (I .LT. 3) THEN
              I1 = JLIST(1)
              JJC1(1) = JJQ1(3,I1)
              JJC2(1) = JJQ2(3,I1)
            ELSE
              JJC1(I-1) = JJC1(I-2)
              JJC2(I-1) = JJC2(I-2)
            END IF
            JLIST(I) = IJW
            JA1 = JT1
            IF (JT1 .GE. I) JA1 = JA1+1
            JB1 = I
            JA2 = JT2
            IF (JT2 .GE. I) JA2 = JA2+1
            JB2 = I
            IF (JA1-JB1 .GT. 0) THEN
              JT3 = JB1
              JB1 = JA1
              JA1 = JT3
            END IF
            IF (JA2-JB2 .GT. 0) THEN
              JT3 = JB2
              JB2 = JA2
              JA2 = JT3
            END IF
            NPEEL = NPEEL+1
            IF(IDQ.NE.0) THEN
              IF(IDQG.NE.40) THEN
                IF(IDQG.NE.2) THEN
                   WRITE(99,995)
                   RETURN
                ELSE
!
!     TARP KONFIGURACIJU
!                        KAI      N'1=N1+1
!                                 N'2=N2-1
!                                 N'3=N3
                  CALL EL3INT(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI,INTERACT)
!                  print*, "EL3INT 3"
                  IF(INTERACT .NE. 0) THEN
                    int_CSF = 1
                    RETURN
                  END IF 
                END IF
              ELSE
                WRITE(99,994)
                RETURN
              END IF
            END IF
!
!     TARP TU PACIU KONFIGURACIJU
            CALL EL1INT(JA,JB,JA1,JB1,1,ICOLBREI,INTERACT)
!            print*, "EL1INT   3"
            IF(INTERACT .NE. 0) THEN
              int_CSF = 1
              RETURN
            END IF 
            NPEEL = NPEEL-1
            NPEELM = NPEEL-1
            DO I = 1,NPEEL
               JLIST(I) = JLIS(I)
            END DO
            DO I = 1,NPEELM
               JJC1(I)  = JC1S(I)
               JJC2(I)  = JC2S(I)
            END DO
         END DO
      END DO
      RETURN
  300 FORMAT ('Interact_CSF: Error.')
  994 FORMAT('   rie zymes 38?? atv N=N-N  !!!!!!')
  995 FORMAT('   rie zymes 38 atv N=N-N  !!!!!!')
  996 FORMAT('   rie zymes 45 atv N=N-N  !!!!!!')
      END SUBROUTINE Interact_CSF
