!***********************************************************************
!                                                                      *
      SUBROUTINE ENDSUM 
!                                                                      *
!   Generates the last part of  rscf92.sum  (on stream 24).            *
!                                                                      *
!   Call(s) to: [LIB92]: ENGOUT, RINT.                                 *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNW*NCF),                   *
!                                  JCUPA(NNNW*NCF)                     *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE def_C
      USE eigv_C
      USE orb_C
      USE scf_C
!gg      USE syma_C
      USE wave_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rint_I 
      USE engoutgg_I 
      USE csfwgt_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, MODE
      REAL(DOUBLE) :: WA, WB, WC, WD, WE
!-----------------------------------------------
!
!   Write out the orbital properties
!
      WRITE (24, 301) 
      DO I = 1, NW 
         WRITE (24, 302) NP(I),NH(I),E(I),PZ(I),GAMA(I),PF(2,I),QF(2,I), &
            SCNSTY(I), MF(I) 
      END DO 
!
      WRITE (24, 303) 
      DO I = 1, NW 
         WA = RINT(I,I,-1) 
         WB = RINT(I,I,1) 
         WC = RINT(I,I,2) 
         WD = RINT (I,I, 4)
         WE = 0.d0
         IF (NH(I) /= 's ' .AND. NH(i) /= 'p-') then
            WE = RINT(I,I,-3)
         END IF
         WRITE (24,304) NP(I),NH(I),WE,WA,WB,WC,WD, UCF(I)
      END DO 
!
      IF (NCMIN /= 0) THEN 
         MODE = 0 
         CALL ENGOUTGG (EVAL,ICCMIN,NCMIN,MODE)
!GG         CALL ENGOUT (EVAL, IATJPO, IASPAR, ICCMIN, NCMIN, MODE) 
         CALL CSFWGT (.FALSE.) 
      ENDIF 
!
      CLOSE(24) 
!
      RETURN  
!
  301 FORMAT(/,'Radial wavefunction summary:'/,/,67X,'Self'/,'Subshell',6X,'e',&
         13X,'p0',5X,'gamma',5X,'P(2)',7X,'Q(2)',3X,'Consistency',' MTP'/) 
  302 FORMAT(1X,I2,A2,1X,1P,D17.10,1P,D11.3,0P,F6.2,1P,3(D11.3),I5) 
  303 FORMAT (/18X,'-3',14X,'-1',29X,'2',14x,'4',5X,'Generalised'     &
              /'Subshell',4X,'<  r  >',8X,'<  r  >',8X,'<  r  >',8X,  &
               '<  r  >',8X,'<  r  >',6X,'occupation'/)
  304 FORMAT (1X,1I2,1A2,1X,1P,6D15.5)
      RETURN  
!
      END SUBROUTINE ENDSUM 
