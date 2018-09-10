!***********************************************************************
      SUBROUTINE CSFWGT(LSTDIO) 
!                                                                      *
!   Print  the  weights of the largest five CSF contributors to each   *
!   ASF.                                                               *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC, ISPAR, ITJPO.          *
!                                                                      *
!                                          Last updated: 21 Dec 1992   *
!                                          Last updated: 24 Feb 1997   *
!   Midified by G. Gaigalas                              05 Feb 2017   *
!      It was deleted the arrays:  JQSA(3*NNNWP*NCF),                  *
!                                  JCUPA(NNNWP*NCF)                    *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:06:32   1/ 3/07  
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE 
      USE def_C
      USE eigv_C
      USE jlabl_C, LABJ=>JLBR, LABP=>JLBP
      USE syma_C,          ONLY: JPGG
      USE iounit_C
      USE hblock_C
      USE pos_C
      USE blkidx_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL, INTENT(IN) :: LSTDIO 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER, DIMENSION(5) :: ICONF 
      INTEGER :: IV, JBLOCK, NCF, NCFPAT, NCMINPAT, NEVECPAT, NEVECOFF, NELT, &
         ICF, IVTJPO, IVSPAR, I, J, ITEMP, IREM, IP 
      REAL(DOUBLE) :: TEMP, W
      REAL(DOUBLE), DIMENSION(5) :: COEFF
      CHARACTER :: RECORD*256, CNUM*8 
!-----------------------------------------------
!
 
      IF (LSTDIO) THEN 
         WRITE (ISTDO, 300) 
      ELSE 
         WRITE (24, 300) 
      ENDIF 
 
      DO IV = 1, NCMIN 
!         loop over eigenvectors
 
         JBLOCK = IDXBLK(IV) 
         NCF = NCFBLK(JBLOCK) 
         NCFPAT = NCFPAST(JBLOCK) 
         NCMINPAT = NCMINPAST(JBLOCK) 
         NEVECPAT = NEVECPAST(JBLOCK) 
         NEVECOFF = NEVECPAT + (IV - NCMINPAT - 1)*NCF 
         NELT = MIN(5,NCF)                       !Find maximum 5 ... within block 
         ICF = ICCMIN(IV) 
!GGGG
         ivtjpo = IABS(JPGG(jblock))
         IF(JPGG(jblock) .GE. 0) THEN
            ivspar = 2
         ELSE
            ivspar = 1
         END IF
!GG         IVTJPO = IATJPO(IV)                     ! j-value related 
!GG         IVSPAR = IASPAR(IV)                     ! parity related 
 
 
         DO I = 1, NELT 
            COEFF(I) = EVEC(NEVECOFF + I) 
            ICONF(I) = I 
         END DO 
 
!         sort the first nelt in decreasing order
         DO I = 1, NELT 
            DO J = I + 1, NELT 
               IF (ABS(COEFF(J)) <= ABS(COEFF(I))) CYCLE  
               TEMP = COEFF(I) 
               COEFF(I) = COEFF(J) 
               COEFF(J) = TEMP 
               ITEMP = ICONF(I) 
               ICONF(I) = ICONF(J) 
               ICONF(J) = ITEMP 
            END DO 
         END DO 
 
         L20: DO I = NELT + 1, NCF 
            W = EVEC(NEVECOFF + I) 
            IF (W==0.D0 .OR. ABS(W)<=ABS(COEFF(NELT))) CYCLE  L20 
!               we have a non-zero value larger than the largest so far
            DO J = 1, NELT 
               IF (ABS(W) <= ABS(COEFF(J))) CYCLE  
               COEFF(NELT:1+J:(-1)) = COEFF(NELT-1:J:(-1)) 
               ICONF(NELT:1+J:(-1)) = ICONF(NELT-1:J:(-1)) 
               COEFF(J) = W 
               ICONF(J) = I 
               CYCLE  L20 
            END DO 
         END DO L20 
 
!GG         IP = (IASPAR(IV) + 3)/2 
           ip =  ivspar
 
         IF (LSTDIO) THEN 
            WRITE (ISTDO, 320) JBLOCK, ICF, LABJ(IVTJPO), LABP(IP), (COEFF(I),I&
               =1,NELT) 
            WRITE (ISTDO, 330) (ICONF(I),I=1,NELT) 
         ELSE 
 
            WRITE (24, 320) JBLOCK, ICF, LABJ(IVTJPO), LABP(IP), (COEFF(I),I=1,&
               NELT) 
            WRITE (24, 330) (ICONF(I),I=1,NELT) 
         ENDIF 
      END DO 
 
  300 FORMAT(/,'Weights of major contributors to ASF:'/,/,&
         'Block Level J Parity      CSF contributions'/) 
  310 FORMAT(1X,A14,80A) 
  320 FORMAT(I3,1X,I5,2X,2A4,5(3X,F8.4)) 
  330 FORMAT(19X,5(3X,I8)) 
 
      RETURN  
      END SUBROUTINE CSFWGT 
