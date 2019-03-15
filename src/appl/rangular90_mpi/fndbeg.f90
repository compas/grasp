!***********************************************************************
!                                                                      *
      SUBROUTINE FNDBEG(JASTRT, JBSTRT, INDEX, LLISTT, LLISTV)
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:01:42   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
!-----------------------------------------------
!   C o m m o n   B l o c k s
!-----------------------------------------------
      USE mcp_C
      USE orb_C
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(INOUT) :: JASTRT
      INTEGER, INTENT(INOUT) :: JBSTRT
      INTEGER, INTENT(INOUT) :: INDEX
      INTEGER, INTENT(OUT) :: LLISTT
      INTEGER, DIMENSION(:), pointer :: LLISTV
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IOS, ICREAD, IRREAD, INDX, NREC, I, LABEL
      REAL(DOUBLE) :: COEFF
      CHARACTER :: SRTLAB*8, MCPLAB*3
!-----------------------------------------------
!
!   Begin by examining file 30; this is the last to be updated by
!   SUBROUTINE MCP
!
!   Read and check the character part of the file header
!
      REWIND (30)
      READ (30) MCPLAB, SRTLAB
      IF (SRTLAB == '  SORTED') THEN
         JASTRT = NCF + 1
         JBSTRT = NCF + 1
         LLISTT = 0
         LLISTV(:KMAX) = 0
         GO TO 8
      ENDIF
!
      READ (30)
      READ (30)
!
!   Read as many records as possible
!
    2 CONTINUE
      READ (30, IOSTAT=IOS) ICREAD, IRREAD, INDX
!
      IF (IOS == 0) THEN
!
!   No errors or end-of-file; keep reading
!
         JASTRT = ICREAD
         JBSTRT = IRREAD
         INDEX = INDX
         GO TO 2
!
      ELSE
!
         IF (JASTRT==NCF .AND. JBSTRT==NCF) THEN
!
!   All coefficients have been generated; sorting may still
!   be necessary; force this option
!
            JASTRT = NCF + 1
            JBSTRT = NCF + 1
!
         ELSE
!
!   Some coefficients remain to be generated; reposition all files
!   for augmentation of lists by SUBROUTINE MCP; update JBSTRT and,
!   if appropriate, JASTRT
!
            DO K = 31, 32 + KMAX
               REWIND (K)
               NREC = 3
               DO I = 1, NREC
                  READ (K)
               END DO
    4          CONTINUE
               READ (K, IOSTAT=IOS) INDX, LABEL, COEFF
               IF (IOS==0 .AND. INDX<=INDEX) THEN
                  NREC = NREC + 1
                  GO TO 4
               ELSE
                  REWIND (K)
                  DO I = 1, NREC
                     READ (K)
                  END DO
                  IF (K > 31) THEN
                     LLISTV(K-32) = NREC - 3
                  ELSE
                     LLISTT = NREC - 3
                  ENDIF
               ENDIF
            END DO
!
!  Now, reposition the sms file. This file should contain the
!  same number of data records as file 33.

            REWIND (20)
            DO I = 1, LLISTV(1)
               READ (20)
            END DO

            JBSTRT = JBSTRT + 1
            IF (JBSTRT > NCF) THEN
               JASTRT = JASTRT + 1
               JBSTRT = JASTRT
            ENDIF
!
         ENDIF
!
      ENDIF
!
    8 CONTINUE
      RETURN
      END SUBROUTINE FNDBEG
