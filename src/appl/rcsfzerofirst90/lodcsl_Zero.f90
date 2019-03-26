!***********************************************************************
!                                                                      *
      SUBROUTINE LODCSL_Zero(NEXT_BLOCK)
!                                                                      *
!   Loads the data from the  .csl  file. A number of checks are made   *
!   to ensure correctness and consistency.                             *
!                                                                      *
!                                                                      *
!   Written by  G. Gaigalas                        Vilnius, May 2016   *
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
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL,  INTENT(OUT) :: NEXT_BLOCK
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER            :: IOS, NCF, NCFD
      CHARACTER(LEN=256) :: RECORD
!-----------------------------------------------
!
!   Initial allocation for arrays with a dimension dependent
!   on the number of CSFs; the initial allocation must be
!   greater than 1
!
      NEXT_BLOCK = .TRUE.
      NCFD = NUM_in_BLK(NBLOCK+1)
      allocate (Found(NCFD))
      allocate (C_shell(NCFD))
      allocate (C_quant(NCFD))
      allocate (C_coupl(NCFD))
      NCF = 0
      DO
         NCF = NCF + 1
         READ (21, '(A)', IOSTAT=IOS) RECORD
         IF (IOS == 0) THEN
            IF (RECORD(1:2) == ' *') THEN
               NBLOCK = NBLOCK + 1
               NCFBLK(NBLOCK) = NCF - 1
               NotFound = NCFBLK(NBLOCK)
               Found(1:NotFound) = 0
               RETURN
            ENDIF
            C_shell(NCF) = RECORD
!
!   Read the J_sub and v quantum numbers
!
            READ (21, '(A)', IOSTAT=IOS) RECORD
            C_quant(NCF) = RECORD
!
!   Read the X, J, and (sign of) P quantum numbers
!
            READ (21, '(A)') RECORD
            C_coupl(NCF) = RECORD
            WRITE(22,'(A)') TRIM(C_shell(NCF))
            WRITE(22,'(A)') TRIM(C_quant(NCF))
            WRITE(22,'(A)') TRIM(C_coupl(NCF))
         ELSE
            EXIT
         ENDIF
      END DO
      NBLOCK = NBLOCK + 1
      NCFBLK(NBLOCK) = NCF - 1
      NotFound = NCFBLK(NBLOCK)
      Found(1:NotFound) = 0
      NEXT_BLOCK = .FALSE.
      RETURN
!
      END SUBROUTINE LODCSL_Zero
