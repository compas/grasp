!
!***********************************************************************
!                                                                      *
      MODULE stat_C
!
!***********************************************************************
!                                                                      *
      USE vast_kind_param, ONLY:  DOUBLE, BYTE
!...Created by Pacific-Sierra Research 77to90  4.3E  07:21:55   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      INTEGER, DIMENSION(:,:,:), pointer :: jqsar
      INTEGER(BYTE), DIMENSION(:,:,:), pointer :: jqsa
      INTEGER, DIMENSION(:,:), pointer :: jcupar
      INTEGER(BYTE), DIMENSION(:,:), pointer :: jcupa
      END MODULE stat_C
