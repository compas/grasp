!*******************************************************************
!                                                                  *
      SUBROUTINE setmcp_mem(NBLOCK,NFILE)
!                                                                  *
!   Written by  G. Gaigalas               Vilnius, September 2021  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE, LONG
      USE memory_man
      USE rmcdhf_mem_C
      USE iounit_C
      USE mpi_C
      USE hmat_C,          ONLY: NELMNT
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use systemmem_I
      use systemfreemem_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: NBLOCK, NFILE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      CHARACTER :: MCPLAB*3
      REAL(DOUBLE) :: COEFF_dumm
      REAL(DOUBLE) :: SYSMemT, SYSMemF
      INTEGER :: I, JBLOCK, JBLOCKT, NCFT, NCOEFF, IOS, NCF
      INTEGER(LONG) :: NCONTR_tot, NCONTR_tot2
      INTEGER :: ICLMN_dumm, LAB_dumm, NCONTR_dumm, INDX_dumm
      INTEGER :: NELMNTGG, IENDCDUM, IROW_dumm
!-----------------------------------------------

      NCONTR_tot = 0
      NCONTR_tot2 = 0
      DO JBLOCK = 1, NBLOCK
        IF(NFILE == 30) THEN
          !*** Read in IROW from file mcp.30 ***
          READ (30) MCPLAB, JBLOCKT, NCF
          IF (JBLOCKT /= JBLOCK) THEN
            WRITE (ISTDE, *) 'setmcp_mem', '1: jblockt .NE. jblock'
            STOP
          ENDIF
          READ (30) NELMNTGG
          NELMNT = INT8(NELMNTGG)
          IF(JBLOCK == 1) THEN
            CALL ALLOC (IMIN1_30,NBLOCK,'IMIN1_30','setmcp_mem')
            CALL ALLOC (IMIN2_30,NBLOCK,'IMIN2_30','setmcp_mem')
            CALL ALLOC (IMAX1_30,NBLOCK,'IMAX1_30','setmcp_mem')
            CALL ALLOC (IMAX2_30,NBLOCK,'IMAX2_30','setmcp_mem')
            IMIN1_30(JBLOCK) = 0
            IMAX1_30(JBLOCK) = NCF
            IMIN2_30(JBLOCK) = 0
            IMAX2_30(JBLOCK) = NELMNT
          ELSE
            IMIN1_30(JBLOCK) = IMAX1_30(JBLOCK-1) + 0
            IMAX1_30(JBLOCK) = IMAX1_30(JBLOCK-1) + NCF
            IMIN2_30(JBLOCK) = IMAX2_30(JBLOCK-1) + 0
            IMAX2_30(JBLOCK) = IMAX2_30(JBLOCK-1) + NELMNT
          END IF
          READ (30) (IENDCDUM,I=MYID + 1,NCF,NPROCS),                  &
                                                 (IROW_dumm,I=1,NELMNT)
        ELSE
          READ (NFILE) MCPLAB, JBLOCKT, NCFT, NCOEFF
          IF (JBLOCKT /= JBLOCK) THEN
            WRITE (ISTDE, *) 'setmcp_mem', ':2 jblockt .NE. jblock'
            STOP
          ENDIF

          NCONTR_tot2 = NCONTR_tot2 + 1
          READ (NFILE, IOSTAT=IOS) LAB_dumm, NCONTR_dumm
          DO WHILE (LAB_dumm /= 0 .OR. NCONTR_dumm /= 0)
            !*** Read column index, sparse matrix index, and
            !*** coefficient for all contributions from this integral
            READ (NFILE) (ICLMN_dumm,INDX_dumm,COEFF_dumm,             &
                                                       I=1,NCONTR_dumm)
            NCONTR_tot = NCONTR_tot + NCONTR_dumm
            NCONTR_tot2 = NCONTR_tot2 + 1
            READ (NFILE, IOSTAT=IOS) LAB_dumm, NCONTR_dumm
          END DO
        END IF
      END DO

      if(NFILE == 30) then
         ! Total physical memory
         SYSMemT = SYSTEMMem()  ! GB
         SYSMemF = SYSTEMFreeMem()  ! GB
         IF(MYID == 0) WRITE(*,'(A53)')                                &
         '-----------------------------------------------------'
         IF(MYID == 0) WRITE(*,'(A53)')                                &
         'Spin-angular coefficient are putting into the memory:'
         IF(MYID == 0) WRITE(*,'(A53)')                                &
         '-----------------------------------------------------'
         IF(MYID == 0) WRITE(*,'(A1)') ''
         IF(MYID == 0) WRITE(*,'(A30,F10.2,A3)')                       &
                          'Total memory on computer:',SYSMemT,' Gb'
         IF(MYID == 0) WRITE(*,'(A30,F10.2,A3)')                       &
                          'Free  memory on computer:',SYSMemF,' Gb'
         IF(MYID == 0) WRITE(*,'(A1)') ''
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.30:'
         CALL ALLOC (IENDC_30,IMAX1_30(NBLOCK),'IENDC_30','setmcp_mem')
         CALL ALLOC (IROW_30, IMAX2_30(NBLOCK),'IROW_30', 'setmcp_mem')
      else if(NFILE == 31) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.31:'
         CALL ALLOC (LAB_31,   NCONTR_tot2,'LAB_31',   'setmcp_mem')
         CALL ALLOC (NCONTR_31,NCONTR_tot2,'NCONTR_31','setmcp_mem')
         CALL ALLOC (ICLMN_31, NCONTR_tot, 'ICLMN_31', 'setmcp_mem')
         CALL ALLOC (INDX_31,  NCONTR_tot, 'INDX_31',  'setmcp_mem')
         CALL ALLOC (COEFF_31, NCONTR_tot, 'COEFF_31', 'setmcp_mem')
      else if(NFILE == 32) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.32:'
         CALL ALLOC (LAB_32,   NCONTR_tot2,'LAB_32',   'setmcp_mem')
         CALL ALLOC (NCONTR_32,NCONTR_tot2,'NCONTR_32','setmcp_mem')
         CALL ALLOC (ICLMN_32, NCONTR_tot, 'ICLMN_32', 'setmcp_mem')
         CALL ALLOC (INDX_32,  NCONTR_tot, 'INDX_32',  'setmcp_mem')
         CALL ALLOC (COEFF_32, NCONTR_tot, 'COEFF_32', 'setmcp_mem')
      else if(NFILE == 33) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.33:'
         CALL ALLOC (LAB_33,   NCONTR_tot2,'LAB_33',   'setmcp_mem')
         CALL ALLOC (NCONTR_33,NCONTR_tot2,'NCONTR_33','setmcp_mem')
         CALL ALLOC (ICLMN_33, NCONTR_tot, 'ICLMN_33', 'setmcp_mem')
         CALL ALLOC (INDX_33,  NCONTR_tot, 'INDX_33',  'setmcp_mem')
         CALL ALLOC (COEFF_33, NCONTR_tot, 'COEFF_33', 'setmcp_mem')
      else if(NFILE == 34) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.34:'
         CALL ALLOC (LAB_34,   NCONTR_tot2,'LAB_34',   'setmcp_mem')
         CALL ALLOC (NCONTR_34,NCONTR_tot2,'NCONTR_34','setmcp_mem')
         CALL ALLOC (ICLMN_34, NCONTR_tot, 'ICLMN_34', 'setmcp_mem')
         CALL ALLOC (INDX_34,  NCONTR_tot, 'INDX_34',  'setmcp_mem')
         CALL ALLOC (COEFF_34, NCONTR_tot, 'COEFF_34', 'setmcp_mem')
      else if(NFILE == 35) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.35:'
         CALL ALLOC (LAB_35,   NCONTR_tot2,'LAB_35',   'setmcp_mem')
         CALL ALLOC (NCONTR_35,NCONTR_tot2,'NCONTR_35','setmcp_mem')
         CALL ALLOC (ICLMN_35, NCONTR_tot, 'ICLMN_35', 'setmcp_mem')
         CALL ALLOC (INDX_35,  NCONTR_tot, 'INDX_35',  'setmcp_mem')
         CALL ALLOC (COEFF_35, NCONTR_tot, 'COEFF_35', 'setmcp_mem')
      else if(NFILE == 36) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.36:'
         CALL ALLOC (LAB_36,   NCONTR_tot2,'LAB_36',   'setmcp_mem')
         CALL ALLOC (NCONTR_36,NCONTR_tot2,'NCONTR_36','setmcp_mem')
         CALL ALLOC (ICLMN_36, NCONTR_tot, 'ICLMN_36', 'setmcp_mem')
         CALL ALLOC (INDX_36,  NCONTR_tot, 'INDX_36',  'setmcp_mem')
         CALL ALLOC (COEFF_36, NCONTR_tot, 'COEFF_36', 'setmcp_mem')
      else if(NFILE == 37) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.37:'
         CALL ALLOC (LAB_37,   NCONTR_tot2,'LAB_37',   'setmcp_mem')
         CALL ALLOC (NCONTR_37,NCONTR_tot2,'NCONTR_37','setmcp_mem')
         CALL ALLOC (ICLMN_37, NCONTR_tot, 'ICLMN_37', 'setmcp_mem')
         CALL ALLOC (INDX_37,  NCONTR_tot, 'INDX_37',  'setmcp_mem')
         CALL ALLOC (COEFF_37, NCONTR_tot, 'COEFF_37', 'setmcp_mem')
      else if(NFILE == 38) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.38:'
         CALL ALLOC (LAB_38,   NCONTR_tot2,'LAB_38',   'setmcp_mem')
         CALL ALLOC (NCONTR_38,NCONTR_tot2,'NCONTR_38','setmcp_mem')
         CALL ALLOC (ICLMN_38, NCONTR_tot, 'ICLMN_38', 'setmcp_mem')
         CALL ALLOC (INDX_38,  NCONTR_tot, 'INDX_38',  'setmcp_mem')
         CALL ALLOC (COEFF_38, NCONTR_tot, 'COEFF_38', 'setmcp_mem')
      else if(NFILE == 39) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.39:'
         CALL ALLOC (LAB_39,   NCONTR_tot2,'LAB_39',   'setmcp_mem')
         CALL ALLOC (NCONTR_39,NCONTR_tot2,'NCONTR_39','setmcp_mem')
         CALL ALLOC (ICLMN_39, NCONTR_tot, 'ICLMN_39', 'setmcp_mem')
         CALL ALLOC (INDX_39,  NCONTR_tot, 'INDX_39',  'setmcp_mem')
         CALL ALLOC (COEFF_39, NCONTR_tot, 'COEFF_39', 'setmcp_mem')
      else if(NFILE == 40) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.40:'
         CALL ALLOC (LAB_40,   NCONTR_tot2,'LAB_40',   'setmcp_mem')
         CALL ALLOC (NCONTR_40,NCONTR_tot2,'NCONTR_40','setmcp_mem')
         CALL ALLOC (ICLMN_40, NCONTR_tot, 'ICLMN_40', 'setmcp_mem')
         CALL ALLOC (INDX_40,  NCONTR_tot, 'INDX_40',  'setmcp_mem')
         CALL ALLOC (COEFF_40, NCONTR_tot, 'COEFF_40', 'setmcp_mem')
      else if(NFILE == 41) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.41:'
         CALL ALLOC (LAB_41,   NCONTR_tot2,'LAB_41',   'setmcp_mem')
         CALL ALLOC (NCONTR_41,NCONTR_tot2,'NCONTR_41','setmcp_mem')
         CALL ALLOC (ICLMN_41, NCONTR_tot, 'ICLMN_41', 'setmcp_mem')
         CALL ALLOC (INDX_41,  NCONTR_tot, 'INDX_41',  'setmcp_mem')
         CALL ALLOC (COEFF_41, NCONTR_tot, 'COEFF_41', 'setmcp_mem')
      else if(NFILE == 42) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.42:'
         CALL ALLOC (LAB_42,   NCONTR_tot2,'LAB_42',   'setmcp_mem')
         CALL ALLOC (NCONTR_42,NCONTR_tot2,'NCONTR_42','setmcp_mem')
         CALL ALLOC (ICLMN_42, NCONTR_tot, 'ICLMN_42', 'setmcp_mem')
         CALL ALLOC (INDX_42,  NCONTR_tot, 'INDX_42',  'setmcp_mem')
         CALL ALLOC (COEFF_42, NCONTR_tot, 'COEFF_42', 'setmcp_mem')
      else if(NFILE == 43) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.43:'
         CALL ALLOC (LAB_43,   NCONTR_tot2,'LAB_43',   'setmcp_mem')
         CALL ALLOC (NCONTR_43,NCONTR_tot2,'NCONTR_43','setmcp_mem')
         CALL ALLOC (ICLMN_43, NCONTR_tot, 'ICLMN_43', 'setmcp_mem')
         CALL ALLOC (INDX_43,  NCONTR_tot, 'INDX_43',  'setmcp_mem')
         CALL ALLOC (COEFF_43, NCONTR_tot, 'COEFF_43', 'setmcp_mem')
      else if(NFILE == 44) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.44:'
         CALL ALLOC (LAB_44,   NCONTR_tot2,'LAB_44',   'setmcp_mem')
         CALL ALLOC (NCONTR_44,NCONTR_tot2,'NCONTR_44','setmcp_mem')
         CALL ALLOC (ICLMN_44, NCONTR_tot, 'ICLMN_44', 'setmcp_mem')
         CALL ALLOC (INDX_44,  NCONTR_tot, 'INDX_44',  'setmcp_mem')
         CALL ALLOC (COEFF_44, NCONTR_tot, 'COEFF_44', 'setmcp_mem')
      else if(NFILE == 45) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.45:'
         CALL ALLOC (LAB_45,   NCONTR_tot2,'LAB_45',   'setmcp_mem')
         CALL ALLOC (NCONTR_45,NCONTR_tot2,'NCONTR_45','setmcp_mem')
         CALL ALLOC (ICLMN_45, NCONTR_tot, 'ICLMN_45', 'setmcp_mem')
         CALL ALLOC (INDX_45,  NCONTR_tot, 'INDX_45',  'setmcp_mem')
         CALL ALLOC (COEFF_45, NCONTR_tot, 'COEFF_45', 'setmcp_mem')
      else if(NFILE == 46) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.46:'
         CALL ALLOC (LAB_46,   NCONTR_tot2,'LAB_46',   'setmcp_mem')
         CALL ALLOC (NCONTR_46,NCONTR_tot2,'NCONTR_46','setmcp_mem')
         CALL ALLOC (ICLMN_46, NCONTR_tot, 'ICLMN_46', 'setmcp_mem')
         CALL ALLOC (INDX_46,  NCONTR_tot, 'INDX_46',  'setmcp_mem')
         CALL ALLOC (COEFF_46, NCONTR_tot, 'COEFF_46', 'setmcp_mem')
      else if(NFILE == 47) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.47:'
         CALL ALLOC (LAB_47,   NCONTR_tot2,'LAB_47',   'setmcp_mem')
         CALL ALLOC (NCONTR_47,NCONTR_tot2,'NCONTR_47','setmcp_mem')
         CALL ALLOC (ICLMN_47, NCONTR_tot, 'ICLMN_47', 'setmcp_mem')
         CALL ALLOC (INDX_47,  NCONTR_tot, 'INDX_47',  'setmcp_mem')
         CALL ALLOC (COEFF_47, NCONTR_tot, 'COEFF_47', 'setmcp_mem')
      else if(NFILE == 48) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.48:'
         CALL ALLOC (LAB_48,   NCONTR_tot2,'LAB_48',   'setmcp_mem')
         CALL ALLOC (NCONTR_48,NCONTR_tot2,'NCONTR_48','setmcp_mem')
         CALL ALLOC (ICLMN_48, NCONTR_tot, 'ICLMN_48', 'setmcp_mem')
         CALL ALLOC (INDX_48,  NCONTR_tot, 'INDX_48',  'setmcp_mem')
         CALL ALLOC (COEFF_48, NCONTR_tot, 'COEFF_48', 'setmcp_mem')
      else if(NFILE == 49) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.49:'
         CALL ALLOC (LAB_49,   NCONTR_tot2,'LAB_49',   'setmcp_mem')
         CALL ALLOC (NCONTR_49,NCONTR_tot2,'NCONTR_49','setmcp_mem')
         CALL ALLOC (ICLMN_49, NCONTR_tot, 'ICLMN_49', 'setmcp_mem')
         CALL ALLOC (INDX_49,  NCONTR_tot, 'INDX_49',  'setmcp_mem')
         CALL ALLOC (COEFF_49, NCONTR_tot, 'COEFF_49', 'setmcp_mem')
      else if(NFILE == 50) then
         IF(MYID == 0) WRITE(*,'(A27)') 'Allocation for mcp.50:'
         CALL ALLOC (LAB_50,   NCONTR_tot2,'LAB_50',   'setmcp_mem')
         CALL ALLOC (NCONTR_50,NCONTR_tot2,'NCONTR_50','setmcp_mem')
         CALL ALLOC (ICLMN_50, NCONTR_tot, 'ICLMN_50', 'setmcp_mem')
         CALL ALLOC (INDX_50,  NCONTR_tot, 'INDX_50',  'setmcp_mem')
         CALL ALLOC (COEFF_50, NCONTR_tot, 'COEFF_50', 'setmcp_mem')
      else
         print*, "Error in setmcp_mem"
         stop
      end if
      REWIND (NFILE)
      IF(NFILE == 30) THEN
         READ (NFILE)
         READ (NFILE)
         READ (NFILE)
      END IF
      READ (NFILE)
      READ (NFILE)
      READ (NFILE)
!
      NCONTR_tot = 0
      NCONTR_tot2 = 0
      DO JBLOCK = 1, NBLOCK
        IF(NFILE == 30) THEN
          !*** Read in IROW from file mcp.30 ***
          READ (30) MCPLAB, JBLOCKT, NCF
          IF (JBLOCKT /= JBLOCK) THEN
            WRITE (ISTDE, *) 'setmcp_mem', '3: jblockt .NE. jblock'
            STOP
          ENDIF
        ELSE
          READ (NFILE) MCPLAB, JBLOCKT, NCFT, NCOEFF
          IF (JBLOCKT /= JBLOCK) THEN
            WRITE (ISTDE, *) 'setmcp_mem', ':4 jblockt .NE. jblock'
            STOP
          ENDIF
          NCONTR_tot2 = NCONTR_tot2 + 1
        ENDIF

         if(NFILE == 30) then
            READ (30) NELMNTGG
            NELMNT = INT8(NELMNTGG)
            READ (30) (IENDC_30(I),I=                                  &
                     IMIN1_30(JBLOCK)+MYID+1,IMAX1_30(JBLOCK),NPROCS), &
                     (IROW_30(I),I=IMIN2_30(JBLOCK)+1,IMAX2_30(JBLOCK))
         else if(NFILE == 31) then
            READ (NFILE, IOSTAT=IOS) LAB_31(NCONTR_tot2),              &
                                  NCONTR_31(NCONTR_tot2)
            DO WHILE (LAB_31(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_31(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_31(NCONTR_tot+I),                   &
                              INDX_31(NCONTR_tot+I),                   &
                             COEFF_31(NCONTR_tot+I),                   &
                        I=1,NCONTR_31(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_31(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_31(NCONTR_tot2),           &
                                     NCONTR_31(NCONTR_tot2)
            END DO
         else if(NFILE == 32) then
            READ (NFILE, IOSTAT=IOS) LAB_32(NCONTR_tot2),              &
                                  NCONTR_32(NCONTR_tot2)
            DO WHILE (LAB_32(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_32(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_32(NCONTR_tot+I),                   &
                              INDX_32(NCONTR_tot+I),                   &
                             COEFF_32(NCONTR_tot+I),                   &
                        I=1,NCONTR_32(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_32(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_32(NCONTR_tot2),           &
                                     NCONTR_32(NCONTR_tot2)
            END DO
         else if(NFILE == 33) then
            READ (NFILE, IOSTAT=IOS) LAB_33(NCONTR_tot2),              &
                                  NCONTR_33(NCONTR_tot2)
            DO WHILE (LAB_33(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_33(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_33(NCONTR_tot+I),                   &
                              INDX_33(NCONTR_tot+I),                   &
                             COEFF_33(NCONTR_tot+I),                   &
                        I=1,NCONTR_33(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_33(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_33(NCONTR_tot2),           &
                                     NCONTR_33(NCONTR_tot2)
            END DO
         else if(NFILE == 34) then
            READ (NFILE, IOSTAT=IOS) LAB_34(NCONTR_tot2),              &
                                  NCONTR_34(NCONTR_tot2)
            DO WHILE (LAB_34(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_34(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_34(NCONTR_tot+I),                   &
                              INDX_34(NCONTR_tot+I),                   &
                             COEFF_34(NCONTR_tot+I),                   &
                        I=1,NCONTR_34(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_34(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_34(NCONTR_tot2),           &
                                     NCONTR_34(NCONTR_tot2)
            END DO
         else if(NFILE == 35) then
            READ (NFILE, IOSTAT=IOS) LAB_35(NCONTR_tot2),              &
                                  NCONTR_35(NCONTR_tot2)
            DO WHILE (LAB_35(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_35(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_35(NCONTR_tot+I),                   &
                              INDX_35(NCONTR_tot+I),                   &
                             COEFF_35(NCONTR_tot+I),                   &
                        I=1,NCONTR_35(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_35(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_35(NCONTR_tot2),           &
                                     NCONTR_35(NCONTR_tot2)
            END DO
         else if(NFILE == 36) then
            READ (NFILE, IOSTAT=IOS) LAB_36(NCONTR_tot2),              &
                                  NCONTR_36(NCONTR_tot2)
            DO WHILE (LAB_36(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_36(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_36(NCONTR_tot+I),                   &
                              INDX_36(NCONTR_tot+I),                   &
                             COEFF_36(NCONTR_tot+I),                   &
                        I=1,NCONTR_36(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_36(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_36(NCONTR_tot2),           &
                                     NCONTR_36(NCONTR_tot2)
            END DO
         else if(NFILE == 37) then
            READ (NFILE, IOSTAT=IOS) LAB_37(NCONTR_tot2),              &
                                  NCONTR_37(NCONTR_tot2)
            DO WHILE (LAB_37(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_37(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_37(NCONTR_tot+I),                   &
                              INDX_37(NCONTR_tot+I),                   &
                             COEFF_37(NCONTR_tot+I),                   &
                        I=1,NCONTR_37(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_37(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_37(NCONTR_tot2),           &
                                     NCONTR_37(NCONTR_tot2)
            END DO
         else if(NFILE == 38) then
            READ (NFILE, IOSTAT=IOS) LAB_38(NCONTR_tot2),              &
                                  NCONTR_38(NCONTR_tot2)
            DO WHILE (LAB_38(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_38(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_38(NCONTR_tot+I),                   &
                              INDX_38(NCONTR_tot+I),                   &
                             COEFF_38(NCONTR_tot+I),                   &
                        I=1,NCONTR_38(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_38(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_38(NCONTR_tot2),           &
                                     NCONTR_38(NCONTR_tot2)
            END DO
         else if(NFILE == 39) then
            READ (NFILE, IOSTAT=IOS) LAB_39(NCONTR_tot2),              &
                                  NCONTR_39(NCONTR_tot2)
            DO WHILE (LAB_39(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_39(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_39(NCONTR_tot+I),                   &
                              INDX_39(NCONTR_tot+I),                   &
                             COEFF_39(NCONTR_tot+I),                   &
                        I=1,NCONTR_39(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_39(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_39(NCONTR_tot2),           &
                                     NCONTR_39(NCONTR_tot2)
            END DO
         else if(NFILE == 40) then
            READ (NFILE, IOSTAT=IOS) LAB_40(NCONTR_tot2),              &
                                  NCONTR_40(NCONTR_tot2)
            DO WHILE (LAB_40(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_40(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_40(NCONTR_tot+I),                   &
                              INDX_40(NCONTR_tot+I),                   &
                             COEFF_40(NCONTR_tot+I),                   &
                        I=1,NCONTR_40(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_40(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_40(NCONTR_tot2),           &
                                     NCONTR_40(NCONTR_tot2)
            END DO
         else if(NFILE == 41) then
            READ (NFILE, IOSTAT=IOS) LAB_41(NCONTR_tot2),              &
                                  NCONTR_41(NCONTR_tot2)
            DO WHILE (LAB_41(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_41(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_41(NCONTR_tot+I),                   &
                              INDX_41(NCONTR_tot+I),                   &
                             COEFF_41(NCONTR_tot+I),                   &
                        I=1,NCONTR_41(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_41(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_41(NCONTR_tot2),           &
                                     NCONTR_41(NCONTR_tot2)
            END DO
         else if(NFILE == 42) then
            READ (NFILE, IOSTAT=IOS) LAB_42(NCONTR_tot2),              &
                                  NCONTR_42(NCONTR_tot2)
            DO WHILE (LAB_42(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_42(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_42(NCONTR_tot+I),                   &
                              INDX_42(NCONTR_tot+I),                   &
                             COEFF_42(NCONTR_tot+I),                   &
                        I=1,NCONTR_42(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_42(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_42(NCONTR_tot2),           &
                                     NCONTR_42(NCONTR_tot2)
            END DO
         else if(NFILE == 43) then
            READ (NFILE, IOSTAT=IOS) LAB_43(NCONTR_tot2),              &
                                  NCONTR_43(NCONTR_tot2)
            DO WHILE (LAB_43(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_43(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_43(NCONTR_tot+I),                   &
                              INDX_43(NCONTR_tot+I),                   &
                             COEFF_43(NCONTR_tot+I),                   &
                        I=1,NCONTR_43(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_43(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_43(NCONTR_tot2),           &
                                     NCONTR_43(NCONTR_tot2)
            END DO
         else if(NFILE == 44) then
            READ (NFILE, IOSTAT=IOS) LAB_44(NCONTR_tot2),              &
                                  NCONTR_44(NCONTR_tot2)
            DO WHILE (LAB_44(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_44(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_44(NCONTR_tot+I),                   &
                              INDX_44(NCONTR_tot+I),                   &
                             COEFF_44(NCONTR_tot+I),                   &
                        I=1,NCONTR_44(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_44(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_44(NCONTR_tot2),           &
                                     NCONTR_44(NCONTR_tot2)
            END DO
         else if(NFILE == 45) then
            READ (NFILE, IOSTAT=IOS) LAB_45(NCONTR_tot2),              &
                                  NCONTR_45(NCONTR_tot2)
            DO WHILE (LAB_45(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_45(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_45(NCONTR_tot+I),                   &
                              INDX_45(NCONTR_tot+I),                   &
                             COEFF_45(NCONTR_tot+I),                   &
                        I=1,NCONTR_45(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_45(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_45(NCONTR_tot2),           &
                                     NCONTR_45(NCONTR_tot2)
            END DO
         else if(NFILE == 46) then
            READ (NFILE, IOSTAT=IOS) LAB_46(NCONTR_tot2),              &
                                  NCONTR_46(NCONTR_tot2)
            DO WHILE (LAB_46(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_46(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_46(NCONTR_tot+I),                   &
                              INDX_46(NCONTR_tot+I),                   &
                             COEFF_46(NCONTR_tot+I),                   &
                        I=1,NCONTR_46(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_46(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_46(NCONTR_tot2),           &
                                     NCONTR_46(NCONTR_tot2)
            END DO
         else if(NFILE == 47) then
            READ (NFILE, IOSTAT=IOS) LAB_47(NCONTR_tot2),              &
                                  NCONTR_47(NCONTR_tot2)
            DO WHILE (LAB_47(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_47(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_47(NCONTR_tot+I),                   &
                              INDX_47(NCONTR_tot+I),                   &
                             COEFF_47(NCONTR_tot+I),                   &
                        I=1,NCONTR_47(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_47(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_47(NCONTR_tot2),           &
                                     NCONTR_47(NCONTR_tot2)
            END DO
         else if(NFILE == 48) then
            READ (NFILE, IOSTAT=IOS) LAB_48(NCONTR_tot2),              &
                                  NCONTR_48(NCONTR_tot2)
            DO WHILE (LAB_48(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_48(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_48(NCONTR_tot+I),                   &
                              INDX_48(NCONTR_tot+I),                   &
                             COEFF_48(NCONTR_tot+I),                   &
                        I=1,NCONTR_48(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_48(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_48(NCONTR_tot2),           &
                                     NCONTR_48(NCONTR_tot2)
            END DO
         else if(NFILE == 49) then
            READ (NFILE, IOSTAT=IOS) LAB_49(NCONTR_tot2),              &
                                  NCONTR_49(NCONTR_tot2)
            DO WHILE (LAB_49(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_49(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_49(NCONTR_tot+I),                   &
                              INDX_49(NCONTR_tot+I),                   &
                             COEFF_49(NCONTR_tot+I),                   &
                        I=1,NCONTR_49(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_49(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_49(NCONTR_tot2),           &
                                     NCONTR_49(NCONTR_tot2)
            END DO
         else if(NFILE == 50) then
            READ (NFILE, IOSTAT=IOS) LAB_50(NCONTR_tot2),              &
                                  NCONTR_50(NCONTR_tot2)
            DO WHILE (LAB_50(NCONTR_tot2) /= 0 .OR.                    &
                   NCONTR_50(NCONTR_tot2) /= 0)
               !*** Read column index, sparse matrix index, and
               !*** coefficient for all contributions from this integral
               READ (NFILE) (ICLMN_50(NCONTR_tot+I),                   &
                              INDX_50(NCONTR_tot+I),                   &
                             COEFF_50(NCONTR_tot+I),                   &
                        I=1,NCONTR_50(NCONTR_tot2))
               NCONTR_tot  = NCONTR_tot + NCONTR_50(NCONTR_tot2)
               NCONTR_tot2 = NCONTR_tot2 + 1
               READ (NFILE, IOSTAT=IOS) LAB_50(NCONTR_tot2),           &
                                     NCONTR_50(NCONTR_tot2)
            END DO
         else
            print*, "Error in setmcp_mem"
            stop
         end if

      END DO
      SYSMemF = SYSTEMFreeMem()  ! GB
      IF(MYID == 0)                                                    &
      WRITE(*,'(A30,F10.2,A3)') 'Free memory on computer',SYSMemF,' Gb'
      CLOSE(NFILE)

      RETURN
      END SUBROUTINE setmcp_mem
