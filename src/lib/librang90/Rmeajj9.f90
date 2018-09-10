!*******************************************************************
!                                                                  *
      SUBROUTINE RMEAJJ9(IT,LQ,J,ITS,LQS,J1S,COEF)
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY: DOUBLE
      USE CONS_C,          ONLY: ZERO
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)  :: IT, LQ, J, ITS, LQS, J1S
      REAL(DOUBLE), INTENT(OUT) :: COEF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: JT, JTS
      INTEGER, DIMENSION(20) :: IDS3,IDS4,IDS5,IDS6,IDS7,IDV7,IDS8,&
      IDV8,IDS9,IDS10,IDS11,IDV11,IDS12,IDV12,IDS13,IDV13,IDS14,   &
      IDS15,IDV15,IDS16,IDV16,IDS17,IDV17,IDS18,IDV18
!-----------------------------------------------
      DATA IDS3/2*0,660,0,1664,0,1650,130,0,816,0,1680,8*0/
      DATA IDS4/2*0,-55770,-20592,6292,-13260,0,12740,-19110,     &
      79968,-47880,-25410,-19278,7*0/
      DATA IDS5/0,-9009,-4056,-260,585,16575,0,25200,4200,-7140,  &
      -475,504,-14280,13566,-4250,5*0/
      DATA IDS6/0,-4992,2028,0,-1920,0,-12870,7350,0,8160,        &
      0,1512,0,-3648,0,-9000,4*0/
      DATA IDS7/-2184,-63,-59904,-302460,5265,848691,0,-145152,   &
      217728,1049580,287337,5184,261120,691866,-750,0,-76608,3*0/
      DATA IDV7/253,23,19481,19481,161,253253,1,36179,            &
      36179,2*36179,299,3*36179,1,3289,3*1/
      DATA IDS8/1224,2652,188598,-31824,204,-12996,0,25500,-38250,&
      -768,81396,3213,-60543,-3066144,727776,207,-41553,3*0/
      DATA IDV8/1265,115,13915,2783,23,13915,1,3*2783,13915,115,  &
      2783,236555,47311,17,21505,3*1/
      DATA IDS9/3380,-2340,5460,-1400,3150,-3570,0,9000,1500,2550,&
      -5320,-10710,1050,1140,6300,8550,-330,5750,2*0/
      DATA IDS10/0,2160,6240,0,-32,0,-21450,-9610,0,2688,         &
      0,7140,0,20160,0,6156,2*0,-10164,0/
      DATA IDS11/0,132,-52728,196,-50,-24990,0,-21160,31740,26250,&
      12920,-357,9583,-344988,5700,171,-15,-4830,-84,0/
      DATA IDV11/1,5,6655,1331,11,1331,1,4*1331,55,1331,6655,     &
      1331,11,2*121,11,1/
      DATA IDS12/2*0,-209950,77520,12920,-25688,0,-4522,2261,     &
      -48640,-285144,931,273885,-112908,-2138580,137781,6654375,  &
      -59616,284089,0/
      DATA IDV12/2*1,2*9317,231,9317,1,3993,1331,22627,9317,44,   &
      90508,429913,158389,3740,156332,39083,17765,1/
      DATA IDS13/2*0,1530,13056,29376,720,0,-1890,-315,101124,    &
      -4560,-13965,-35131,13500,-685900,-1197,-28875,-5060,759,0/
      DATA IDV13/2*1,77,3*1001,1,2*143,2431,1001,572,9724,2431,   &
      17017,2*884,221,17,1/
      DATA IDS14/4*0,22848,0,-121550,4590,0,-32832,0,45220,0,     &
      -31680,0,82764,2*0,144716,0/
      DATA IDS15/5*0,2128,0,-1938,2907,-4860,-17136,-1309,-8505,  &
      15876,420,1287,-6075,-132,-253,-650/
      DATA IDV15/5*1,143,1,3*143,2717,52,572,247,13,20,988,19,95, &
      19/
      DATA IDS16/7*0,570,95,4104,504,-1463,-39501,-60516,-1596,   &
      3933,621,-840,-805,390/
      DATA IDV16/7*1,2*13,221,65,260,884,1105,221,68,3740,2*17,   &
      11/
      DATA IDS17/9*0,-5796,17664,-5313,1771,-16632,-88,693,       &
      -94269,192500,30030,24570/
      DATA IDV17/9*1,221,1235,65,221,1615,2*17,1615,7429,323,437/
      DATA IDS18/13*0,-15000,280,-1170,-48750,-15600,3510,-33930/
      DATA IDV18/13*1,323,2*17,3553,437,19,253/
!
      COEF=ZERO
      IF(IT < ITS) THEN
        JT=IT-25
        JTS=ITS-45
      ELSE
        JT=ITS-25
        JTS=IT-45
      ENDIF
      IF(JTS == 1) THEN
        IF(JT == 7) COEF=-DSQRT(DBLE(60))
      ELSEIF(JTS == 2) THEN
        IF(JT == 8) COEF=-DSQRT(DBLE(12))
        IF(JT == 9) COEF=-DSQRT(DBLE(8))
      ELSEIF(JTS == 3)THEN
        COEF=-DSQRT(DBLE(IDS3(JT))/DBLE(33))
      ELSEIF(JTS == 4)THEN
        IF(IDS4(JT) >= 0) THEN
           COEF=DSQRT(DBLE(IDS4(JT))/DBLE(3003))
        ELSE
           COEF=-DSQRT(-DBLE(IDS4(JT))/DBLE(3003))
        ENDIF
      ELSEIF(JTS == 5)THEN
        IF(IDS5(JT) >= 0) THEN
           COEF=DSQRT(DBLE(IDS5(JT))/DBLE(715))
        ELSE
           COEF=-DSQRT(-DBLE(IDS5(JT))/DBLE(715))
        ENDIF
      ELSEIF(JTS == 6)THEN
        IF(IDS6(JT) >= 0) THEN
           COEF=DSQRT(DBLE(IDS6(JT))/DBLE(143))
        ELSE
           COEF=-DSQRT(-DBLE(IDS6(JT))/DBLE(143))
        ENDIF
      ELSEIF(JTS == 7)THEN
        IF(IDS7(JT) >= 0) THEN
           COEF=DSQRT(DBLE(IDS7(JT))/DBLE(IDV7(JT)))
        ELSE
           COEF=-DSQRT(-DBLE(IDS7(JT))/DBLE(IDV7(JT)))
        ENDIF
      ELSEIF(JTS == 8)THEN
        IF(IDS8(JT) >= 0) THEN
           COEF=DSQRT(DBLE(IDS8(JT))/DBLE(IDV8(JT)))
        ELSE
           COEF=-DSQRT(-DBLE(IDS8(JT))/DBLE(IDV8(JT)))
        ENDIF
      ELSEIF(JTS == 9)THEN
        IF(IDS9(JT) >= 0) THEN
           COEF=DSQRT(DBLE(IDS9(JT))/DBLE(325))
        ELSE
           COEF=-DSQRT(-DBLE(IDS9(JT))/DBLE(325))
        ENDIF
      ELSEIF(JTS == 10)THEN
        IF(IDS10(JT) >= 0) THEN
           COEF=DSQRT(DBLE(IDS10(JT))/DBLE(165))
        ELSE
           COEF=-DSQRT(-DBLE(IDS10(JT))/DBLE(165))
        ENDIF
      ELSEIF(JTS == 11)THEN
        IF(IDS11(JT) >= 0) THEN
           COEF=DSQRT(DBLE(IDS11(JT))/DBLE(IDV11(JT)))
        ELSE
           COEF=-DSQRT(-DBLE(IDS11(JT))/DBLE(IDV11(JT)))
        ENDIF
      ELSEIF(JTS == 12)THEN
        IF(IDS12(JT) >= 0) THEN
           COEF=DSQRT(DBLE(IDS12(JT))/DBLE(IDV12(JT)))
        ELSE
           COEF=-DSQRT(-DBLE(IDS12(JT))/DBLE(IDV12(JT)))
        ENDIF
      ELSEIF(JTS == 13)THEN
        IF(IDS13(JT) >= 0) THEN
           COEF=DSQRT(DBLE(IDS13(JT))/DBLE(IDV13(JT)))
        ELSE
           COEF=-DSQRT(-DBLE(IDS13(JT))/DBLE(IDV13(JT)))
        ENDIF
      ELSEIF(JTS == 14)THEN
        IF(IDS14(JT) >= 0) THEN
           COEF=DSQRT(DBLE(IDS14(JT))/DBLE(715))
        ELSE
           COEF=-DSQRT(-DBLE(IDS14(JT))/DBLE(715))
        ENDIF
      ELSEIF(JTS == 15)THEN
        IF(IDS15(JT) >= 0) THEN
           COEF=DSQRT(DBLE(IDS15(JT))/DBLE(IDV15(JT)))
        ELSE
           COEF=-DSQRT(-DBLE(IDS15(JT))/DBLE(IDV15(JT)))
        ENDIF
      ELSEIF(JTS == 16)THEN
        IF(IDS16(JT) >= 0) THEN
           COEF=DSQRT(DBLE(IDS16(JT))/DBLE(IDV16(JT)))
        ELSE
           COEF=-DSQRT(-DBLE(IDS16(JT))/DBLE(IDV16(JT)))
        ENDIF
      ELSEIF(JTS == 17)THEN
        IF(IDS17(JT) >= 0) THEN
           COEF=DSQRT(DBLE(IDS17(JT))/DBLE(IDV17(JT)))
        ELSE
           COEF=-DSQRT(-DBLE(IDS17(JT))/DBLE(IDV17(JT)))
        ENDIF
      ELSEIF(JTS == 18)THEN
        IF(IDS18(JT) >= 0) THEN
           COEF=DSQRT(DBLE(IDS18(JT))/DBLE(IDV18(JT)))
        ELSE
           COEF=-DSQRT(-DBLE(IDS18(JT))/DBLE(IDV18(JT)))
        ENDIF
      ELSE
        WRITE(0,'(A,4I5)') ' IT ITS JT JTS= ',IT,ITS,JT,JTS
        WRITE(0,'(A)') ' ERROR IN SUB. RMEAJJ9 '
        STOP
      ENDIF
      RETURN
      END SUBROUTINE RMEAJJ9
