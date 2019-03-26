!*******************************************************************
!                                                                  *
      SUBROUTINE RMEW7BJJ(J1,J2,K1,K2,COEF)
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
      USE ribojj_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rumtjj_I
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER,      INTENT(IN)  :: J1, J2, K1, K2
      REAL(DOUBLE), INTENT(OUT) :: COEF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LL, IQ1, LV1, IL1, IQ2, LV2, IL2, JI1, JI2, IFAZ, J, L
      INTEGER, DIMENSION(6)  :: IPR1
      INTEGER, DIMENSION(8)  :: IPR2
      INTEGER, DIMENSION(21) :: I12PS,I12PV,I03PS,I03PV,I14PS, &
      I14PV,I05PS,I05PV,I16PS,I16PV,I07PS,I07PV
      INTEGER, DIMENSION(36) :: I12AS,I12AV,I03AS,I03AV,I14AS, &
      I14AV,I05AS,I05AV,I16AS,I16AV,I07AS,I07AV
!-----------------------------------------------
      DATA IPR1/0,5,9,12,14,15/
      DATA IPR2/0,7,13,18,22,25,27,28/
      DATA I12PS/-252,1056,144,3*0,507,-88,-325,2*0,200,-520,80,  &
      0,-6125,-560,0,390,-3840,2040/
      DATA I12PV/10,70,7,3*1,10,1,14,2*1,3,21,2*1,462,11,1,539,   &
      2*49/
      DATA I12AS/0,-50,6*0,-1280,990,2640,1950,4*0,-480,4*0,-360, &
      1872,42,390,3*0,48,2*0,234,0,1040,340,0/
      DATA I12AV/8*1,4*49,4*1,49,4*1,539,49,1,11,6*1,7,1,11,7,1/
      DATA I03PS/-1188,-196,0,-234,2*0,189,0,-1911,1470,0,-56,3*0,&
      394805,5250,53760,-78,17408,12920/
      DATA I03PV/70,10,1,7,2*1,110,1,242,121,5*1,22022,121,1573,  &
      2*847,1001/
      DATA I03AS/8*0,110,0,-240,4*0,-1920,0,52,-224,2*0,32490,2*0,&
      -7644,0,12,224,2*0,-52,0,-2040,-364,0,1292/
      DATA I03AV/8*1,7,1,7,4*1,77,1,7,11,2*1,847,2*1,121,1,70,10, &
      2*1,70,1,77,121,1,77/
      DATA I14PS/0,54,-528,546,-378,0,-4335,-96,11271,3822,0,120, &
      12000,-624,-960,-30345,-210,-228480,580476,146880,-627912/
      DATA I14PV/1,2*7,2*11,1,110,11,1694,121,2*1,77,2*11,3146,   &
      121,1573,2*5929,7007/
      DATA I14AS/3*0,-90,4*0,2640,480,-360,20592,-42,390,2*0,     &
      468180,2*0,21840,0,-359424,-6750,0,917280,-10710,2*0,36,2*0,&
      -858,0,-5304,-69768,0/
      DATA I14AV/8*1,2*49,2*539,1,11,2*1,5929,2*1,847,1,5929,539, &
      1,5929,121,2*1,11,2*1,7,1,121,847,1/
      DATA I05PS/3*0,144,14,0,975,0,-50421,336,-7000,-88,3*0,     &
      -6845,-3360,14280,-1836,-103360,-28424/
      DATA I05PV/3*1,7,2*1,22,1,2002,11,143,4*1,26026,143,1859,77,&
      1001,143143/
      DATA I05AS/10*0,390,2*0,-70,3*0,32,14,2*0,-576,2*0,-210,0,  &
      63888,-154,0,-8568,-176,0,1938,1088,0,-28424/
      DATA I05AV/10*1,7,6*1,7,3*1,77,2*1,11,1,1183,13,1,169,7,1,  &
      91,11,1,1183/
      DATA I16PS/3*0,576,390,-144,0,520,-49,-12480,-408,520,      &
      -1960,-1664,3264,552250,-43520,-38760,-2652,15504,38760/
      DATA I16PV/3*1,11,77,7,1,11,2*121,11,3,33,2*11,4719,847,    &
      11011,2*121,143/
      DATA I16AS/6*0,-130,3*0,390,-48,-234,1040,340,0,3120,2*0,   &
      -1170,0,18720,-36,858,-5304,-69768,2*0,1020,4*0,-33592,     &
      31654,0/
      DATA I16AV/10*1,11,1,7,11,7,1,121,2*1,121,1,121,11,7,121,   &
      847,2*1,11,4*1,2*121,1/
      DATA I07PS/4*0,162,272,2*0,11025,4624,-1632,-120,3*0,       &
      1224510,306000,12558240,-6460,77520,-297160/
      DATA I07PV/4*1,2*7,2*1,1573,121,143,4*1,20449,11011,143143, &
      121,1573,1859/
      DATA I07AS/13*0,60,4*0,1600,0,2040,4410,2*0,18360,0,-18816, &
      -11016,0,11628,34,0,7752,9690,0,222870/
      DATA I07AV/18*1,77,1,77,121,2*1,121,1,845,455,1,1183,5,1,   &
      143,121,1,1859/
!
      COEF = ZERO
      LL = 7
      CALL RUMTJJ(J1,LL,IQ1,LV1,IL1)
      CALL RUMTJJ(J2,LL,IQ2,LV2,IL2)
      IF(J1 > J2) THEN
        JI1=J2
        JI2=J1
        IFAZ=IL2-IL1+IQ2-IQ1
        ELSE
          JI1=J1
          JI2=J2
          IFAZ=4
      ENDIF
      IF(J1 > 17) THEN
        JI1=JI1-17
        JI2=JI2-17
        J=IPR2(JI1)+JI2
        L=2
      ELSE
        JI1=JI1-11
        JI2=JI2-11
        L=1
        J=IPR1(JI1)+JI2
      ENDIF
      IF(K1 == 1 .AND. K2 == 2) THEN
        IF(L == 1) THEN
           IF(I12PS(J) >= 0) THEN
              COEF=DSQRT(DBLE(I12PS(J))/DBLE(I12PV(J)))
           ELSE
              COEF=-DSQRT(-DBLE(I12PS(J))/DBLE(I12PV(J)))
           ENDIF
        ELSE
           IF(I12AS(J) >= 0) THEN
              COEF=DSQRT(DBLE(I12AS(J))/DBLE(I12AV(J)))
           ELSE
              COEF=-DSQRT(-DBLE(I12AS(J))/DBLE(I12AV(J)))
           ENDIF
        ENDIF
      ELSEIF(K1 == 0 .AND. K2 == 3) THEN
        IF(L == 1) THEN
           IF(I03PS(J) >= 0) THEN
              COEF=DSQRT(DBLE(I03PS(J))/DBLE(I03PV(J)))
           ELSE
              COEF=-DSQRT(-DBLE(I03PS(J))/DBLE(I03PV(J)))
           ENDIF
        ELSE
           IF(I03AS(J) >= 0) THEN
              COEF=DSQRT(DBLE(I03AS(J))/DBLE(I03AV(J)))
           ELSE
              COEF=-DSQRT(-DBLE(I03AS(J))/DBLE(I03AV(J)))
           ENDIF
        ENDIF
      ELSEIF(K1 == 1 .AND. K2 == 4) THEN
        IF(L == 1) THEN
           IF(I14PS(J) >= 0) THEN
              COEF=DSQRT(DBLE(I14PS(J))/DBLE(I14PV(J)))
           ELSE
              COEF=-DSQRT(-DBLE(I14PS(J))/DBLE(I14PV(J)))
           ENDIF
        ELSE
           IF(I14AS(J) >= 0) THEN
              COEF=DSQRT(DBLE(I14AS(J))/DBLE(I14AV(J)))
           ELSE
              COEF=-DSQRT(-DBLE(I14AS(J))/DBLE(I14AV(J)))
           ENDIF
        ENDIF
      ELSEIF(K1 == 0 .AND. K2 == 5) THEN
        IF(L == 1) THEN
           IF(I05PS(J) >= 0) THEN
              COEF=DSQRT(DBLE(I05PS(J))/DBLE(I05PV(J)))
           ELSE
              COEF=-DSQRT(-DBLE(I05PS(J))/DBLE(I05PV(J)))
           ENDIF
        ELSE
           IF(I05AS(J) >= 0) THEN
              COEF=DSQRT(DBLE(I05AS(J))/DBLE(I05AV(J)))
           ELSE
              COEF=-DSQRT(-DBLE(I05AS(J))/DBLE(I05AV(J)))
           ENDIF
        ENDIF
      ELSEIF(K1 == 1 .AND. K2 == 6) THEN
        IF(L == 1) THEN
           IF(I16PS(J) >= 0) THEN
              COEF=DSQRT(DBLE(I16PS(J))/DBLE(I16PV(J)))
           ELSE
              COEF=-DSQRT(-DBLE(I16PS(J))/DBLE(I16PV(J)))
           ENDIF
        ELSE
           IF(I16AS(J) >= 0) THEN
              COEF=DSQRT(DBLE(I16AS(J))/DBLE(I16AV(J)))
           ELSE
              COEF=-DSQRT(-DBLE(I16AS(J))/DBLE(I16AV(J)))
           ENDIF
        ENDIF
      ELSEIF(K1 == 0 .AND. K2 == 7) THEN
        IF(L == 1) THEN
           IF(I07PS(J) >= 0) THEN
              COEF=DSQRT(DBLE(I07PS(J))/DBLE(I07PV(J)))
           ELSE
              COEF=-DSQRT(-DBLE(I07PS(J))/DBLE(I07PV(J)))
           ENDIF
        ELSE
           IF(I07AS(J) >= 0) THEN
              COEF=DSQRT(DBLE(I07AS(J))/DBLE(I07AV(J)))
           ELSE
              COEF=-DSQRT(-DBLE(I07AS(J))/DBLE(I07AV(J)))
           ENDIF
        ENDIF
      ELSE
        WRITE(0,'(A,4I5)') ' J1 J2 = ',J1,J2
        WRITE(0,'(A)') ' ERROR IN SUB. RMEW7JJ '
        STOP
      ENDIF
      IF(MOD(IFAZ,4) /= 0)COEF=-COEF
      RETURN
      END
