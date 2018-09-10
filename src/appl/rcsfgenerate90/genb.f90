!     last edited July 31, 1996
      SUBROUTINE GEN(ANSATS, POSN, POSL, SKAL, CF, FIRST, MINJ, MAXJ, PAR) 
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE kopp1_I 
      USE kopp2_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: SKAL 
      INTEGER , INTENT(INOUT) :: CF 
      INTEGER , INTENT(IN) :: MINJ 
      INTEGER , INTENT(IN) :: MAXJ 
      INTEGER  :: PAR 
      LOGICAL , INTENT(IN) :: FIRST 
      INTEGER , INTENT(IN) :: ANSATS(15,0:10,0:1) 
      INTEGER , INTENT(IN) :: POSN(110) 
      INTEGER , INTENT(IN) :: POSL(110) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: FIL_1 = 7 
      INTEGER, PARAMETER :: FIL_2 = 8 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(0:10,0:11,0:1) :: KOPPL 
      INTEGER , DIMENSION(0:10,0:1,0:5,20) :: JKVANT 
      INTEGER , DIMENSION(0:10,0:1) :: ANTMAX 
      INTEGER :: POS, I, N, L, K 
      INTEGER , DIMENSION(20) :: J, JK, ORBIT, ANTEL 
      INTEGER :: I1, I2, I3, I4, I5, I6, I7, I8, I9, I10, I11, I12, I13, I14, &
         I15, I16, I17, I18, I19, I20 
      INTEGER , DIMENSION(20) :: PLUS, S 
      INTEGER :: RESJ, JK1, JK2, JK3, JK4, JK5, JK6, JK7, JK8, JK9, JK10, JK11&
         , JK12, JK13, JK14, JK15, JK16, JK17, JK18, FIL 
      INTEGER , DIMENSION(20) :: ANTKO 
      INTEGER , DIMENSION(0:10,0:1,0:5,20) :: SENIOR 
      INTEGER :: N1, N10 
      CHARACTER :: RAD1*200, RAD2*200, RAD3*200 
      CHARACTER, DIMENSION(0:10,0:1) :: L1*2 
!-----------------------------------------------
      DATA (L1(I,0),I=0,10)/ 's ', 'p-', 'd-', 'f-', 'g-', 'h-', 'i-', 'k-', &
         'l-', 'm-', 'n-'/  
      DATA (L1(I,1),I=0,10)/ 's ', 'p ', 'd ', 'f ', 'g ', 'h ', 'i ', 'k ', &
         'l ', 'm ', 'n '/  
!     The value of antmax(l-number,x) is the maximum number of electrons
!     in the orbital, x represents +/- coupling of s- and l- number
      DATA (ANTMAX(I,0),I=0,10)/ 2, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20/  
      DATA (ANTMAX(I,1),I=0,10)/ 0, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22/  
!     The value of koppl(l-number,number of electrons,x) is the number of
!     possible couplings for a certain orbital. If the orbital is
!     populated with more than half of the maximal number of electrons
!     the index "number of electrons" should be substituted with
!     "antmax(l-number) - number of electrons".
      DATA (KOPPL(0,I,0),I=0,1)/ 1, 1/  
!     l=0
      DATA (KOPPL(1,I,0),I=0,1)/ 1, 1/  
      DATA (KOPPL(1,I,1),I=0,2)/ 1, 1, 2/  
!     l=1
      DATA (KOPPL(2,I,0),I=0,2)/ 1, 1, 2/  
      DATA (KOPPL(2,I,1),I=0,3)/ 1, 1, 3, 3/  
!     l=2
      DATA (KOPPL(3,I,0),I=0,3)/ 1, 1, 3, 3/  
      DATA (KOPPL(3,I,1),I=0,4)/ 1, 1, 4, 6, 8/  
!     l=3
      DATA (KOPPL(4,I,0),I=0,4)/ 1, 1, 4, 6, 8/  
      DATA (KOPPL(4,I,1),I=0,5)/ 1, 1, 5, 10, 16, 20/  
!     l=4
      DATA (KOPPL(5,I,0),I=0,5)/ 1, 1, 5, 10, 16, 20/  
      DATA (KOPPL(5,I,1),I=0,2)/ 1, 1, 6/  
!     l=5
      DATA (KOPPL(6,I,0),I=0,2)/ 1, 1, 6/  
      DATA (KOPPL(6,I,1),I=0,2)/ 1, 1, 7/  
!     l=6
      DATA (KOPPL(7,I,0),I=0,2)/ 1, 1, 7/  
      DATA (KOPPL(7,I,1),I=0,2)/ 1, 1, 8/  
!     l=7
      DATA (KOPPL(8,I,0),I=0,2)/ 1, 1, 8/  
      DATA (KOPPL(8,I,1),I=0,2)/ 1, 1, 9/  
!     l=8
      DATA (KOPPL(9,I,0),I=0,2)/ 1, 1, 9/  
      DATA (KOPPL(9,I,1),I=0,2)/ 1, 1, 10/  
!     l=9
      DATA (KOPPL(10,I,0),I=0,2)/ 1, 1, 10/  
      DATA (KOPPL(10,I,1),I=0,2)/ 1, 1, 11/  
!     l=10
 
!  JKVANT(l-number, +/-, number of electrons, coupling number) is 2*J-number
 
      DATA JKVANT(0,0,0,1)/ 0/  
!     data  SENIOR(0,0,0,1)           / 0/
      DATA SENIOR(0,0,0,1)/ -1/  
!     l=0 #=0
      DATA JKVANT(0,0,1,1)/ 1/  
!     data  SENIOR(0,0,1,1)           / 1/
      DATA SENIOR(0,0,1,1)/ -1/  
!     l=0 #=1
      DATA JKVANT(1,0,0,1)/ 0/  
!     data  SENIOR(1,0,0,1)           / 0/
      DATA SENIOR(1,0,0,1)/ -1/  
!     l=1 #=0 -
      DATA JKVANT(1,0,1,1)/ 1/  
!     data  SENIOR(1,0,1,1)           / 1/
      DATA SENIOR(1,0,1,1)/ -1/  
!     l=1 #=1 -
      DATA JKVANT(1,1,0,1)/ 0/  
!     data  SENIOR(1,1,0,1)           / 0/
      DATA SENIOR(1,1,0,1)/ -1/  
!     l=1 #=0 +
      DATA JKVANT(1,1,1,1)/ 3/  
!     data  SENIOR(1,1,1,1)           / 1/
      DATA SENIOR(1,1,1,1)/ -1/  
!     l=1 #=1 +
      DATA (JKVANT(1,1,2,I),I=1,2)/ 0, 4/  
!     data (SENIOR(1,1,2,i),i=1,2)    / 0, 2/
      DATA (SENIOR(1,1,2,I),I=1,2)/ -1, -1/  
!     l=1 #=2 +
      DATA JKVANT(2,0,0,1)/ 0/  
!     data  SENIOR(2,0,0,1)           / 0/
      DATA SENIOR(2,0,0,1)/ -1/  
!     l=2 #=0 -
      DATA JKVANT(2,0,1,1)/ 3/  
!     data  SENIOR(2,0,1,1)           / 1/
      DATA SENIOR(2,0,1,1)/ -1/  
!     l=2 #=1 -
      DATA (JKVANT(2,0,2,I),I=1,2)/ 0, 4/  
!     data (SENIOR(2,0,2,i),i=1,2)    / 0, 2/
      DATA (SENIOR(2,0,2,I),I=1,2)/ -1, -1/  
!     l=2 #=2 -
      DATA JKVANT(2,1,0,1)/ 0/  
!     data  SENIOR(2,1,0,1)           / 0/
      DATA SENIOR(2,1,0,1)/ -1/  
!     l=2 #=0 +
      DATA JKVANT(2,1,1,1)/ 5/  
!     data  SENIOR(2,1,1,1)           / 1/
      DATA SENIOR(2,1,1,1)/ -1/  
!     l=2 #=1 +
      DATA (JKVANT(2,1,2,I),I=1,3)/ 0, 4, 8/  
!     data (SENIOR(2,1,2,i),i=1,3)    / 0, 2, 2/
      DATA (SENIOR(2,1,2,I),I=1,3)/ -1, -1, -1/  
!     l=2 #=2 +
      DATA (JKVANT(2,1,3,I),I=1,3)/ 5, 3, 9/  
!     data (SENIOR(2,1,3,i),i=1,3)    / 1, 3, 3/
      DATA (SENIOR(2,1,3,I),I=1,3)/ -1, -1, -1/  
!     l=2 #=3 +
      DATA JKVANT(3,0,0,1)/ 0/  
!     data  SENIOR(3,0,0,1)           / 0/
      DATA SENIOR(3,0,0,1)/ -1/  
!     l=3 #=0 -
      DATA JKVANT(3,0,1,1)/ 5/  
!     data  SENIOR(3,0,1,1)           / 1/
      DATA SENIOR(3,0,1,1)/ -1/  
!     l=3 #=1 -
      DATA (JKVANT(3,0,2,I),I=1,3)/ 0, 4, 8/  
!     data (SENIOR(3,0,2,i),i=1,3)    / 0, 2, 2/
      DATA (SENIOR(3,0,2,I),I=1,3)/ -1, -1, -1/  
!     l=3 #=2 -
      DATA (JKVANT(3,0,3,I),I=1,3)/ 5, 3, 9/  
!     data (SENIOR(3,0,3,i),i=1,3)    / 1, 3, 3/
      DATA (SENIOR(3,0,3,I),I=1,3)/ -1, -1, -1/  
!     l=3 #=3 -
      DATA JKVANT(3,1,0,1)/ 0/  
!     data  SENIOR(3,1,0,1)           / 0/
      DATA SENIOR(3,1,0,1)/ -1/  
!     l=3 #=0 +
      DATA JKVANT(3,1,1,1)/ 7/  
!     data  SENIOR(3,1,1,1)           / 1/
      DATA SENIOR(3,1,1,1)/ -1/  
!     l=3 #=1 +
      DATA (JKVANT(3,1,2,I),I=1,4)/ 0, 4, 8, 12/  
!     data (SENIOR(3,1,2,i),i=1,4)    / 0, 2, 2, 2/
      DATA (SENIOR(3,1,2,I),I=1,4)/ -1, -1, -1, -1/  
!     l=3 #=2 +
      DATA (JKVANT(3,1,3,I),I=1,6)/ 7, 3, 5, 9, 11, 15/  
!     data (SENIOR(3,1,3,i),i=1,6)    / 1, 3, 3, 3, 3, 3/
      DATA (SENIOR(3,1,3,I),I=1,6)/ -1, -1, -1, -1, -1, -1/  
!     l=3 #=3 +
      DATA (JKVANT(3,1,4,I),I=1,8)/ 0, 4, 8, 12, 4, 8, 10, 16/  
!     data (SENIOR(3,1,4,i),i=1,8)    / 0, 2, 2, 2, 4, 4, 4, 4/
      DATA (SENIOR(3,1,4,I),I=1,8)/ -1, 2, 2, -1, 4, 4, -1, -1/  
!     l=3 #=4 +
      DATA JKVANT(4,0,0,1)/ 0/  
!     data  SENIOR(4,0,0,1)           / 0/
      DATA SENIOR(4,0,0,1)/ -1/  
!     l=4 #=0 -
      DATA JKVANT(4,0,1,1)/ 7/  
!     data  SENIOR(4,0,1,1)           / 1/
      DATA SENIOR(4,0,1,1)/ -1/  
!     l=4 #=1 -
      DATA (JKVANT(4,0,2,I),I=1,4)/ 0, 4, 8, 12/  
!     data (SENIOR(4,0,2,i),i=1,4)    / 0, 2, 2, 2/
      DATA (SENIOR(4,0,2,I),I=1,4)/ -1, -1, -1, -1/  
!     l=4 #=2 -
      DATA (JKVANT(4,0,3,I),I=1,6)/ 7, 3, 5, 9, 11, 15/  
!     data (SENIOR(4,0,3,i),i=1,6)    / 1, 3, 3, 3, 3, 3/
      DATA (SENIOR(4,0,3,I),I=1,6)/ -1, -1, -1, -1, -1, -1/  
!     l=4 #=3 -
      DATA (JKVANT(4,0,4,I),I=1,8)/ 0, 4, 8, 12, 4, 8, 10, 16/  
!     data (SENIOR(4,0,4,i),i=1,8)    / 0, 2, 2, 2, 4, 4, 4, 4/
      DATA (SENIOR(4,0,4,I),I=1,8)/ -1, 2, 2, -1, 4, 4, -1, -1/  
!     l=4 #=4 -
      DATA JKVANT(4,1,0,1)/ 0/  
!     data  SENIOR(4,1,0,1)           / 0/
      DATA SENIOR(4,1,0,1)/ -1/  
!     l=4 #=0 +
      DATA JKVANT(4,1,1,1)/ 9/  
!     data  SENIOR(4,1,1,1)           / 1/
      DATA SENIOR(4,1,1,1)/ -1/  
!     l=4 #=1 +
      DATA (JKVANT(4,1,2,I),I=1,5)/ 0, 4, 8, 12, 16/  
!     data (SENIOR(4,1,2,i),i=1,5)    / 0, 2, 2, 2, 2/
      DATA (SENIOR(4,1,2,I),I=1,5)/ -1, -1, -1, -1, -1/  
!     l=4 #=2 +
      DATA (JKVANT(4,1,3,I),I=1,10)/ 9, 3, 5, 7, 9, 11, 13, 15, 17, 21/  
!     data (SENIOR(4,1,3,i),i=1,10)   / 1, 3, 3, 3, 3, 3, 3, 3, 3, 3/
      DATA (SENIOR(4,1,3,I),I=1,10)/ 1, -1, -1, -1, 3, -1, -1, -1, -1, -1/  
!     l=4 #=3 +
      DATA (JKVANT(4,1,4,I),I=1,16)/ 0, 4, 8, 12, 16, 0, 4, 6, 8, 10, 12, 14, &
         16, 18, 20, 24/  
!     data (SENIOR(4,1,4,i),i=1,16)   / 0, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4,
!    :                                  4, 4, 4, 4, 4/
      DATA (SENIOR(4,1,4,I),I=1,16)/ 0, 2, 2, 2, 2, 4, 4, -1, 4, -1, 4, -1, 4, &
         -1, -1, -1/  
!     l=4 #=4 +
      DATA (JKVANT(4,1,5,I),I=1,20)/ 9, 3, 5, 7, 9, 11, 13, 15, 17, 21, 1, 5, 7&
         , 9, 11, 13, 15, 17, 19, 25/  
!     data (SENIOR(4,1,5,i),i=1,20)   / 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 5,
!    :                                  5, 5, 5, 5, 5, 5, 5, 5, 5/
      DATA (SENIOR(4,1,5,I),I=1,20)/ 1, -1, 3, 3, 3, 3, 3, 3, 3, -1, -1, 5, 5, &
         5, 5, 5, 5, 5, -1, -1/  
!     l=4 #=5 +
      DATA JKVANT(5,0,0,1)/ 0/  
!     data  SENIOR(5,0,0,1)           / 0/
      DATA SENIOR(5,0,0,1)/ -1/  
!     l=5 #=0 -
      DATA JKVANT(5,0,1,1)/ 9/  
!     data  SENIOR(5,0,1,1)           / 1/
      DATA SENIOR(5,0,1,1)/ -1/  
!     l=5 #=1 -
      DATA (JKVANT(5,0,2,I),I=1,5)/ 0, 4, 8, 12, 16/  
!     data (SENIOR(5,0,2,i),i=1,5)    / 0, 2, 2, 2, 2/
      DATA (SENIOR(5,0,2,I),I=1,5)/ -1, -1, -1, -1, -1/  
!     l=5 #=2 -
      DATA (JKVANT(5,0,3,I),I=1,10)/ 9, 3, 5, 7, 9, 11, 13, 15, 17, 21/  
!     data (SENIOR(5,0,3,i),i=1,10)   / 1, 3, 3, 3, 3, 3, 3, 3, 3, 3/
      DATA (SENIOR(5,0,3,I),I=1,10)/ 1, -1, -1, -1, 3, -1, -1, -1, -1, -1/  
!     l=5 #=3 -
      DATA (JKVANT(5,0,4,I),I=1,16)/ 0, 4, 8, 12, 16, 0, 4, 6, 8, 10, 12, 14, &
         16, 18, 20, 24/  
!     data (SENIOR(5,0,4,i),i=1,16)   / 0, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4,
!    :                                  4, 4, 4, 4, 4/
      DATA (SENIOR(5,0,4,I),I=1,16)/ 0, 2, 2, 2, 2, 4, 4, -1, 4, -1, 4, -1, 4, &
         -1, -1, -1/  
!     l=5 #=4 -
      DATA (JKVANT(5,0,5,I),I=1,20)/ 9, 3, 5, 7, 9, 11, 13, 15, 17, 21, 1, 5, 7&
         , 9, 11, 13, 15, 17, 19, 25/  
!     data (SENIOR(5,0,5,i),i=1,20)   / 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 5,
!    :                                  5, 5, 5, 5, 5, 5, 5, 5, 5/
      DATA (SENIOR(5,0,5,I),I=1,20)/ 1, -1, 3, 3, 3, 3, 3, 3, 3, -1, -1, 5, 5, &
         5, 5, 5, 5, 5, -1, -1/  
!     l=5 #=5 -
      DATA JKVANT(5,1,0,1)/ 0/  
!     data  SENIOR(5,1,0,1)           / 0/
      DATA SENIOR(5,1,0,1)/ -1/  
!     l=5 #=0 +
      DATA JKVANT(5,1,1,1)/ 11/  
!     data  SENIOR(5,1,1,1)           / 1/
      DATA SENIOR(5,1,1,1)/ -1/  
!     l=5 #=1 +
      DATA (JKVANT(5,1,2,I),I=1,6)/ 0, 4, 8, 12, 16, 20/  
!     data (SENIOR(5,1,2,i),i=1,6)    / 0, 2, 2, 2, 2, 2/
      DATA (SENIOR(5,1,2,I),I=1,6)/ -1, -1, -1, -1, -1, -1/  
!     l=5 #=2 +
      DATA JKVANT(6,0,0,1)/ 0/  
!     data  SENIOR(6,0,0,1)           / 0/
      DATA SENIOR(6,0,0,1)/ -1/  
!     l=6 #=0 -
      DATA JKVANT(6,0,1,1)/ 11/  
!     data  SENIOR(6,0,1,1)           / 1/
      DATA SENIOR(6,0,1,1)/ -1/  
!     l=6 #=1 -
      DATA (JKVANT(6,0,2,I),I=1,6)/ 0, 4, 8, 12, 16, 20/  
!     data (SENIOR(6,0,2,i),i=1,6)    / 0, 2, 2, 2, 2, 2/
      DATA (SENIOR(6,0,2,I),I=1,6)/ -1, -1, -1, -1, -1, -1/  
!     l=6 #=2 -
      DATA JKVANT(6,1,0,1)/ 0/  
!     data  SENIOR(6,1,0,1)           / 0/
      DATA SENIOR(6,1,0,1)/ -1/  
!     l=6 #=0 +
      DATA JKVANT(6,1,1,1)/ 13/  
!     data  SENIOR(6,1,1,1)           / 1/
      DATA SENIOR(6,1,1,1)/ -1/  
!     l=6 #=1 +
      DATA (JKVANT(6,1,2,I),I=1,7)/ 0, 4, 8, 12, 16, 20, 24/  
!     data (SENIOR(6,1,2,i),i=1,7)    / 0, 2, 2, 2, 2, 2, 2/
      DATA (SENIOR(6,1,2,I),I=1,7)/ -1, -1, -1, -1, -1, -1, -1/  
!     l=6 #=2 +
      DATA JKVANT(7,0,0,1)/ 0/  
!     data  SENIOR(7,0,0,1)           / 0/
      DATA SENIOR(7,0,0,1)/ -1/  
!     l=7 #=0 -
      DATA JKVANT(7,0,1,1)/ 13/  
!     data  SENIOR(7,0,1,1)           / 1/
      DATA SENIOR(7,0,1,1)/ -1/  
!     l=7 #=1 -
      DATA (JKVANT(7,0,2,I),I=1,7)/ 0, 4, 8, 12, 16, 20, 24/  
!     data (SENIOR(7,0,2,i),i=1,7)    / 0, 2, 2, 2, 2, 2, 2/
      DATA (SENIOR(7,0,2,I),I=1,7)/ -1, -1, -1, -1, -1, -1, -1/  
!     l=7 #=2 -
      DATA JKVANT(7,1,0,1)/ 0/  
!     data  SENIOR(7,1,0,1)           / 0/
      DATA SENIOR(7,1,0,1)/ -1/  
!     l=7 #=0 +
      DATA JKVANT(7,1,1,1)/ 15/  
!     data  SENIOR(7,1,1,1)           / 1/
      DATA SENIOR(7,1,1,1)/ -1/  
!     l=7 #=1 +
      DATA (JKVANT(7,1,2,I),I=1,8)/ 0, 4, 8, 12, 16, 20, 24, 28/  
!     data (SENIOR(7,1,2,i),i=1,8)    / 0, 2, 2, 2, 2, 2, 2, 2/
      DATA (SENIOR(7,1,2,I),I=1,8)/ -1, -1, -1, -1, -1, -1, -1, -1/  
!     l=7 #=2 +
      DATA JKVANT(8,0,0,1)/ 0/  
!     data  SENIOR(8,0,0,1)           / 0/
      DATA SENIOR(8,0,0,1)/ -1/  
!     l=8 #=0 -
      DATA JKVANT(8,0,1,1)/ 15/  
!     data  SENIOR(8,0,1,1)           / 1/
      DATA SENIOR(8,0,1,1)/ -1/  
!     l=8 #=1 -
      DATA (JKVANT(8,0,2,I),I=1,8)/ 0, 4, 8, 12, 16, 20, 24, 28/  
!     data (SENIOR(8,0,2,i),i=1,8)    / 0, 2, 2, 2, 2, 2, 2, 2/
      DATA (SENIOR(8,0,2,I),I=1,8)/ -1, -1, -1, -1, -1, -1, -1, -1/  
!     l=8 #=2 -
      DATA JKVANT(8,1,0,1)/ 0/  
!     data  SENIOR(8,1,0,1)           / 0/
      DATA SENIOR(8,1,0,1)/ -1/  
!     l=8 #=0 +
      DATA JKVANT(8,1,1,1)/ 17/  
!     data  SENIOR(8,1,1,1)           / 1/
      DATA SENIOR(8,1,1,1)/ -1/  
!     l=8 #=1 +
      DATA (JKVANT(8,1,2,I),I=1,9)/ 0, 4, 8, 12, 16, 20, 24, 28, 32/  
!     data (SENIOR(8,1,2,i),i=1,9)    / 0, 2, 2, 2, 2, 2, 2, 2, 2/
      DATA (SENIOR(8,1,2,I),I=1,9)/ -1, -1, -1, -1, -1, -1, -1, -1, -1/  
!     l=8 #=2 +
      DATA JKVANT(9,0,0,1)/ 0/  
!     data  SENIOR(9,0,0,1)           / 0/
      DATA SENIOR(9,0,0,1)/ -1/  
!     l=9 #=0 -
      DATA JKVANT(9,0,1,1)/ 17/  
!     data  SENIOR(9,0,1,1)           / 1/
      DATA SENIOR(9,0,1,1)/ -1/  
!     l=9 #=1 -
      DATA (JKVANT(9,0,2,I),I=1,9)/ 0, 4, 8, 12, 16, 20, 24, 28, 32/  
!     data (SENIOR(9,0,2,i),i=1,9)    / 0, 2, 2, 2, 2, 2, 2, 2, 2/
      DATA (SENIOR(9,0,2,I),I=1,9)/ -1, -1, -1, -1, -1, -1, -1, -1, -1/  
!     l=9 #=2 -
      DATA JKVANT(9,1,0,1)/ 0/  
!     data  SENIOR(9,1,0,1)           / 0/
      DATA SENIOR(9,1,0,1)/ -1/  
!     l=9 #=0 +
      DATA JKVANT(9,1,1,1)/ 19/  
!     data  SENIOR(9,1,1,1)           / 1/
      DATA SENIOR(9,1,1,1)/ -1/  
!     l=9 #=1 +
      DATA (JKVANT(9,1,2,I),I=1,10)/ 0, 4, 8, 12, 16, 20, 24, 28, 32, 36/  
!     data (SENIOR(9,1,2,i),i=1,10)   / 0, 2, 2, 2, 2, 2, 2, 2, 2, 2/
      DATA (SENIOR(9,1,2,I),I=1,10)/ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1/  
!     l=9 #=2 +
      DATA JKVANT(10,0,0,1)/ 0/  
!     data  SENIOR(10,0,0,1)          / 0/
      DATA SENIOR(10,0,0,1)/ -1/  
!     l=10 #=0 -
      DATA JKVANT(10,0,1,1)/ 19/  
!     data  SENIOR(10,0,1,1)          / 1/
      DATA SENIOR(10,0,1,1)/ -1/  
!     l=10 #=1 -
      DATA (JKVANT(10,0,2,I),I=1,10)/ 0, 4, 8, 12, 16, 20, 24, 28, 32, 36/  
!     data (SENIOR(10,0,2,i),i=1,10)   / 0, 2, 2, 2, 2, 2, 2, 2, 2, 2/
      DATA (SENIOR(10,0,2,I),I=1,10)/ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1/  
!     l=10 #=2 -
      DATA JKVANT(10,1,0,1)/ 0/  
!     data  SENIOR(10,1,0,1)          / 0/
      DATA SENIOR(10,1,0,1)/ -1/  
!     l=10 #=0 +
      DATA JKVANT(10,1,1,1)/ 21/  
!     data  SENIOR(10,1,1,1)          / 1/
      DATA SENIOR(10,1,1,1)/ -1/  
!     l=10 #=1 +
      DATA (JKVANT(10,1,2,I),I=1,11)/ 0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40/  
!     data (SENIOR(10,1,2,i),i=1,11)   / 0, 2, 2, 2, 2, 2, 2, 2, 2, 2,
!    :                                   2/
      DATA (SENIOR(10,1,2,I),I=1,11)/ -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, &
         -1/  
!     l=10 #=2 +
      IF (FIRST) THEN 
         FIL = FIL_1 
      ELSE 
         FIL = FIL_2 
      ENDIF 
      ANTKO = 1 
      POS = 0 
      DO I = 1, 110 
         N = POSN(I) 
         L = POSL(I) 
!Jacek mailed the fix 98-10-29
         IF (N < 10) THEN 
            DO K = 0, MIN(L,1) 
         !do 20 k=0,min(n-1,1)
               IF (ANSATS(N,L,K) == 0) CYCLE  
               RAD1(POS*9+1:POS*9+9) = '         ' 
               RAD1(POS*9+3:POS*9+3) = CHAR(48 + N) 
               RAD1(POS*9+4:POS*9+5) = L1(L,K) 
               RAD1(POS*9+6:POS*9+9) = '(  )' 
               IF (ANSATS(N,L,K) >= 10) THEN 
                  RAD1(POS*9+7:POS*9+8) = CHAR(ANSATS(N,L,K)/10+48) 
               ELSE 
                  RAD1(POS*9+7:POS*9+7) = ' ' 
               ENDIF 
               RAD1(POS*9+8:POS*9+8) = CHAR(MOD(ANSATS(N,L,K),10)+48) 
               POS = POS + 1 
               IF (POS > SKAL) THEN 
                  WRITE (*, *) 'More than 20 subshells' 
                  RETURN  
               ENDIF 
               ORBIT(POS) = L 
               ANTEL(POS) = MIN(ANSATS(N,L,K),ANTMAX(L,K)-ANSATS(N,L,K)) 
               ANTKO(POS) = KOPPL(L,ANTEL(POS),K) 
               PLUS(POS) = K 
            END DO 
         ELSE 
            DO K = 0, MIN(L,1) 
         !do 20 k=0,min(n-1,1)
               IF (ANSATS(N,L,K) == 0) CYCLE  
               RAD1(POS*9+1:POS*9+9) = '         ' 
               N1 = MOD(N,10) 
               N10 = N/10 
               RAD1(POS*9+2:POS*9+2) = CHAR(48 + N10) 
               RAD1(POS*9+3:POS*9+3) = CHAR(48 + N1) 
               RAD1(POS*9+4:POS*9+5) = L1(L,K) 
               RAD1(POS*9+6:POS*9+9) = '(  )' 
               IF (ANSATS(N,L,K) >= 10) THEN 
                  RAD1(POS*9+7:POS*9+8) = CHAR(ANSATS(N,L,K)/10+48) 
               ELSE 
                  RAD1(POS*9+7:POS*9+7) = ' ' 
               ENDIF 
               RAD1(POS*9+8:POS*9+8) = CHAR(MOD(ANSATS(N,L,K),10)+48) 
               POS = POS + 1 
               IF (POS > SKAL) THEN 
                  WRITE (*, *) 'More than 20 subshells' 
                  RETURN  
               ENDIF 
               ORBIT(POS) = L 
               ANTEL(POS) = MIN(ANSATS(N,L,K),ANTMAX(L,K)-ANSATS(N,L,K)) 
               ANTKO(POS) = KOPPL(L,ANTEL(POS),K) 
               PLUS(POS) = K 
            END DO 
         ENDIF 
      END DO 
 
      IF (POS == 0) RETURN  
      DO I1 = 1, ANTKO(1) 
         DO I2 = 1, ANTKO(2) 
            DO I3 = 1, ANTKO(3) 
               DO I4 = 1, ANTKO(4) 
                  DO I5 = 1, ANTKO(5) 
                     DO I6 = 1, ANTKO(6) 
                        DO I7 = 1, ANTKO(7) 
                           DO I8 = 1, ANTKO(8) 
                              DO I9 = 1, ANTKO(9) 
                                 DO I10 = 1, ANTKO(10) 
                                    DO I11 = 1, ANTKO(11) 
                                    DO I12 = 1, ANTKO(12) 
                                    DO I13 = 1, ANTKO(13) 
                                    DO I14 = 1, ANTKO(14) 
                                    DO I15 = 1, ANTKO(15) 
                                    DO I16 = 1, ANTKO(16) 
                                    DO I17 = 1, ANTKO(17) 
                                    DO I18 = 1, ANTKO(18) 
                                    DO I19 = 1, ANTKO(19) 
                                    DO I20 = 1, ANTKO(20) 
 
                                    J(1) = JKVANT(ORBIT(1),PLUS(1),ANTEL(1),I1) 
                                    S(1) = SENIOR(ORBIT(1),PLUS(1),ANTEL(1),I1) 
                                    IF (POS == 1) THEN 
                                    IF (J(1)>=MINJ .AND. J(1)<=MAXJ) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, J, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:9) 
                                    WRITE (FIL, 999) RAD2(1:9) 
                                    WRITE (FIL, 999) RAD3(1:11) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    DO RESJ = MINJ, MAXJ, 2 
                                    JK(POS-1) = RESJ 
                                    J(2) = JKVANT(ORBIT(2),PLUS(2),ANTEL(2),I2) 
                                    S(2) = SENIOR(ORBIT(2),PLUS(2),ANTEL(2),I2) 
                                    IF (POS == 2) THEN 
                                    IF (RESJ>=ABS(J(1)-J(2)) .AND. RESJ<=J(1)+J&
                                       (2)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:18) 
                                    WRITE (FIL, 999) RAD2(1:18) 
                                    WRITE (FIL, 999) RAD3(1:20) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(3) = JKVANT(ORBIT(3),PLUS(3),ANTEL(3),I3) 
                                    S(3) = SENIOR(ORBIT(3),PLUS(3),ANTEL(3),I3) 
                                    DO JK1 = ABS(J(1)-J(2)), J(1) + J(2), 2 
                                    JK(1) = JK1 
                                    IF (POS == 3) THEN 
                                    IF (RESJ>=ABS(JK1 - J(3)) .AND. RESJ<=JK1+J&
                                       (3)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:27) 
                                    WRITE (FIL, 999) RAD2(1:27) 
                                    WRITE (FIL, 999) RAD3(1:29) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(4) = JKVANT(ORBIT(4),PLUS(4),ANTEL(4),I4) 
                                    S(4) = SENIOR(ORBIT(4),PLUS(4),ANTEL(4),I4) 
                                    DO JK2 = ABS(JK1 - J(3)), JK1 + J(3), 2 
                                    JK(2) = JK2 
                                    IF (POS == 4) THEN 
                                    IF (RESJ>=ABS(JK2 - J(4)) .AND. RESJ<=JK2+J&
                                       (4)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:36) 
                                    WRITE (FIL, 999) RAD2(1:36) 
                                    WRITE (FIL, 999) RAD3(1:38) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(5) = JKVANT(ORBIT(5),PLUS(5),ANTEL(5),I5) 
                                    S(5) = SENIOR(ORBIT(5),PLUS(5),ANTEL(5),I5) 
                                    DO JK3 = ABS(JK2 - J(4)), JK2 + J(4), 2 
                                    JK(3) = JK3 
                                    IF (POS == 5) THEN 
                                    IF (RESJ>=ABS(JK3 - J(5)) .AND. RESJ<=JK3+J&
                                       (5)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:45) 
                                    WRITE (FIL, 999) RAD2(1:45) 
                                    WRITE (FIL, 999) RAD3(1:47) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(6) = JKVANT(ORBIT(6),PLUS(6),ANTEL(6),I6) 
                                    S(6) = SENIOR(ORBIT(6),PLUS(6),ANTEL(6),I6) 
                                    DO JK4 = ABS(JK3 - J(5)), JK3 + J(5), 2 
                                    JK(4) = JK4 
                                    IF (POS == 6) THEN 
                                    IF (RESJ>=ABS(JK4 - J(6)) .AND. RESJ<=JK4+J&
                                       (6)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:54) 
                                    WRITE (FIL, 999) RAD2(1:54) 
                                    WRITE (FIL, 999) RAD3(1:56) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(7) = JKVANT(ORBIT(7),PLUS(7),ANTEL(7),I7) 
                                    S(7) = SENIOR(ORBIT(7),PLUS(7),ANTEL(7),I7) 
                                    DO JK5 = ABS(JK4 - J(6)), JK4 + J(6), 2 
                                    JK(5) = JK5 
                                    IF (POS == 7) THEN 
                                    IF (RESJ>=ABS(JK5 - J(7)) .AND. RESJ<=JK5+J&
                                       (7)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:63) 
                                    WRITE (FIL, 999) RAD2(1:63) 
                                    WRITE (FIL, 999) RAD3(1:65) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(8) = JKVANT(ORBIT(8),PLUS(8),ANTEL(8),I8) 
                                    S(8) = SENIOR(ORBIT(8),PLUS(8),ANTEL(8),I8) 
                                    DO JK6 = ABS(JK5 - J(7)), JK5 + J(7), 2 
                                    JK(6) = JK6 
                                    IF (POS == 8) THEN 
                                    IF (RESJ>=ABS(JK6 - J(8)) .AND. RESJ<=JK6+J&
                                       (8)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:72) 
                                    WRITE (FIL, 999) RAD2(1:72) 
                                    WRITE (FIL, 999) RAD3(1:74) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(9) = JKVANT(ORBIT(9),PLUS(9),ANTEL(9),I9) 
                                    S(9) = SENIOR(ORBIT(9),PLUS(9),ANTEL(9),I9) 
                                    DO JK7 = ABS(JK6 - J(8)), JK6 + J(8), 2 
                                    JK(7) = JK7 
                                    IF (POS == 9) THEN 
                                    IF (RESJ>=ABS(JK7 - J(9)) .AND. RESJ<=JK7+J&
                                       (9)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:81) 
                                    WRITE (FIL, 999) RAD2(1:81) 
                                    WRITE (FIL, 999) RAD3(1:83) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(10) = JKVANT(ORBIT(10),PLUS(10),ANTEL(10)&
                                       ,I10) 
                                    S(10) = SENIOR(ORBIT(10),PLUS(10),ANTEL(10)&
                                       ,I10) 
                                    DO JK8 = ABS(JK7 - J(9)), JK7 + J(9), 2 
                                    JK(8) = JK8 
                                    IF (POS == 10) THEN 
                                    IF (RESJ>=ABS(JK8 - J(10)) .AND. RESJ<=JK8+&
                                       J(10)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:90) 
                                    WRITE (FIL, 999) RAD2(1:90) 
                                    WRITE (FIL, 999) RAD3(1:92) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(11) = JKVANT(ORBIT(11),PLUS(11),ANTEL(11)&
                                       ,I11) 
                                    S(11) = SENIOR(ORBIT(11),PLUS(11),ANTEL(11)&
                                       ,I11) 
                                    DO JK9 = ABS(JK8 - J(10)), JK8 + J(10), 2 
                                    JK(9) = JK9 
                                    IF (POS == 11) THEN 
                                    IF (RESJ>=ABS(JK9 - J(11)) .AND. RESJ<=JK9+&
                                       J(11)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:99) 
                                    WRITE (FIL, 999) RAD2(1:99) 
                                    WRITE (FIL, 999) RAD3(1:101) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(12) = JKVANT(ORBIT(12),PLUS(12),ANTEL(12)&
                                       ,I12) 
                                    S(12) = SENIOR(ORBIT(12),PLUS(12),ANTEL(12)&
                                       ,I12) 
                                    DO JK10 = ABS(JK9 - J(11)), JK9 + J(11), 2 
                                    JK(10) = JK10 
                                    IF (POS == 12) THEN 
                                    IF (RESJ>=ABS(JK10 - J(12)) .AND. RESJ<=&
                                       JK10+J(12)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:108) 
                                    WRITE (FIL, 999) RAD2(1:108) 
                                    WRITE (FIL, 999) RAD3(1:110) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(13) = JKVANT(ORBIT(13),PLUS(13),ANTEL(13)&
                                       ,I13) 
                                    S(13) = SENIOR(ORBIT(13),PLUS(13),ANTEL(13)&
                                       ,I13) 
                                    DO JK11 = ABS(JK10 - J(12)), JK10 + J(12), &
                                       2 
                                    JK(11) = JK11 
                                    IF (POS == 13) THEN 
                                    IF (RESJ>=ABS(JK11 - J(13)) .AND. RESJ<=&
                                       JK11+J(13)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:117) 
                                    WRITE (FIL, 999) RAD2(1:117) 
                                    WRITE (FIL, 999) RAD3(1:119) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(14) = JKVANT(ORBIT(14),PLUS(14),ANTEL(14)&
                                       ,I14) 
                                    S(14) = SENIOR(ORBIT(14),PLUS(14),ANTEL(14)&
                                       ,I14) 
                                    DO JK12 = ABS(JK11 - J(13)), JK11 + J(13), &
                                       2 
                                    JK(12) = JK12 
                                    IF (POS == 14) THEN 
                                    IF (RESJ>=ABS(JK12 - J(14)) .AND. RESJ<=&
                                       JK12+J(14)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:126) 
                                    WRITE (FIL, 999) RAD2(1:126) 
                                    WRITE (FIL, 999) RAD3(1:128) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(15) = JKVANT(ORBIT(15),PLUS(15),ANTEL(15)&
                                       ,I15) 
                                    S(15) = SENIOR(ORBIT(15),PLUS(15),ANTEL(15)&
                                       ,I15) 
                                    DO JK13 = ABS(JK12 - J(14)), JK12 + J(14), &
                                       2 
                                    JK(13) = JK13 
                                    IF (POS == 15) THEN 
                                    IF (RESJ>=ABS(JK13 - J(15)) .AND. RESJ<=&
                                       JK13+J(15)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:135) 
                                    WRITE (FIL, 999) RAD2(1:135) 
                                    WRITE (FIL, 999) RAD3(1:137) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(16) = JKVANT(ORBIT(16),PLUS(16),ANTEL(16)&
                                       ,I16) 
                                    S(16) = SENIOR(ORBIT(16),PLUS(16),ANTEL(16)&
                                       ,I16) 
                                    DO JK14 = ABS(JK13 - J(15)), JK13 + J(15), &
                                       2 
                                    JK(14) = JK14 
                                    IF (POS == 16) THEN 
                                    IF (RESJ>=ABS(JK14 - J(16)) .AND. RESJ<=&
                                       JK14+J(16)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:144) 
                                    WRITE (FIL, 999) RAD2(1:144) 
                                    WRITE (FIL, 999) RAD3(1:146) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(17) = JKVANT(ORBIT(17),PLUS(17),ANTEL(17)&
                                       ,I17) 
                                    S(17) = SENIOR(ORBIT(17),PLUS(17),ANTEL(17)&
                                       ,I17) 
                                    DO JK15 = ABS(JK14 - J(16)), JK14 + J(16), &
                                       2 
                                    JK(15) = JK15 
                                    IF (POS == 17) THEN 
                                    IF (RESJ>=ABS(JK15 - J(17)) .AND. RESJ<=&
                                       JK15+J(17)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:153) 
                                    WRITE (FIL, 999) RAD2(1:153) 
                                    WRITE (FIL, 999) RAD3(1:155) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(18) = JKVANT(ORBIT(18),PLUS(18),ANTEL(18)&
                                       ,I18) 
                                    S(18) = SENIOR(ORBIT(18),PLUS(18),ANTEL(18)&
                                       ,I18) 
                                    DO JK16 = ABS(JK15 - J(17)), JK15 + J(17), &
                                       2 
                                    JK(16) = JK16 
                                    IF (POS == 18) THEN 
                                    IF (RESJ>=ABS(JK16 - J(18)) .AND. RESJ<=&
                                       JK16+J(18)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:162) 
                                    WRITE (FIL, 999) RAD2(1:162) 
                                    WRITE (FIL, 999) RAD3(1:164) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(19) = JKVANT(ORBIT(19),PLUS(19),ANTEL(19)&
                                       ,I19) 
                                    S(19) = SENIOR(ORBIT(19),PLUS(19),ANTEL(19)&
                                       ,I19) 
                                    DO JK17 = ABS(JK16 - J(18)), JK16 + J(18), &
                                       2 
                                    JK(17) = JK17 
                                    IF (POS == 19) THEN 
                                    IF (RESJ>=ABS(JK17 - J(19)) .AND. RESJ<=&
                                       JK17+J(19)) THEN 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:171) 
                                    WRITE (FIL, 999) RAD2(1:171) 
                                    WRITE (FIL, 999) RAD3(1:173) 
                                    CF = CF + 1 
                                    ENDIF 
                                    ELSE 
 
                                    J(20) = JKVANT(ORBIT(20),PLUS(20),ANTEL(20)&
                                       ,I20) 
                                    S(20) = SENIOR(ORBIT(20),PLUS(20),ANTEL(20)&
                                       ,I20) 
                                    DO JK18 = ABS(JK17 - J(19)), JK17 + J(19), &
                                       2 
                                    IF (RESJ<ABS(JK18 - J(20)) .OR. RESJ>JK18+J&
                                       (20)) CYCLE  
                                    JK(18) = JK18 
                                    CALL KOPP1 (POS, RAD2, J, S, ANTKO) 
                                    CALL KOPP2 (POS, RAD3, JK, J, PAR, ANTKO) 
                                    WRITE (FIL, 999) RAD1(1:180) 
                                    WRITE (FIL, 999) RAD2(1:180) 
                                    WRITE (FIL, 999) RAD3(1:182) 
                                    CF = CF + 1 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    ENDIF 
                                    END DO 
                                    END DO 
                                    END DO 
                                    END DO 
                                    END DO 
                                    END DO 
                                    END DO 
                                    END DO 
                                    END DO 
                                    END DO 
                                 END DO 
                              END DO 
                           END DO 
                        END DO 
                     END DO 
                  END DO 
               END DO 
            END DO 
         END DO 
      END DO 
  999 FORMAT(2A) 
      RETURN  
      END SUBROUTINE GEN 
