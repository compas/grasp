!     last edited Januar 2, 1997
      subroutine blanda(org, varmax, lock, minj, maxj, skal, nmax, low, posn, &
         posl, lim, dubbel, first)
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:44:54  12/27/06
!...Switches:
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use slug_I
      use gen_I
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: varmax
      integer  :: minj
      integer  :: maxj
      integer  :: skal
      integer , intent(in) :: nmax
      logical  :: first
      integer  :: org(15,0:10)
      integer  :: low(15,0:10)
      integer  :: posn(110)
      integer  :: posl(110)
      integer , intent(in) :: lim(15)
      logical  :: lock(15,0:10)
      logical  :: dubbel(15,0:10)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: fil_1 = 7
      integer, parameter :: fil_2 = 8
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(15,0:10) :: antel, start
      integer :: cf
      integer , dimension(15,0:10,0:1) :: ansats
      integer , dimension(15,0:10) :: varupp, varned
      integer :: an10, an20, an21, an30, an31, an32, an40, an41, an42, an43, k&
         , an50, an51, an52, an53, an54, an60, an61, an62, an63, an64, an65, &
         an70, an71, an72, an73, an74, an75, an76
      integer , dimension(15,0:10) :: stopp
      integer :: an80, an81, an82, an83, an84, an85, an86, an87, an90, an91, &
         an92, an93, an94, an95, an96, an97, an98, ana0, ana1, ana2, ana3, ana4&
         , ana5, ana6, ana7, ana8, ana9, plus21, plus31, plus32, plus41, plus42&
         , plus43, plus51, plus52, plus53, plus54, plus61, plus62, plus63, &
         plus64, plus65, plus71, plus72, plus73, plus74, plus75, plus76, plus81&
         , plus82, plus83, plus84, plus85, plus86, plus87, plus91, plus92, &
         plus93, plus94, plus95, plus96, plus97, plus98, plusa1, plusa2, plusa3&
         , plusa4, plusa5, plusa6, plusa7, plusa8, plusa9, par0, par, ress, &
         resl, i, j, antal
      integer , dimension(15,0:10) :: steg
      integer :: dum, ras1, ras3, elar, rasett, rastre
      integer , dimension(15,0:10) :: ras
      integer :: plusba, plusca, plusda, plusea, plusfa, plusb1, plusb2, plusb3&
         , plusb4, plusb5, plusb6, plusb7, plusb8, plusb9, plusc1, plusc2, &
         plusc3, plusc4, plusc5, plusc6, plusc7, plusc8, plusc9, plusd1, plusd2&
         , plusd3, plusd4, plusd5, plusd6, plusd7, plusd8, plusd9, pluse1, &
         pluse2, pluse3, pluse4, pluse5, pluse6, pluse7, pluse8, pluse9, plusf1&
         , plusf2, plusf3, plusf4, plusf5, plusf6, plusf7, plusf8, plusf9, anba&
         , anca, anda, anea, anfa, anb0, anb1, anb2, anb3, anb4, anb5, anb6, &
         anb7, anb8, anb9, anc0, anc1, anc2, anc3, anc4, anc5, anc6, anc7, anc8&
         , anc9, and0, and1, and2, and3, and4, and5, and6, and7, and8, and9, &
         ane0, ane1, ane2, ane3, ane4, ane5, ane6, ane7, ane8, ane9, anf0, anf1&
         , anf2, anf3, anf4, anf5, anf6, anf7, anf8, anf9
      logical :: finns, napp
!-----------------------------------------------

      cf = 0
      antal = 0
      par0 = 0
      finns = .FALSE.
      do i = 1, nmax
         do j = 0, min(10,i - 1)
            if (dubbel(i,j)) then
               steg(i,j) = -2
            else
               steg(i,j) = -1
            endif
            antal = antal + org(i,j)
            par0 = mod(par0 + j*org(i,j),2)
         end do
      end do
      if (nmax < 15) then
         do i = nmax + 1, 15
            steg(i,:min(10,i-1)) = -1
         end do
      endif
!     1s
      call slug (1, 0, varmax, varupp, varned, ansats, org, lock(1,0), dubbel, &
         low, start(1,0), stopp(1,0))
      do an10 = start(1,0), stopp(1,0), steg(1,0)
         antel(1,0) = an10
         if (antel(1,0)>antal .or. antel(1,0)<lim(1)) cycle
         ansats(1,0,0) = an10
!     2s
         call slug (2, 0, varmax, varupp, varned, ansats, org, lock(2,0), &
            dubbel, low, start(2,0), stopp(2,0))
         do an20 = start(2,0), stopp(2,0), steg(2,0)
            antel(2,0) = an20 + antel(1,0)
            if (antel(2,0) > antal) cycle
            ansats(2,0,0) = an20
!     2p
            call slug (2, 1, varmax, varupp, varned, ansats, org, lock(2,1), &
               dubbel, low, start(2,1), stopp(2,1))
            do an21 = start(2,1), stopp(2,1), steg(2,1)
               antel(2,1) = an21 + antel(2,0)
               if (antel(2,1)>antal .or. antel(2,1)<lim(2)) cycle
               do plus21 = min(an21,4), max(an21 - 2,0), -1
                  ansats(2,1,1) = plus21
                  ansats(2,1,0) = an21 - plus21
!     3s
                  call slug (3, 0, varmax, varupp, varned, ansats, org, lock(3,&
                     0), dubbel, low, start(3,0), stopp(3,0))
                  do an30 = start(3,0), stopp(3,0), steg(3,0)
                     antel(3,0) = an30 + antel(2,1)
                     if (antel(3,0) > antal) cycle
                     ansats(3,0,0) = an30
!     3p
                     call slug (3, 1, varmax, varupp, varned, ansats, org, lock&
                        (3,1), dubbel, low, start(3,1), stopp(3,1))
                     do an31 = start(3,1), stopp(3,1), steg(3,1)
                        antel(3,1) = an31 + antel(3,0)
                        if (antel(3,1) > antal) cycle
                        do plus31 = min(an31,4), max(an31 - 2,0), -1
                           ansats(3,1,1) = plus31
                           ansats(3,1,0) = an31 - plus31
!     3d
                           call slug (3, 2, varmax, varupp, varned, ansats, org&
                              , lock(3,2), dubbel, low, start(3,2), stopp(3,2))
                           do an32 = start(3,2), stopp(3,2), steg(3,2)
                              antel(3,2) = an32 + antel(3,1)
                              if (antel(3,2)>antal .or. antel(3,2)<lim(3)) &
                                 cycle
                              do plus32 = min(an32,6), max(an32 - 4,0), -1
                                 ansats(3,2,1) = plus32
                                 ansats(3,2,0) = an32 - plus32
!     4s
                                 call slug (4, 0, varmax, varupp, varned, &
                                    ansats, org, lock(4,0), dubbel, low, start(&
                                    4,0), stopp(4,0))
                                 do an40 = start(4,0), stopp(4,0), steg(4,0)
                                    antel(4,0) = an40 + antel(3,2)
                                    if (antel(4,0) > antal) cycle
                                    ansats(4,0,0) = an40
!     4p
                                    call slug (4, 1, varmax, varupp, varned, &
                                       ansats, org, lock(4,1), dubbel, low, &
                                       start(4,1), stopp(4,1))
                                    do an41 = start(4,1), stopp(4,1), steg(4,1)
                                    antel(4,1) = an41 + antel(4,0)
                                    if (antel(4,1) > antal) cycle
                                    do plus41 = min(an41,4), max(an41 - 2,0), &
                                       -1
                                    ansats(4,1,1) = plus41
                                    ansats(4,1,0) = an41 - plus41
!     4d
                                    call slug (4, 2, varmax, varupp, varned, &
                                       ansats, org, lock(4,2), dubbel, low, &
                                       start(4,2), stopp(4,2))
                                    do an42 = start(4,2), stopp(4,2), steg(4,2)
                                    antel(4,2) = an42 + antel(4,1)
                                    if (antel(4,2) > antal) cycle
                                    do plus42 = min(an42,6), max(an42 - 4,0), &
                                       -1
                                    ansats(4,2,1) = plus42
                                    ansats(4,2,0) = an42 - plus42
!     4f
                                    call slug (4, 3, varmax, varupp, varned, &
                                       ansats, org, lock(4,3), dubbel, low, &
                                       start(4,3), stopp(4,3))
                                    do an43 = start(4,3), stopp(4,3), steg(4,3)
                                    antel(4,3) = an43 + antel(4,2)
                                    if (antel(4,3)>antal .or. antel(4,3)<lim(4)&
                                       ) cycle
                                    do plus43 = min(an43,8), max(an43 - 6,0), &
                                       -1
                                    ansats(4,3,1) = plus43
                                    ansats(4,3,0) = an43 - plus43
!     5s
                                    call slug (5, 0, varmax, varupp, varned, &
                                       ansats, org, lock(5,0), dubbel, low, &
                                       start(5,0), stopp(5,0))
                                    do an50 = start(5,0), stopp(5,0), steg(5,0)
                                    antel(5,0) = an50 + antel(4,3)
                                    if (antel(5,0) > antal) cycle
                                    ansats(5,0,0) = an50
!     5p
                                    call slug (5, 1, varmax, varupp, varned, &
                                       ansats, org, lock(5,1), dubbel, low, &
                                       start(5,1), stopp(5,1))
                                    do an51 = start(5,1), stopp(5,1), steg(5,1)
                                    antel(5,1) = an51 + antel(5,0)
                                    if (antel(5,1) > antal) cycle
                                    do plus51 = min(an51,4), max(an51 - 2,0), &
                                       -1
                                    ansats(5,1,1) = plus51
                                    ansats(5,1,0) = an51 - plus51
!     5d
                                    call slug (5, 2, varmax, varupp, varned, &
                                       ansats, org, lock(5,2), dubbel, low, &
                                       start(5,2), stopp(5,2))
                                    do an52 = start(5,2), stopp(5,2), steg(5,2)
                                    antel(5,2) = an52 + antel(5,1)
                                    if (antel(5,2) > antal) cycle
                                    do plus52 = min(an52,6), max(an52 - 4,0), &
                                       -1
                                    ansats(5,2,1) = plus52
                                    ansats(5,2,0) = an52 - plus52

!     5f
                                    call slug (5, 3, varmax, varupp, varned, &
                                       ansats, org, lock(5,3), dubbel, low, &
                                       start(5,3), stopp(5,3))
                                    do an53 = start(5,3), stopp(5,3), steg(5,3)
                                    antel(5,3) = an53 + antel(5,2)
                                    if (antel(5,3) > antal) cycle
                                    do plus53 = min(an53,8), max(an53 - 6,0), &
                                       -1
                                    ansats(5,3,1) = plus53
                                    ansats(5,3,0) = an53 - plus53
!     5g
                                    call slug (5, 4, varmax, varupp, varned, &
                                       ansats, org, lock(5,4), dubbel, low, &
                                       start(5,4), stopp(5,4))
                                    do an54 = start(5,4), stopp(5,4), steg(5,4)
                                    antel(5,4) = an54 + antel(5,3)
                                    if (antel(5,4)>antal .or. antel(5,4)<lim(5)&
                                       ) cycle
                                    do plus54 = min(an54,10), max(an54 - 8,0), &
                                       -1
                                    ansats(5,4,1) = plus54
                                    ansats(5,4,0) = an54 - plus54
!     6s
                                    call slug (6, 0, varmax, varupp, varned, &
                                       ansats, org, lock(6,0), dubbel, low, &
                                       start(6,0), stopp(6,0))
                                    do an60 = start(6,0), stopp(6,0), steg(6,0)
                                    antel(6,0) = an60 + antel(5,4)
                                    if (antel(6,0)>antal .or. ansats(5,4,1)>2) &
                                       cycle
                                    ansats(6,0,0) = an60
!     6p
                                    call slug (6, 1, varmax, varupp, varned, &
                                       ansats, org, lock(6,1), dubbel, low, &
                                       start(6,1), stopp(6,1))
                                    do an61 = start(6,1), stopp(6,1), steg(6,1)
                                    antel(6,1) = an61 + antel(6,0)
                                    if (antel(6,1) > antal) cycle
                                    do plus61 = min(an61,4), max(an61 - 2,0), &
                                       -1
                                    ansats(6,1,1) = plus61
                                    ansats(6,1,0) = an61 - plus61
!     6d
                                    call slug (6, 2, varmax, varupp, varned, &
                                       ansats, org, lock(6,2), dubbel, low, &
                                       start(6,2), stopp(6,2))
                                    do an62 = start(6,2), stopp(6,2), steg(6,2)
                                    antel(6,2) = an62 + antel(6,1)
                                    if (antel(6,2) > antal) cycle
                                    do plus62 = min(an62,6), max(an62 - 4,0), &
                                       -1
                                    ansats(6,2,1) = plus62
                                    ansats(6,2,0) = an62 - plus62
!     6f
                                    call slug (6, 3, varmax, varupp, varned, &
                                       ansats, org, lock(6,3), dubbel, low, &
                                       start(6,3), stopp(6,3))
                                    do an63 = start(6,3), stopp(6,3), steg(6,3)
                                    antel(6,3) = an63 + antel(6,2)
                                    if (antel(6,3) > antal) cycle
                                    do plus63 = min(an63,8), max(an63 - 6,0), &
                                       -1
                                    ansats(6,3,1) = plus63
                                    ansats(6,3,0) = an63 - plus63
!     6g
                                    call slug (6, 4, varmax, varupp, varned, &
                                       ansats, org, lock(6,4), dubbel, low, &
                                       start(6,4), stopp(6,4))
                                    do an64 = start(6,4), stopp(6,4), steg(6,4)
                                    antel(6,4) = an64 + antel(6,3)
                                    if (antel(6,4) > antal) cycle
                                    do plus64 = min(an64,10), max(an64 - 8,0), &
                                       -1
                                    ansats(6,4,1) = plus64
                                    ansats(6,4,0) = an64 - plus64
!     6h
                                    call slug (6, 5, varmax, varupp, varned, &
                                       ansats, org, lock(6,5), dubbel, low, &
                                       start(6,5), stopp(6,5))
                                    do an65 = start(6,5), stopp(6,5), steg(6,5)
                                    antel(6,5) = an65 + antel(6,4)
                                    if (.not.(antel(6,5)<=antal .and. ansats(6,&
                                       4,1)<=2 .and. antel(6,5)>=lim(6))) &
                                       cycle
                                    do plus65 = min(an65,12), max(an65 - 10,0)&
                                       , -1
                                    ansats(6,5,1) = plus65
                                    ansats(6,5,0) = an65 - plus65
!     7s
                                    call slug (7, 0, varmax, varupp, varned, &
                                       ansats, org, lock(7,0), dubbel, low, &
                                       start(7,0), stopp(7,0))
                                    do an70 = start(7,0), stopp(7,0), steg(7,0)
                                    antel(7,0) = an70 + antel(6,5)
                                    if (.not.(antel(7,0)<=antal .and. ansats(6,&
                                       5,1)<=2 .and. ansats(6,5,0)<=2)) cycle
                                    ansats(7,0,0) = an70
!     7p
                                    call slug (7, 1, varmax, varupp, varned, &
                                       ansats, org, lock(7,1), dubbel, low, &
                                       start(7,1), stopp(7,1))
                                    do an71 = start(7,1), stopp(7,1), steg(7,1)
                                    antel(7,1) = an71 + antel(7,0)
                                    if (antel(7,1) > antal) cycle
                                    do plus71 = min(an71,4), max(an71 - 2,0), &
                                       -1
                                    ansats(7,1,1) = plus71
                                    ansats(7,1,0) = an71 - plus71
!     7d
                                    call slug (7, 2, varmax, varupp, varned, &
                                       ansats, org, lock(7,2), dubbel, low, &
                                       start(7,2), stopp(7,2))
                                    do an72 = start(7,2), stopp(7,2), steg(7,2)
                                    antel(7,2) = an72 + antel(7,1)
                                    if (antel(7,2) > antal) cycle
                                    do plus72 = min(an72,6), max(an72 - 4,0), &
                                       -1
                                    ansats(7,2,1) = plus72
                                    ansats(7,2,0) = an72 - plus72
!     7f
                                    call slug (7, 3, varmax, varupp, varned, &
                                       ansats, org, lock(7,3), dubbel, low, &
                                       start(7,3), stopp(7,3))
                                    do an73 = start(7,3), stopp(7,3), steg(7,3)
                                    antel(7,3) = an73 + antel(7,2)
                                    if (antel(7,3) > antal) cycle
                                    do plus73 = min(an73,8), max(an73 - 6,0), &
                                       -1
                                    ansats(7,3,1) = plus73
                                    ansats(7,3,0) = an73 - plus73
!     7g
                                    call slug (7, 4, varmax, varupp, varned, &
                                       ansats, org, lock(7,4), dubbel, low, &
                                       start(7,4), stopp(7,4))
                                    do an74 = start(7,4), stopp(7,4), steg(7,4)
                                    antel(7,4) = an74 + antel(7,3)
                                    if (antel(7,4) > antal) cycle
                                    do plus74 = min(an74,10), max(an74 - 8,0), &
                                       -1
                                    ansats(7,4,1) = plus74
                                    ansats(7,4,0) = an74 - plus74
!     7h
                                    call slug (7, 5, varmax, varupp, varned, &
                                       ansats, org, lock(7,5), dubbel, low, &
                                       start(7,5), stopp(7,5))
                                    do an75 = start(7,5), stopp(7,5), steg(7,5)
                                    antel(7,5) = an75 + antel(7,4)
                                    if (antel(7,5)>antal .or. ansats(7,4,1)>2) &
                                       cycle
                                    do plus75 = min(an75,12), max(an75 - 10,0)&
                                       , -1
                                    ansats(7,5,1) = plus75
                                    ansats(7,5,0) = an75 - plus75
!     7i
                                    call slug (7, 6, varmax, varupp, varned, &
                                       ansats, org, lock(7,6), dubbel, low, &
                                       start(7,6), stopp(7,6))
                                    do an76 = start(7,6), stopp(7,6), steg(7,6)
                                    antel(7,6) = an76 + antel(7,5)
                                    if (.not.(antel(7,6)<=antal .and. ansats(7,&
                                       5,1)<=2 .and. ansats(7,5,0)<=2 .and. &
                                       antel(7,6)>=lim(7))) cycle
                                    do plus76 = min(an76,14), max(an76 - 12,0)&
                                       , -1
                                    ansats(7,6,1) = plus76
                                    ansats(7,6,0) = an76 - plus76
!     8s
                                    call slug (8, 0, varmax, varupp, varned, &
                                       ansats, org, lock(8,0), dubbel, low, &
                                       start(8,0), stopp(8,0))
                                    do an80 = start(8,0), stopp(8,0), steg(8,0)
                                    antel(8,0) = an80 + antel(7,6)
                                    if (.not.(antel(8,0)<=antal .and. ansats(7,&
                                       6,1)<=2 .and. ansats(7,6,0)<=2)) cycle
                                    ansats(8,0,0) = an80
!     8p
                                    call slug (8, 1, varmax, varupp, varned, &
                                       ansats, org, lock(8,1), dubbel, low, &
                                       start(8,1), stopp(8,1))
                                    do an81 = start(8,1), stopp(8,1), steg(8,1)
                                    antel(8,1) = an81 + antel(8,0)
                                    if (antel(8,1) > antal) cycle
                                    do plus81 = min(an81,4), max(an81 - 2,0), &
                                       -1
                                    ansats(8,1,1) = plus81
                                    ansats(8,1,0) = an81 - plus81
!     8d
                                    call slug (8, 2, varmax, varupp, varned, &
                                       ansats, org, lock(8,2), dubbel, low, &
                                       start(8,2), stopp(8,2))
                                    do an82 = start(8,2), stopp(8,2), steg(8,2)
                                    antel(8,2) = an82 + antel(8,1)
                                    if (antel(8,2) > antal) cycle
                                    do plus82 = min(an82,6), max(an82 - 4,0), &
                                       -1
                                    ansats(8,2,1) = plus82
                                    ansats(8,2,0) = an82 - plus82
!     8f
                                    call slug (8, 3, varmax, varupp, varned, &
                                       ansats, org, lock(8,3), dubbel, low, &
                                       start(8,3), stopp(8,3))
                                    do an83 = start(8,3), stopp(8,3), steg(8,3)
                                    antel(8,3) = an83 + antel(8,2)
                                    if (antel(8,3) > antal) cycle
                                    do plus83 = min(an83,8), max(an83 - 6,0), &
                                       -1
                                    ansats(8,3,1) = plus83
                                    ansats(8,3,0) = an83 - plus83
!     8g
                                    call slug (8, 4, varmax, varupp, varned, &
                                       ansats, org, lock(8,4), dubbel, low, &
                                       start(8,4), stopp(8,4))

                                    do an84 = start(8,4), stopp(8,4), steg(8,4)
                                    antel(8,4) = an84 + antel(8,3)
                                    if (antel(8,4) > antal) cycle
                                    do plus84 = min(an84,10), max(an84 - 8,0), &
                                       -1
                                    ansats(8,4,1) = plus84
                                    ansats(8,4,0) = an84 - plus84
!     8h
                                    call slug (8, 5, varmax, varupp, varned, &
                                       ansats, org, lock(8,5), dubbel, low, &
                                       start(8,5), stopp(8,5))
                                    do an85 = start(8,5), stopp(8,5), steg(8,5)
                                    antel(8,5) = an85 + antel(8,4)
                                    if (antel(8,5)>antal .or. ansats(8,4,1)>2) &
                                       cycle
                                    do plus85 = min(an85,12), max(an85 - 10,0)&
                                       , -1
                                    ansats(8,5,1) = plus85
                                    ansats(8,5,0) = an85 - plus85
!     8i
                                    call slug (8, 6, varmax, varupp, varned, &
                                       ansats, org, lock(8,6), dubbel, low, &
                                       start(8,6), stopp(8,6))
                                    do an86 = start(8,6), stopp(8,6), steg(8,6)
                                    antel(8,6) = an86 + antel(8,5)
                                    if (.not.(antel(8,6)<=antal .and. ansats(8,&
                                       5,1)<=2 .and. ansats(8,5,0)<=2)) cycle
                                    do plus86 = min(an86,14), max(an86 - 12,0)&
                                       , -1
                                    ansats(8,6,1) = plus86
                                    ansats(8,6,0) = an86 - plus86
!     8k
                                    call slug (8, 7, varmax, varupp, varned, &
                                       ansats, org, lock(8,7), dubbel, low, &
                                       start(8,7), stopp(8,7))
                                    do an87 = start(8,7), stopp(8,7), steg(8,7)
                                    antel(8,7) = an87 + antel(8,6)
                                    if (.not.(antel(8,7)<=antal .and. ansats(8,&
                                       6,1)<=2 .and. ansats(8,6,0)<=2 .and. &
                                       antel(8,7)>=lim(8))) cycle
                                    do plus87 = min(an87,16), max(an87 - 14,0)&
                                       , -1
                                    ansats(8,7,1) = plus87
                                    ansats(8,7,0) = an87 - plus87
!     9s
                                    call slug (9, 0, varmax, varupp, varned, &
                                       ansats, org, lock(9,0), dubbel, low, &
                                       start(9,0), stopp(9,0))
                                    do an90 = start(9,0), stopp(9,0), steg(9,0)
                                    antel(9,0) = an90 + antel(8,7)
                                    if (.not.(antel(9,0)<=antal .and. ansats(8,&
                                       7,1)<=2 .and. ansats(8,7,0)<=2)) cycle
                                    ansats(9,0,0) = an90
!     9p
                                    call slug (9, 1, varmax, varupp, varned, &
                                       ansats, org, lock(9,1), dubbel, low, &
                                       start(9,1), stopp(9,1))
                                    do an91 = start(9,1), stopp(9,1), steg(9,1)
                                    antel(9,1) = an91 + antel(9,0)
                                    if (antel(9,1) > antal) cycle
                                    do plus91 = min(an91,4), max(an91 - 2,0), &
                                       -1
                                    ansats(9,1,1) = plus91
                                    ansats(9,1,0) = an91 - plus91
!     9d
                                    call slug (9, 2, varmax, varupp, varned, &
                                       ansats, org, lock(9,2), dubbel, low, &
                                       start(9,2), stopp(9,2))
                                    do an92 = start(9,2), stopp(9,2), steg(9,2)
                                    antel(9,2) = an92 + antel(9,1)
                                    if (antel(9,2) > antal) cycle
                                    do plus92 = min(an92,6), max(an92 - 4,0), &
                                       -1
                                    ansats(9,2,1) = plus92
                                    ansats(9,2,0) = an92 - plus92
!     9f
                                    call slug (9, 3, varmax, varupp, varned, &
                                       ansats, org, lock(9,3), dubbel, low, &
                                       start(9,3), stopp(9,3))
                                    do an93 = start(9,3), stopp(9,3), steg(9,3)
                                    antel(9,3) = an93 + antel(9,2)
                                    if (antel(9,3) > antal) cycle
                                    do plus93 = min(an93,8), max(an93 - 6,0), &
                                       -1
                                    ansats(9,3,1) = plus93
                                    ansats(9,3,0) = an93 - plus93
!     9g
                                    call slug (9, 4, varmax, varupp, varned, &
                                       ansats, org, lock(9,4), dubbel, low, &
                                       start(9,4), stopp(9,4))
                                    do an94 = start(9,4), stopp(9,4), steg(9,4)
                                    antel(9,4) = an94 + antel(9,3)
                                    if (antel(9,4) > antal) cycle
                                    do plus94 = min(an94,10), max(an94 - 8,0), &
                                       -1
                                    ansats(9,4,1) = plus94
                                    ansats(9,4,0) = an94 - plus94
!     9h
                                    call slug (9, 5, varmax, varupp, varned, &
                                       ansats, org, lock(9,5), dubbel, low, &
                                       start(9,5), stopp(9,5))
                                    do an95 = start(9,5), stopp(9,5), steg(9,5)
                                    antel(9,5) = an95 + antel(9,4)
                                    if (antel(9,5)>antal .or. ansats(9,4,1)>2) &
                                       cycle
                                    do plus95 = min(an95,12), max(an95 - 10,0)&
                                       , -1
                                    ansats(9,5,1) = plus95
                                    ansats(9,5,0) = an95 - plus95
!     9i
                                    call slug (9, 6, varmax, varupp, varned, &
                                       ansats, org, lock(9,6), dubbel, low, &
                                       start(9,6), stopp(9,6))
                                    do an96 = start(9,6), stopp(9,6), steg(9,6)
                                    antel(9,6) = an96 + antel(9,5)
                                    if (.not.(antel(9,6)<=antal .and. ansats(9,&
                                       5,1)<=2 .and. ansats(9,5,0)<=2)) cycle
                                    do plus96 = min(an96,14), max(an96 - 12,0)&
                                       , -1
                                    ansats(9,6,1) = plus96
                                    ansats(9,6,0) = an96 - plus96
!     9k
                                    call slug (9, 7, varmax, varupp, varned, &
                                       ansats, org, lock(9,7), dubbel, low, &
                                       start(9,7), stopp(9,7))
                                    do an97 = start(9,7), stopp(9,7), steg(9,7)
                                    antel(9,7) = an97 + antel(9,6)
                                    if (.not.(antel(9,7)<=antal .and. ansats(9,&
                                       6,1)<=2 .and. ansats(9,6,0)<=2)) cycle
                                    do plus97 = min(an97,16), max(an97 - 14,0)&
                                       , -1
                                    ansats(9,7,1) = plus97
                                    ansats(9,7,0) = an97 - plus97
!     9l
                                    call slug (9, 8, varmax, varupp, varned, &
                                       ansats, org, lock(9,8), dubbel, low, &
                                       start(9,8), stopp(9,8))
                                    do an98 = start(9,8), stopp(9,8), steg(9,8)
                                    antel(9,8) = an98 + antel(9,7)
                                    if (.not.(antel(9,8)<=antal .and. ansats(9,&
                                       7,1)<=2 .and. ansats(9,7,0)<=2 .and. &
                                       antel(9,8)>=lim(9))) cycle
                                    do plus98 = min(an98,18), max(an98 - 16,0)&
                                       , -1
                                    ansats(9,8,1) = plus98
                                    ansats(9,8,0) = an98 - plus98
!     10s
                                    call slug (10, 0, varmax, varupp, varned, &
                                       ansats, org, lock(10,0), dubbel, low, &
                                       start(10,0), stopp(10,0))
                                    do ana0 = start(10,0), stopp(10,0), steg(10&
                                       ,0)
                                    antel(10,0) = ana0 + antel(9,8)
                                    if (.not.(antel(10,0)<=antal .and. ansats(9&
                                       ,8,1)<=2 .and. ansats(9,8,0)<=2)) cycle
                                    ansats(10,0,0) = ana0
!     10p
                                    call slug (10, 1, varmax, varupp, varned, &
                                       ansats, org, lock(10,1), dubbel, low, &
                                       start(10,1), stopp(10,1))
                                    do ana1 = start(10,1), stopp(10,1), steg(10&
                                       ,1)
                                    antel(10,1) = ana1 + antel(10,0)
                                    if (antel(10,1) > antal) cycle
                                    do plusa1 = min(ana1,4), max(ana1 - 2,0), &
                                       -1
                                    ansats(10,1,1) = plusa1
                                    ansats(10,1,0) = ana1 - plusa1
!     10d
                                    call slug (10, 2, varmax, varupp, varned, &
                                       ansats, org, lock(10,2), dubbel, low, &
                                       start(10,2), stopp(10,2))
                                    do ana2 = start(10,2), stopp(10,2), steg(10&
                                       ,2)
                                    antel(10,2) = ana2 + antel(10,1)
                                    if (antel(10,2) > antal) cycle
                                    do plusa2 = min(ana2,6), max(ana2 - 4,0), &
                                       -1
                                    ansats(10,2,1) = plusa2
                                    ansats(10,2,0) = ana2 - plusa2
!     10f
                                    call slug (10, 3, varmax, varupp, varned, &
                                       ansats, org, lock(10,3), dubbel, low, &
                                       start(10,3), stopp(10,3))
                                    do ana3 = start(10,3), stopp(10,3), steg(10&
                                       ,3)
                                    antel(10,3) = ana3 + antel(10,2)
                                    if (antel(10,3) > antal) cycle
                                    do plusa3 = min(ana3,8), max(ana3 - 6,0), &
                                       -1
                                    ansats(10,3,1) = plusa3
                                    ansats(10,3,0) = ana3 - plusa3
!     10g
                                    call slug (10, 4, varmax, varupp, varned, &
                                       ansats, org, lock(10,4), dubbel, low, &
                                       start(10,4), stopp(10,4))
                                    do ana4 = start(10,4), stopp(10,4), steg(10&
                                       ,4)
                                    antel(10,4) = ana4 + antel(10,3)
                                    if (antel(10,4) > antal) cycle
                                    do plusa4 = min(ana4,10), max(ana4 - 8,0), &
                                       -1
                                    ansats(10,4,1) = plusa4
                                    ansats(10,4,0) = ana4 - plusa4
!     10h
                                    call slug (10, 5, varmax, varupp, varned, &
                                       ansats, org, lock(10,5), dubbel, low, &
                                       start(10,5), stopp(10,5))
                                    do ana5 = start(10,5), stopp(10,5), steg(10&
                                       ,5)
                                    antel(10,5) = ana5 + antel(10,4)
                                    if (antel(10,5)>antal .or. ansats(10,4,1)>2&
                                       ) cycle
                                    do plusa5 = min(ana5,12), max(ana5 - 10,0)&
                                       , -1
                                    ansats(10,5,1) = plusa5
                                    ansats(10,5,0) = ana5 - plusa5
!     10i
                                    call slug (10, 6, varmax, varupp, varned, &
                                       ansats, org, lock(10,6), dubbel, low, &
                                       start(10,6), stopp(10,6))
                                    do ana6 = start(10,6), stopp(10,6), steg(10&
                                       ,6)
                                    antel(10,6) = ana6 + antel(10,5)
                                    if (.not.(antel(10,6)<=antal .and. ansats(&
                                       10,5,1)<=2 .and. ansats(10,5,0)<=2)) &
                                       cycle
                                    do plusa6 = min(ana6,14), max(ana6 - 12,0)&
                                       , -1
                                    ansats(10,6,1) = plusa6
                                    ansats(10,6,0) = ana6 - plusa6
!     10k
                                    call slug (10, 7, varmax, varupp, varned, &
                                       ansats, org, lock(10,7), dubbel, low, &
                                       start(10,7), stopp(10,7))
                                    do ana7 = start(10,7), stopp(10,7), steg(10&
                                       ,7)
                                    antel(10,7) = ana7 + antel(10,6)
                                    if (.not.(antel(10,7)<=antal .and. ansats(&
                                       10,6,1)<=2 .and. ansats(10,6,0)<=2)) &
                                       cycle
                                    do plusa7 = min(ana7,16), max(ana7 - 14,0)&
                                       , -1
                                    ansats(10,7,1) = plusa7
                                    ansats(10,7,0) = ana7 - plusa7
!     10l
                                    call slug (10, 8, varmax, varupp, varned, &
                                       ansats, org, lock(10,8), dubbel, low, &
                                       start(10,8), stopp(10,8))
                                    do ana8 = start(10,8), stopp(10,8), steg(10&
                                       ,8)
                                    antel(10,8) = ana8 + antel(10,7)
                                    if (.not.(antel(10,8)<=antal .and. ansats(&
                                       10,7,1)<=2 .and. ansats(10,7,0)<=2)) &
                                       cycle
                                    do plusa8 = min(ana8,18), max(ana8 - 16,0)&
                                       , -1
                                    ansats(10,8,1) = plusa8
                                    ansats(10,8,0) = ana8 - plusa8
!     10m
                                    call slug (10, 9, varmax, varupp, varned, &
                                       ansats, org, lock(10,9), dubbel, low, &
                                       start(10,9), stopp(10,9))
                                    do ana9 = start(10,9), stopp(10,9), steg(10&
                                       ,9)
                                    antel(10,9) = ana9 + antel(10,8)
                                    if (.not.(antel(10,9)<=antal .and. ansats(&
                                       10,8,1)<=2 .and. ansats(10,8,0)<=2&
                                        .and. antel(10,9)>=lim(10))) cycle
                                    do plusa9 = min(ana9,20), max(ana9 - 18,0)&
                                       , -1
                                    ansats(10,9,1) = plusa9
                                    ansats(10,9,0) = ana9 - plusa9
!     11s
                                    call slug (11, 0, varmax, varupp, varned, &
                                       ansats, org, lock(11,0), dubbel, low, &
                                       start(11,0), stopp(11,0))
                                    do anb0 = start(11,0), stopp(11,0), steg(11&
                                       ,0)
                                    antel(11,0) = anb0 + antel(10,9)
                                    if (.not.(antel(11,0)<=antal .and. ansats(&
                                       10,9,1)<=2 .and. ansats(10,9,0)<=2)) &
                                       cycle
                                    ansats(11,0,0) = anb0
!     11p
                                    call slug (11, 1, varmax, varupp, varned, &
                                       ansats, org, lock(11,1), dubbel, low, &
                                       start(11,1), stopp(11,1))
                                    do anb1 = start(11,1), stopp(11,1), steg(11&
                                       ,1)
                                    antel(11,1) = anb1 + antel(11,0)
                                    if (antel(11,1) > antal) cycle
                                    do plusb1 = min(anb1,4), max(anb1 - 2,0), &
                                       -1
                                    ansats(11,1,1) = plusb1
                                    ansats(11,1,0) = anb1 - plusb1
!     11d
                                    call slug (11, 2, varmax, varupp, varned, &
                                       ansats, org, lock(11,2), dubbel, low, &
                                       start(11,2), stopp(11,2))
                                    do anb2 = start(11,2), stopp(11,2), steg(11&
                                       ,2)
                                    antel(11,2) = anb2 + antel(11,1)
                                    if (antel(11,2) > antal) cycle
                                    do plusb2 = min(anb2,6), max(anb2 - 4,0), &
                                       -1
                                    ansats(11,2,1) = plusb2
                                    ansats(11,2,0) = anb2 - plusb2
!     11f
                                    call slug (11, 3, varmax, varupp, varned, &
                                       ansats, org, lock(11,3), dubbel, low, &
                                       start(11,3), stopp(11,3))
                                    do anb3 = start(11,3), stopp(11,3), steg(11&
                                       ,3)
                                    antel(11,3) = anb3 + antel(11,2)
                                    if (antel(11,3) > antal) cycle
                                    do plusb3 = min(anb3,8), max(anb3 - 6,0), &
                                       -1
                                    ansats(11,3,1) = plusb3
                                    ansats(11,3,0) = anb3 - plusb3
!     11g
                                    call slug (11, 4, varmax, varupp, varned, &
                                       ansats, org, lock(11,4), dubbel, low, &
                                       start(11,4), stopp(11,4))
                                    do anb4 = start(11,4), stopp(11,4), steg(11&
                                       ,4)
                                    antel(11,4) = anb4 + antel(11,3)
                                    if (antel(11,4) > antal) cycle
                                    do plusb4 = min(anb4,10), max(anb4 - 8,0), &
                                       -1
                                    ansats(11,4,1) = plusb4
                                    ansats(11,4,0) = anb4 - plusb4
!     11h
                                    call slug (11, 5, varmax, varupp, varned, &
                                       ansats, org, lock(11,5), dubbel, low, &
                                       start(11,5), stopp(11,5))
                                    do anb5 = start(11,5), stopp(11,5), steg(11&
                                       ,5)
                                    antel(11,5) = anb5 + antel(11,4)
                                    if (antel(11,5)>antal .or. ansats(11,4,1)>2&
                                       ) cycle
                                    do plusb5 = min(anb5,12), max(anb5 - 10,0)&
                                       , -1
                                    ansats(11,5,1) = plusb5
                                    ansats(11,5,0) = anb5 - plusb5
!     11i
                                    call slug (11, 6, varmax, varupp, varned, &
                                       ansats, org, lock(11,6), dubbel, low, &
                                       start(11,6), stopp(11,6))
                                    do anb6 = start(11,6), stopp(11,6), steg(11&
                                       ,6)
                                    antel(11,6) = anb6 + antel(11,5)
                                    if (.not.(antel(11,6)<=antal .and. ansats(&
                                       11,5,1)<=2 .and. ansats(11,5,0)<=2)) &
                                       cycle
                                    do plusb6 = min(anb6,14), max(anb6 - 12,0)&
                                       , -1
                                    ansats(11,6,1) = plusb6
                                    ansats(11,6,0) = anb6 - plusb6
!     11k
                                    call slug (11, 7, varmax, varupp, varned, &
                                       ansats, org, lock(11,7), dubbel, low, &
                                       start(11,7), stopp(11,7))
                                    do anb7 = start(11,7), stopp(11,7), steg(11&
                                       ,7)
                                    antel(11,7) = anb7 + antel(11,6)
                                    if (.not.(antel(11,7)<=antal .and. ansats(&
                                       11,6,1)<=2 .and. ansats(11,6,0)<=2)) &
                                       cycle
                                    do plusb7 = min(anb7,16), max(anb7 - 14,0)&
                                       , -1
                                    ansats(11,7,1) = plusb7
                                    ansats(11,7,0) = anb7 - plusb7
!     11l
                                    call slug (11, 8, varmax, varupp, varned, &
                                       ansats, org, lock(11,8), dubbel, low, &
                                       start(11,8), stopp(11,8))
                                    do anb8 = start(11,8), stopp(11,8), steg(11&
                                       ,8)
                                    antel(11,8) = anb8 + antel(11,7)
                                    if (.not.(antel(11,8)<=antal .and. ansats(&
                                       11,7,1)<=2 .and. ansats(11,7,0)<=2)) &
                                       cycle
                                    do plusb8 = min(anb8,18), max(anb8 - 16,0)&
                                       , -1
                                    ansats(11,8,1) = plusb8
                                    ansats(11,8,0) = anb8 - plusb8
!     11m
                                    call slug (11, 9, varmax, varupp, varned, &
                                       ansats, org, lock(11,9), dubbel, low, &
                                       start(11,9), stopp(11,9))
                                    do anb9 = start(11,9), stopp(11,9), steg(11&
                                       ,9)
                                    antel(11,9) = anb9 + antel(11,8)
                                    if (.not.(antel(11,9)<=antal .and. ansats(&
                                       11,8,1)<=2 .and. ansats(11,8,0)<=2)) &
                                       cycle
                                    do plusb9 = min(anb9,20), max(anb9 - 18,0)&
                                       , -1
                                    ansats(11,9,1) = plusb9
                                    ansats(11,9,0) = anb9 - plusb9
!     11n
                                    call slug (11, 10, varmax, varupp, varned, &
                                       ansats, org, lock(11,10), dubbel, low, &
                                       start(11,10), stopp(11,10))
                                    do anba = start(11,10), stopp(11,10), steg(&
                                       11,10)
                                    antel(11,10) = anba + antel(11,9)
                                    if (.not.(antel(11,10)<=antal .and. ansats(&
                                       11,9,1)<=2 .and. ansats(11,9,0)<=2&
                                        .and. antel(11,10)>=lim(11))) cycle
                                    do plusba = min(anba,22), max(anba - 20,0)&
                                       , -1
                                    ansats(11,10,1) = plusba
                                    ansats(11,10,0) = anba - plusba
!     12s
                                    call slug (12, 0, varmax, varupp, varned, &
                                       ansats, org, lock(12,0), dubbel, low, &
                                       start(12,0), stopp(12,0))
                                    do anc0 = start(12,0), stopp(12,0), steg(12&
                                       ,0)
                                    antel(12,0) = anc0 + antel(11,10)
                                    if (.not.(antel(12,0)<=antal .and. ansats(&
                                       11,10,1)<=2 .and. ansats(11,10,0)<=2)) &
                                       cycle
                                    ansats(12,0,0) = anc0
!     12p
                                    call slug (12, 1, varmax, varupp, varned, &
                                       ansats, org, lock(12,1), dubbel, low, &
                                       start(12,1), stopp(12,1))
                                    do anc1 = start(12,1), stopp(12,1), steg(12&
                                       ,1)
                                    antel(12,1) = anc1 + antel(12,0)
                                    if (antel(12,1) > antal) cycle
                                    do plusc1 = min(anc1,4), max(anc1 - 2,0), &
                                       -1
                                    ansats(12,1,1) = plusc1
                                    ansats(12,1,0) = anc1 - plusc1
!     12d
                                    call slug (12, 2, varmax, varupp, varned, &
                                       ansats, org, lock(12,2), dubbel, low, &
                                       start(12,2), stopp(12,2))
                                    do anc2 = start(12,2), stopp(12,2), steg(12&
                                       ,2)
                                    antel(12,2) = anc2 + antel(12,1)
                                    if (antel(12,2) > antal) cycle
                                    do plusc2 = min(anc2,6), max(anc2 - 4,0), &
                                       -1
                                    ansats(12,2,1) = plusc2
                                    ansats(12,2,0) = anc2 - plusc2
!     12f
                                    call slug (12, 3, varmax, varupp, varned, &
                                       ansats, org, lock(12,3), dubbel, low, &
                                       start(12,3), stopp(12,3))
                                    do anc3 = start(12,3), stopp(12,3), steg(12&
                                       ,3)
                                    antel(12,3) = anc3 + antel(12,2)
                                    if (antel(12,3) > antal) cycle
                                    do plusc3 = min(anc3,8), max(anc3 - 6,0), &
                                       -1
                                    ansats(12,3,1) = plusc3
                                    ansats(12,3,0) = anc3 - plusc3
!     12g
                                    call slug (12, 4, varmax, varupp, varned, &
                                       ansats, org, lock(12,4), dubbel, low, &
                                       start(12,4), stopp(12,4))
                                    do anc4 = start(12,4), stopp(12,4), steg(12&
                                       ,4)
                                    antel(12,4) = anc4 + antel(12,3)
                                    if (antel(12,4) > antal) cycle
                                    do plusc4 = min(anc4,10), max(anc4 - 8,0), &
                                       -1
                                    ansats(12,4,1) = plusc4
                                    ansats(12,4,0) = anc4 - plusc4
!     12h
                                    call slug (12, 5, varmax, varupp, varned, &
                                       ansats, org, lock(12,5), dubbel, low, &
                                       start(12,5), stopp(12,5))
                                    do anc5 = start(12,5), stopp(12,5), steg(12&
                                       ,5)
                                    antel(12,5) = anc5 + antel(12,4)
                                    if (antel(12,5)>antal .or. ansats(12,4,1)>2&
                                       ) cycle
                                    do plusc5 = min(anc5,12), max(anc5 - 10,0)&
                                       , -1
                                    ansats(12,5,1) = plusc5
                                    ansats(12,5,0) = anc5 - plusc5
!     12i
                                    call slug (12, 6, varmax, varupp, varned, &
                                       ansats, org, lock(12,6), dubbel, low, &
                                       start(12,6), stopp(12,6))
                                    do anc6 = start(12,6), stopp(12,6), steg(12&
                                       ,6)
                                    antel(12,6) = anc6 + antel(12,5)
                                    if (.not.(antel(12,6)<=antal .and. ansats(&
                                       12,5,1)<=2 .and. ansats(12,5,0)<=2)) &
                                       cycle
                                    do plusc6 = min(anc6,14), max(anc6 - 12,0)&
                                       , -1
                                    ansats(12,6,1) = plusc6
                                    ansats(12,6,0) = anc6 - plusc6
!     12k
                                    call slug (12, 7, varmax, varupp, varned, &
                                       ansats, org, lock(12,7), dubbel, low, &
                                       start(12,7), stopp(12,7))
                                    do anc7 = start(12,7), stopp(12,7), steg(12&
                                       ,7)
                                    antel(12,7) = anc7 + antel(12,6)
                                    if (.not.(antel(12,7)<=antal .and. ansats(&
                                       12,6,1)<=2 .and. ansats(12,6,0)<=2)) &
                                       cycle
                                    do plusc7 = min(anc7,16), max(anc7 - 14,0)&
                                       , -1
                                    ansats(12,7,1) = plusc7
                                    ansats(12,7,0) = anc7 - plusc7
!     12l
                                    call slug (12, 8, varmax, varupp, varned, &
                                       ansats, org, lock(12,8), dubbel, low, &
                                       start(12,8), stopp(12,8))
                                    do anc8 = start(12,8), stopp(12,8), steg(12&
                                       ,8)
                                    antel(12,8) = anc8 + antel(12,7)
                                    if (.not.(antel(12,8)<=antal .and. ansats(&
                                       12,7,1)<=2 .and. ansats(12,7,0)<=2)) &
                                       cycle
                                    do plusc8 = min(anc8,18), max(anc8 - 16,0)&
                                       , -1
                                    ansats(12,8,1) = plusc8
                                    ansats(12,8,0) = anc8 - plusc8
!     12m
                                    call slug (12, 9, varmax, varupp, varned, &
                                       ansats, org, lock(12,9), dubbel, low, &
                                       start(12,9), stopp(12,9))
                                    do anc9 = start(12,9), stopp(12,9), steg(12&
                                       ,9)
                                    antel(12,9) = anc9 + antel(12,8)
                                    if (.not.(antel(12,9)<=antal .and. ansats(&
                                       12,8,1)<=2 .and. ansats(12,8,0)<=2)) &
                                       cycle
                                    do plusc9 = min(anc9,20), max(anc9 - 18,0)&
                                       , -1
                                    ansats(12,9,1) = plusc9
                                    ansats(12,9,0) = anc9 - plusc9
!     12n
                                    call slug (12, 10, varmax, varupp, varned, &
                                       ansats, org, lock(12,10), dubbel, low, &
                                       start(12,10), stopp(12,10))
                                    do anca = start(12,10), stopp(12,10), steg(&
                                       12,10)
                                    antel(12,10) = anca + antel(12,9)
                                    if (.not.(antel(12,10)<=antal .and. ansats(&
                                       12,9,1)<=2 .and. ansats(12,9,0)<=2&
                                        .and. antel(12,10)>=lim(12))) cycle
                                    do plusca = min(anca,22), max(anca - 20,0)&
                                       , -1
                                    ansats(12,10,1) = plusca
                                    ansats(12,10,0) = anca - plusca
!     13s
                                    call slug (13, 0, varmax, varupp, varned, &
                                       ansats, org, lock(13,0), dubbel, low, &
                                       start(13,0), stopp(13,0))
                                    do and0 = start(13,0), stopp(13,0), steg(13&
                                       ,0)
                                    antel(13,0) = and0 + antel(12,10)
                                    if (.not.(antel(13,0)<=antal .and. ansats(&
                                       12,10,1)<=2 .and. ansats(12,10,0)<=2)) &
                                       cycle
                                    ansats(13,0,0) = and0
!     13p
                                    call slug (13, 1, varmax, varupp, varned, &
                                       ansats, org, lock(13,1), dubbel, low, &
                                       start(13,1), stopp(13,1))
                                    do and1 = start(13,1), stopp(13,1), steg(13&
                                       ,1)
                                    antel(13,1) = and1 + antel(13,0)
                                    if (antel(13,1) > antal) cycle
                                    do plusd1 = min(and1,4), max(and1 - 2,0), &
                                       -1
                                    ansats(13,1,1) = plusd1
                                    ansats(13,1,0) = and1 - plusd1
!     13d
                                    call slug (13, 2, varmax, varupp, varned, &
                                       ansats, org, lock(13,2), dubbel, low, &
                                       start(13,2), stopp(13,2))
                                    do and2 = start(13,2), stopp(13,2), steg(13&
                                       ,2)
                                    antel(13,2) = and2 + antel(13,1)
                                    if (antel(13,2) > antal) cycle
                                    do plusd2 = min(and2,6), max(and2 - 4,0), &
                                       -1
                                    ansats(13,2,1) = plusd2
                                    ansats(13,2,0) = and2 - plusd2
!     13f
                                    call slug (13, 3, varmax, varupp, varned, &
                                       ansats, org, lock(13,3), dubbel, low, &
                                       start(13,3), stopp(13,3))
                                    do and3 = start(13,3), stopp(13,3), steg(13&
                                       ,3)
                                    antel(13,3) = and3 + antel(13,2)
                                    if (antel(13,3) > antal) cycle
                                    do plusd3 = min(and3,8), max(and3 - 6,0), &
                                       -1
                                    ansats(13,3,1) = plusd3
                                    ansats(13,3,0) = and3 - plusd3
!     13g
                                    call slug (13, 4, varmax, varupp, varned, &
                                       ansats, org, lock(13,4), dubbel, low, &
                                       start(13,4), stopp(13,4))
                                    do and4 = start(13,4), stopp(13,4), steg(13&
                                       ,4)
                                    antel(13,4) = and4 + antel(13,3)
                                    if (antel(13,4) > antal) cycle
                                    do plusd4 = min(and4,10), max(and4 - 8,0), &
                                       -1
                                    ansats(13,4,1) = plusd4
                                    ansats(13,4,0) = and4 - plusd4
!     13h
                                    call slug (13, 5, varmax, varupp, varned, &
                                       ansats, org, lock(13,5), dubbel, low, &
                                       start(13,5), stopp(13,5))
                                    do and5 = start(13,5), stopp(13,5), steg(13&
                                       ,5)
                                    antel(13,5) = and5 + antel(13,4)
                                    if (antel(13,5)>antal .or. ansats(13,4,1)>2&
                                       ) cycle
                                    do plusd5 = min(and5,12), max(and5 - 10,0)&
                                       , -1
                                    ansats(13,5,1) = plusd5
                                    ansats(13,5,0) = and5 - plusd5
!     13i
                                    call slug (13, 6, varmax, varupp, varned, &
                                       ansats, org, lock(13,6), dubbel, low, &
                                       start(13,6), stopp(13,6))
                                    do and6 = start(13,6), stopp(13,6), steg(13&
                                       ,6)
                                    antel(13,6) = and6 + antel(13,5)
                                    if (.not.(antel(13,6)<=antal .and. ansats(&
                                       13,5,1)<=2 .and. ansats(13,5,0)<=2)) &
                                       cycle
                                    do plusd6 = min(and6,14), max(and6 - 12,0)&
                                       , -1
                                    ansats(13,6,1) = plusd6
                                    ansats(13,6,0) = and6 - plusd6
!     13k
                                    call slug (13, 7, varmax, varupp, varned, &
                                       ansats, org, lock(13,7), dubbel, low, &
                                       start(13,7), stopp(13,7))
                                    do and7 = start(13,7), stopp(13,7), steg(13&
                                       ,7)
                                    antel(13,7) = and7 + antel(13,6)
                                    if (.not.(antel(13,7)<=antal .and. ansats(&
                                       13,6,1)<=2 .and. ansats(13,6,0)<=2)) &
                                       cycle
                                    do plusd7 = min(and7,16), max(and7 - 14,0)&
                                       , -1
                                    ansats(13,7,1) = plusd7
                                    ansats(13,7,0) = and7 - plusd7
!     13l
                                    call slug (13, 8, varmax, varupp, varned, &
                                       ansats, org, lock(13,8), dubbel, low, &
                                       start(13,8), stopp(13,8))
                                    do and8 = start(13,8), stopp(13,8), steg(13&
                                       ,8)
                                    antel(13,8) = and8 + antel(13,7)
                                    if (.not.(antel(13,8)<=antal .and. ansats(&
                                       13,7,1)<=2 .and. ansats(13,7,0)<=2)) &
                                       cycle
                                    do plusd8 = min(and8,18), max(and8 - 16,0)&
                                       , -1
                                    ansats(13,8,1) = plusd8
                                    ansats(13,8,0) = and8 - plusd8
!     13m
                                    call slug (13, 9, varmax, varupp, varned, &
                                       ansats, org, lock(13,9), dubbel, low, &
                                       start(13,9), stopp(13,9))
                                    do and9 = start(13,9), stopp(13,9), steg(13&
                                       ,9)
                                    antel(13,9) = and9 + antel(13,8)
                                    if (.not.(antel(13,9)<=antal .and. ansats(&
                                       13,8,1)<=2 .and. ansats(13,8,0)<=2)) &
                                       cycle
                                    do plusd9 = min(and9,20), max(and9 - 18,0)&
                                       , -1
                                    ansats(13,9,1) = plusd9
                                    ansats(13,9,0) = and9 - plusd9
!     13n
                                    call slug (13, 10, varmax, varupp, varned, &
                                       ansats, org, lock(13,10), dubbel, low, &
                                       start(13,10), stopp(13,10))
                                    do anda = start(13,10), stopp(13,10), steg(&
                                       13,10)
                                    antel(13,10) = anda + antel(13,9)
                                    if (.not.(antel(13,10)<=antal .and. ansats(&
                                       13,9,1)<=2 .and. ansats(13,9,0)<=2&
                                        .and. antel(13,10)>=lim(13))) cycle
                                    do plusda = min(anda,22), max(anda - 20,0)&
                                       , -1
                                    ansats(13,10,1) = plusda
                                    ansats(13,10,0) = anda - plusda
!     14s
                                    call slug (14, 0, varmax, varupp, varned, &
                                       ansats, org, lock(14,0), dubbel, low, &
                                       start(14,0), stopp(14,0))
                                    do ane0 = start(14,0), stopp(14,0), steg(14&
                                       ,0)
                                    antel(14,0) = ane0 + antel(13,10)
                                    if (.not.(antel(14,0)<=antal .and. ansats(&
                                       13,10,1)<=2 .and. ansats(13,10,0)<=2)) &
                                       cycle
                                    ansats(14,0,0) = ane0
!     14p
                                    call slug (14, 1, varmax, varupp, varned, &
                                       ansats, org, lock(14,1), dubbel, low, &
                                       start(14,1), stopp(14,1))
                                    do ane1 = start(14,1), stopp(14,1), steg(14&
                                       ,1)
                                    antel(14,1) = ane1 + antel(14,0)
                                    if (antel(14,1) > antal) cycle
                                    do pluse1 = min(ane1,4), max(ane1 - 2,0), &
                                       -1
                                    ansats(14,1,1) = pluse1
                                    ansats(14,1,0) = ane1 - pluse1
!     14d
                                    call slug (14, 2, varmax, varupp, varned, &
                                       ansats, org, lock(14,2), dubbel, low, &
                                       start(14,2), stopp(14,2))
                                    do ane2 = start(14,2), stopp(14,2), steg(14&
                                       ,2)
                                    antel(14,2) = ane2 + antel(14,1)
                                    if (antel(14,2) > antal) cycle
                                    do pluse2 = min(ane2,6), max(ane2 - 4,0), &
                                       -1
                                    ansats(14,2,1) = pluse2
                                    ansats(14,2,0) = ane2 - pluse2
!     14f
                                    call slug (14, 3, varmax, varupp, varned, &
                                       ansats, org, lock(14,3), dubbel, low, &
                                       start(14,3), stopp(14,3))
                                    do ane3 = start(14,3), stopp(14,3), steg(14&
                                       ,3)
                                    antel(14,3) = ane3 + antel(14,2)
                                    if (antel(14,3) > antal) cycle
                                    do pluse3 = min(ane3,8), max(ane3 - 6,0), &
                                       -1
                                    ansats(14,3,1) = pluse3
                                    ansats(14,3,0) = ane3 - pluse3
!     14g
                                    call slug (14, 4, varmax, varupp, varned, &
                                       ansats, org, lock(14,4), dubbel, low, &
                                       start(14,4), stopp(14,4))
                                    do ane4 = start(14,4), stopp(14,4), steg(14&
                                       ,4)
                                    antel(14,4) = ane4 + antel(14,3)
                                    if (antel(14,4) > antal) cycle
                                    do pluse4 = min(ane4,10), max(ane4 - 8,0), &
                                       -1
                                    ansats(14,4,1) = pluse4
                                    ansats(14,4,0) = ane4 - pluse4
!     14h
                                    call slug (14, 5, varmax, varupp, varned, &
                                       ansats, org, lock(14,5), dubbel, low, &
                                       start(14,5), stopp(14,5))
                                    do ane5 = start(14,5), stopp(14,5), steg(14&
                                       ,5)
                                    antel(14,5) = ane5 + antel(14,4)
                                    if (antel(14,5)>antal .or. ansats(14,4,1)>2&
                                       ) cycle
                                    do pluse5 = min(ane5,12), max(ane5 - 10,0)&
                                       , -1
                                    ansats(14,5,1) = pluse5
                                    ansats(14,5,0) = ane5 - pluse5
!     14i
                                    call slug (14, 6, varmax, varupp, varned, &
                                       ansats, org, lock(14,6), dubbel, low, &
                                       start(14,6), stopp(14,6))
                                    do ane6 = start(14,6), stopp(14,6), steg(14&
                                       ,6)
                                    antel(14,6) = ane6 + antel(14,5)
                                    if (.not.(antel(14,6)<=antal .and. ansats(&
                                       14,5,1)<=2 .and. ansats(14,5,0)<=2)) &
                                       cycle
                                    do pluse6 = min(ane6,14), max(ane6 - 12,0)&
                                       , -1
                                    ansats(14,6,1) = pluse6
                                    ansats(14,6,0) = ane6 - pluse6
!     14k
                                    call slug (14, 7, varmax, varupp, varned, &
                                       ansats, org, lock(14,7), dubbel, low, &
                                       start(14,7), stopp(14,7))
                                    do ane7 = start(14,7), stopp(14,7), steg(14&
                                       ,7)
                                    antel(14,7) = ane7 + antel(14,6)
                                    if (.not.(antel(14,7)<=antal .and. ansats(&
                                       14,6,1)<=2 .and. ansats(14,6,0)<=2)) &
                                       cycle
                                    do pluse7 = min(ane7,16), max(ane7 - 14,0)&
                                       , -1
                                    ansats(14,7,1) = pluse7
                                    ansats(14,7,0) = ane7 - pluse7
!     14l
                                    call slug (14, 8, varmax, varupp, varned, &
                                       ansats, org, lock(14,8), dubbel, low, &
                                       start(14,8), stopp(14,8))
                                    do ane8 = start(14,8), stopp(14,8), steg(14&
                                       ,8)
                                    antel(14,8) = ane8 + antel(14,7)
                                    if (.not.(antel(14,8)<=antal .and. ansats(&
                                       14,7,1)<=2 .and. ansats(14,7,0)<=2)) &
                                       cycle
                                    do pluse8 = min(ane8,18), max(ane8 - 16,0)&
                                       , -1
                                    ansats(14,8,1) = pluse8
                                    ansats(14,8,0) = ane8 - pluse8
!     14m
                                    call slug (14, 9, varmax, varupp, varned, &
                                       ansats, org, lock(14,9), dubbel, low, &
                                       start(14,9), stopp(14,9))
                                    do ane9 = start(14,9), stopp(14,9), steg(14&
                                       ,9)
                                    antel(14,9) = ane9 + antel(14,8)
                                    if (.not.(antel(14,9)<=antal .and. ansats(&
                                       14,8,1)<=2 .and. ansats(14,8,0)<=2)) &
                                       cycle
                                    do pluse9 = min(ane9,20), max(ane9 - 18,0)&
                                       , -1
                                    ansats(14,9,1) = pluse9
                                    ansats(14,9,0) = ane9 - pluse9
!     14n
                                    call slug (14, 10, varmax, varupp, varned, &
                                       ansats, org, lock(14,10), dubbel, low, &
                                       start(14,10), stopp(14,10))
                                    do anea = start(14,10), stopp(14,10), steg(&
                                       14,10)
                                    antel(14,10) = anea + antel(14,9)
                                    if (.not.(antel(14,10)<=antal .and. ansats(&
                                       14,9,1)<=2 .and. ansats(14,9,0)<=2&
                                        .and. antel(14,10)>=lim(14))) cycle
                                    do plusea = min(anea,22), max(anea - 20,0)&
                                       , -1
                                    ansats(14,10,1) = plusea
                                    ansats(14,10,0) = anea - plusea
!     15s
                                    call slug (15, 0, varmax, varupp, varned, &
                                       ansats, org, lock(15,0), dubbel, low, &
                                       start(15,0), stopp(15,0))
                                    do anf0 = start(15,0), stopp(15,0), steg(15&
                                       ,0)
                                    antel(15,0) = anf0 + antel(14,10)
                                    if (.not.(antel(15,0)<=antal .and. ansats(&
                                       14,10,1)<=2 .and. ansats(14,10,0)<=2)) &
                                       cycle
                                    ansats(15,0,0) = anf0
!     15p
                                    call slug (15, 1, varmax, varupp, varned, &
                                       ansats, org, lock(15,1), dubbel, low, &
                                       start(15,1), stopp(15,1))
                                    do anf1 = start(15,1), stopp(15,1), steg(15&
                                       ,1)
                                    antel(15,1) = anf1 + antel(15,0)
                                    if (antel(15,1) > antal) cycle
                                    do plusf1 = min(anf1,4), max(anf1 - 2,0), &
                                       -1
                                    ansats(15,1,1) = plusf1
                                    ansats(15,1,0) = anf1 - plusf1
!     15d
                                    call slug (15, 2, varmax, varupp, varned, &
                                       ansats, org, lock(15,2), dubbel, low, &
                                       start(15,2), stopp(15,2))
                                    do anf2 = start(15,2), stopp(15,2), steg(15&
                                       ,2)
                                    antel(15,2) = anf2 + antel(15,1)
                                    if (antel(15,2) > antal) cycle
                                    do plusf2 = min(anf2,6), max(anf2 - 4,0), &
                                       -1
                                    ansats(15,2,1) = plusf2
                                    ansats(15,2,0) = anf2 - plusf2
!     15f
                                    call slug (15, 3, varmax, varupp, varned, &
                                       ansats, org, lock(15,3), dubbel, low, &
                                       start(15,3), stopp(15,3))
                                    do anf3 = start(15,3), stopp(15,3), steg(15&
                                       ,3)
                                    antel(15,3) = anf3 + antel(15,2)
                                    if (antel(15,3) > antal) cycle
                                    do plusf3 = min(anf3,8), max(anf3 - 6,0), &
                                       -1
                                    ansats(15,3,1) = plusf3
                                    ansats(15,3,0) = anf3 - plusf3
!     15g
                                    call slug (15, 4, varmax, varupp, varned, &
                                       ansats, org, lock(15,4), dubbel, low, &
                                       start(15,4), stopp(15,4))
                                    do anf4 = start(15,4), stopp(15,4), steg(15&
                                       ,4)
                                    antel(15,4) = anf4 + antel(15,3)
                                    if (antel(15,4) > antal) cycle
                                    do plusf4 = min(anf4,10), max(anf4 - 8,0), &
                                       -1
                                    ansats(15,4,1) = plusf4
                                    ansats(15,4,0) = anf4 - plusf4
!     15h
                                    call slug (15, 5, varmax, varupp, varned, &
                                       ansats, org, lock(15,5), dubbel, low, &
                                       start(15,5), stopp(15,5))
                                    do anf5 = start(15,5), stopp(15,5), steg(15&
                                       ,5)
                                    antel(15,5) = anf5 + antel(15,4)
                                    if (antel(15,5)>antal .or. ansats(15,4,1)>2&
                                       ) cycle
                                    do plusf5 = min(anf5,12), max(anf5 - 10,0)&
                                       , -1
                                    ansats(15,5,1) = plusf5
                                    ansats(15,5,0) = anf5 - plusf5
!     15i
                                    call slug (15, 6, varmax, varupp, varned, &
                                       ansats, org, lock(15,6), dubbel, low, &
                                       start(15,6), stopp(15,6))
                                    do anf6 = start(15,6), stopp(15,6), steg(15&
                                       ,6)
                                    antel(15,6) = anf6 + antel(15,5)
                                    if (.not.(antel(15,6)<=antal .and. ansats(&
                                       15,5,1)<=2 .and. ansats(15,5,0)<=2)) &
                                       cycle
                                    do plusf6 = min(anf6,14), max(anf6 - 12,0)&
                                       , -1
                                    ansats(15,6,1) = plusf6
                                    ansats(15,6,0) = anf6 - plusf6
!     15k
                                    call slug (15, 7, varmax, varupp, varned, &
                                       ansats, org, lock(15,7), dubbel, low, &
                                       start(15,7), stopp(15,7))
                                    do anf7 = start(15,7), stopp(15,7), steg(15&
                                       ,7)
                                    antel(15,7) = anf7 + antel(15,6)
                                    if (.not.(antel(15,7)<=antal .and. ansats(&
                                       15,6,1)<=2 .and. ansats(15,6,0)<=2)) &
                                       cycle
                                    do plusf7 = min(anf7,16), max(anf7 - 14,0)&
                                       , -1
                                    ansats(15,7,1) = plusf7
                                    ansats(15,7,0) = anf7 - plusf7
!     15l
                                    call slug (15, 8, varmax, varupp, varned, &
                                       ansats, org, lock(15,8), dubbel, low, &
                                       start(15,8), stopp(15,8))
                                    do anf8 = start(15,8), stopp(15,8), steg(15&
                                       ,8)
                                    antel(15,8) = anf8 + antel(15,7)
                                    if (.not.(antel(15,8)<=antal .and. ansats(&
                                       15,7,1)<=2 .and. ansats(15,7,0)<=2)) &
                                       cycle
                                    do plusf8 = min(anf8,18), max(anf8 - 16,0)&
                                       , -1
                                    ansats(15,8,1) = plusf8
                                    ansats(15,8,0) = anf8 - plusf8
!     15m
                                    call slug (15, 9, varmax, varupp, varned, &
                                       ansats, org, lock(15,9), dubbel, low, &
                                       start(15,9), stopp(15,9))
                                    do anf9 = start(15,9), stopp(15,9), steg(15&
                                       ,9)
                                    antel(15,9) = anf9 + antel(15,8)
                                    if (.not.(antel(15,9)<=antal .and. ansats(&
                                       15,8,1)<=2 .and. ansats(15,8,0)<=2)) &
                                       cycle
                                    do plusf9 = min(anf9,20), max(anf9 - 18,0)&
                                       , -1
                                    ansats(15,9,1) = plusf9
                                    ansats(15,9,0) = anf9 - plusf9
!     15n
                                    call slug (15, 10, varmax, varupp, varned, &
                                       ansats, org, lock(15,10), dubbel, low, &
                                       start(15,10), stopp(15,10))
                                    do anfa = start(15,10), stopp(15,10), steg(&
                                       15,10)
                                    antel(15,10) = anfa + antel(15,9)
                                    if (.not.(antel(15,10)==antal .and. ansats(&
                                       15,9,1)<=2 .and. ansats(15,9,0)<=2)) &
                                       cycle
                                    do plusfa = min(anfa,22), max(anfa - 20,0)&
                                       , -1
                                    ansats(15,10,1) = plusfa
                                    ansats(15,10,0) = anfa - plusfa
                                    if (ansats(15,10,1)>2 .or. ansats(15,10,0)>&
                                       2) cycle
                                    par = 0
                                    elar = 0
                                    do i = 1, 15
                                    do j = 0, min(10,i - 1)
                                    do k = 0, min(j,1)
                                    elar = elar + ansats(i,j,k)
                                    par = mod(par + j*ansats(i,j,k),2)
                                    end do
                                    end do
                                    end do
                                    if (par /= par0) cycle
                                    if (elar == antal) then
                                    call gen (ansats, posn, posl, skal, cf, &
                                       first, minj, maxj, par0)
                                    else
                                    write (*, *) 'FEL'
                                    endif
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      if (first) then
         rewind (fil_1)
      else
         rewind (fil_2)
      endif
      if (cf == 0) then
         write (*, 1005) 'No configuration state has been generated.'
      else if (cf == 1) then
         write (*, 1005) 'One configuration state has been generated.'
      else if (cf < 10) then
         write (*, 1001) cf, ' configuration states have been generated.'
      else if (cf < 100) then
         write (*, 1002) cf, ' configuration states have been generated.'
      else if (cf < 1000) then
         write (*, 1003) cf, ' configuration states have been generated.'
      else if (cf < 10000) then
         write (*, 1004) cf, ' configuration states have been generated.'
      else if (cf < 100000) then
         write (*, 1006) cf, ' configuration states have been generated.'
      else
         write (*, *) cf, ' configuration states have been generated.'
      endif
! 1000 format(A)
 1001 format(' ',i1,a)
 1002 format(' ',i2,a)
 1003 format(' ',i3,a)
 1004 format(' ',i4,a)
 1005 format(' ',a)
 1006 format(' ',i5,a)
 5000 format(11i2)
      return
      end subroutine blanda
