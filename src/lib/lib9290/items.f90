!***********************************************************************
!
      subroutine items(ncmin, ncf, record, ierr)
!
! Purpose:
!   Parse a list of levels stored in char string record. Formats
!   accepted are: 1, 2, 5, 7-10
!
! Input:
!   ncf, record
!
! Output:
!   ierr = 0 - normal
!          1 - no state
!         -1 - cannot decode
!         -2 - out of range (1,ncf)
!         -3 - cannot decode format like 8-3
!
! Input/Output:
!   ncmin - on input, length of array iccmin() (0 then not allocated);
!           on outpu, new length of array iccmin().
!
! Note:
!   This complicated ncmin is due to consideration for probable
!   multiple calls of this routine (block version).
!   The caller is responsible for deallocation of pccmin
!
!   Extracted and modified from rscf92/getold.f
!   by  Xinghong He                                        Jul 17 1997
!   Parameter ncd removed for simpler structure            Jun 10 1998
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  19:51:17   2/16/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE memory_man
      USE def_C, ONLY: ICCMIN
      USE iounit_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt_I
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: ncmin
      integer, intent(in) :: ncf
      integer, intent(out) :: ierr
      CHARACTER (LEN = *), INTENT(IN) :: RECORD
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: length_of_record, ncmin_in, ncd, ifirst, istart, i, iend, &
         isize, lenth, ios, level, level1, number, j
      character :: reci, form*7, cnum*3
!-----------------------------------------------
!
!

      length_of_record = len_trim(record)

      ncmin_in = ncmin
      ncd = ncmin                           ! Current length of array iccmin()

      if (ncd == 0) then
         ncd = 1
         call alloc (iccmin, ncd, 'ICCMIN', 'ITEMS')
      endif
!
!   parse record from left to right
!
      ifirst = 0
      istart = 1
      i = 1
!
!   .. skip the blanks and commas(this implementation allows input to
!      start with blanks
!
    2 continue
      reci = record(i:i)
      if (reci/=' ' .and. reci/=',') then
         istart = i
      else
         i = i + 1
         if (i <= length_of_record) then
            go to 2
         else
            go to 4
         endif
      endif
!
!   .. search for end of string (blank, comma, or dash)
!
    3 continue
      reci = record(i:i)
      if (reci/=' ' .and. reci/=',' .and. reci/='-') then
         i = i + 1
         if (i <= length_of_record) go to 3
      endif
!
!     ... read integer
!
      iend = i - 1
      isize = iend - istart + 1
      call convrt (isize, cnum, lenth)
      form = '(1i'//cnum(1:lenth)//')'
      read (record(istart:iend), form, iostat=ios) level
      if (ios /= 0) then
         write (istde, *) 'items: unable to decode '//record(istart:iend)//';'
         ierr = -1
         return
      endif

      if (ifirst == 0) then
!       .. this is the either the first or an isolated level
         ncmin = ncmin + 1
         if (ncmin > ncd) then
            call ralloc (iccmin, ncmin, 'ICCMIN', 'ITEMS')
            ncd = ncmin
         endif

         if (level<1 .or. level>ncf) then
            write (istde, *) 'items: serial numbers must be', &
               ' in the range [1,', ncf, '];'
            ierr = -2
            return
         endif

         iccmin(ncmin) = level
         i = i + 1
         if (reci == '-') ifirst = ncmin
         go to 2
      else
!        .. the previous level was the beginning of a range
         level1 = iccmin(ncmin)
         number = level - level1

         if (number < 0) then
            write (istde, *) level1, '-', level, ' not allowed'
            ierr = -3
            return
         endif

         ncmin = ncmin + number
         if (ncmin > ncd) then
            call ralloc (iccmin, ncmin, 'ICCMIN', 'ITEMS')
            ncd = ncmin
         endif
         do j = 1, number
            iccmin(ifirst+j) = level1 + j
         end do
         i = i + 1
         ifirst = 0
         go to 2
      endif

!   at least one level must be requested

    4 continue
      if (ncmin == ncmin_in) then
         ierr = 1
         return
      endif
!
!   trim array to exactly the correct size
!
      if (ncmin /= ncd) then
         write (istde, *) 'items: ncmin .ne. ncd'
         call ralloc (iccmin, ncmin, 'ICCMIN', 'ITEMS')
      endif

      ierr = 0
      return
      end subroutine items
