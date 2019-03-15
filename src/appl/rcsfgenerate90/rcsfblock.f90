!***********************************************************************
      SUBROUTINE RCSFBLOCK
!***********************************************************************
! Purpose:
!  Read in a CSF file 'rcsf.inp' and then:
!   1. Group into blocks
!   2. Detect and remove possible dubplicates within a block (Optional)
!   3. Generate a new CSF file 'rcsf.out'
! Limit:
!  Maximum number of blocks: 100
!  Maximum number of CSF's for a block: 1,000,000
!  Maximum line length of the CSL list: 100
!
! Written by Charlotte F Fischer
! Updated by Xinghong He on 98-10-29
!***********************************************************************
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      PARAMETER (MAX_CSF = 1000000)
      CHARACTER*300 LINE1, LINE2, LINE3
      CHARACTER*1500 string(5),string4
      CHARACTER(LEN=*), PARAMETER:: input='rcsf.out', output='rcsf.out2'
      CHARACTER*5  blkid(100), jblk
      INTEGER indx(100),norb(300)
      INTEGER ncfblk(100)
      INTEGER nblock
      DATA blkid /100*'     '/
      DATA ncfblk /100*0/, indx /100*0/

      !For checking duplicates
      LOGICAL CheckDuplicate, getyn
      CHARACTER*180 a1(MAX_CSF), b1(MAX_CSF), c1(MAX_CSF)
      CHARACTER*180 aa1(MAX_CSF), bb1(MAX_CSF), cc1(MAX_CSF),ct
!***********************************************************************

      write(*,*)
      write(*,*) 'Group CSFs into symmetry blocks'
      write(*,*)

!      WRITE (0,*) 'Perform duplicate check and remove them ?'
!      CheckDuplicate = getyn()
      CheckDuplicate = .false.

      OPEN (19, FILE = input, STATUS = 'OLD', FORM = 'FORMATTED',   &
                ACTION='READWRITE')
!      OPEN (20, FILE = output, STATUS = 'UNKNOWN', FORM = 'FORMATTED')

      nblock = 0

!***************************************************
! Headers, read and write without modification
!***************************************************

      DO i = 1, 5
         READ (19,'(A)') string(i)
      ENDDO

! norb(i) = 0 indicates that an on line 4 is not
! present in the CSFs

      do i = 1,300
         norb(i) = 0
      end do

      NLENGTH = LEN_TRIM(string(4)) + 1

!***************************************************
! Read in all CSF's and place them in different scratch files
! for different blocks
!***************************************************

   10 CONTINUE
      READ (19, '(A)', END = 99) LINE1
      IF (LINE1(1:2) .EQ. ' *') READ (19, '(A)') LINE1

! Check what orbitals are present in current CSF

      do i = 1,nlength/5
         n = index(line1,string(4)(5*(i-1)+1:5*i))
         if (n.ne.0) then
            norb(i) = 1
         end if
      end do

      READ (19, '(A)') LINE2
      READ (19, '(A)') LINE3

      i = LEN_TRIM(line3)
      jblk = line3(i-4:i)

      iunit = 0
      DO i = 1, nblock
         IF (jblk .EQ. blkid(i)) THEN
            iunit = 6000 + i
            ithisblk = i
         ENDIF
      ENDDO

      IF (iunit .EQ. 0) THEN
         !... we have a new block
         nblock = nblock + 1
         ithisblk = nblock
         blkid(nblock) = jblk
         iunit = 6000 + nblock
         OPEN (iunit, STATUS = 'SCRATCH', FORM = 'FORMATTED')
      ENDIF

      WRITE (iunit,'(A)') TRIM(line1)
      WRITE (iunit,'(A)') TRIM(line2)
      WRITE (iunit,'(A)') TRIM(line3)

      ncfblk(ithisblk) = ncfblk(ithisblk) + 1
      GOTO 10

   99 CONTINUE

      write(*,'(i2,a)')  nblock, ' blocks were created'
      write(*,*)
      write(*,*) '      block  J/P            NCSF'

!************************************************************
! Sort the the block order so that it always comes out like
!  0, 1/2, 1, 3/2, 2, 5/2 ...
!************************************************************

      CALL INDEXA (nblock, blkid, .TRUE., indx)


!************************************************************
!  No detection for the possible duplicates, directly write to
!  the rcsl.out file
!************************************************************

      REWIND(19)

      DO ii = 1, 3
          WRITE (19,'(A)') TRIM(string(ii))
      ENDDO

! Write the orbitals that have been found in the CSFs

      do ii = 1,1500
         string4(ii:ii) = ' '
      end do

      iii = 0
      do ii = 1,nlength/5
         if (norb(ii).eq.1) then
            iii = iii + 1
            string4(5*(iii-1)+1:5*iii) = string(4)(5*(ii-1)+1:5*ii)
         end if
      end do
      write (19,'(A)') trim(string4)

      WRITE (19,'(A)') TRIM(string(5))


      DO ii = nblock, 1, -1
         idx = indx(ii)
         PRINT *, nblock-ii+1, blkid(idx),'   ', ncfblk(idx)
         iunit = 6000 + idx
         REWIND (iunit)
  200    READ (iunit, '(A)', END = 9999) line1
         WRITE (19,'(A)') TRIM(line1)
         GOTO 200
 9999    CLOSE (iunit)

         IF (ii .NE. 1) THEN
            WRITE (19,'(A)') ' *'
         END IF
      ENDDO

      END
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE indexa(n,a,ldown,indx)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IMPLICIT NONE
!Revision 1.2  2004/04/23 20:23:01  cff
!
!Revision 1.1.1.1  2003/01/04 21:45:39  georgio
!
!     CHARACTER*(*)    RCSID
!     PARAMETER        ( RCSID
!    & ='$Id: jsplit.f,v 1.2 2004/04/23 20:23:01 cff Exp $'
!    & )
      LOGICAL          ldown              ! .TRUE. then Big ---> Small
      INTEGER          n, indx(n)

      !DOUBLE PRECISION a(n), aimx
      CHARACTER*(*) a(n)
      CHARACTER (LEN=LEN_TRIM(a(1))) aimx

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Sort out the order of array a and store the index in indx (a pointer)
! The input array a is unchanged
! written in the bases of UpDown
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INTEGER          i, j, ipos, jpos, jhere

      ! Initialize the index array
      DO i = 1, n
         indx(i) = i
      ENDDO

      IF (ldown) THEN
         DO i = 1, n
            ipos = indx(i)
            aimx = a(ipos)
            jhere = i
            DO j = i+1, n
               jpos = indx(j)
               IF(a(jpos) .GT. aimx) THEN
                  aimx = a(jpos)
                  jhere = j
               ENDIF
            ENDDO
            indx(i) = indx(jhere)
            indx(jhere) = ipos
         ENDDO
      ELSE
         DO i = 1, n
            ipos = indx(i)
            aimx = a(ipos)
            jhere = i
            DO j = i+1, n
               jpos = indx(j)
               IF(a(jpos) .LT. aimx) THEN
	               aimx = a(jpos)
                  jhere = j
               ENDIF
            ENDDO
            indx(i) = indx(jhere)
            indx(jhere) = ipos
         ENDDO
      ENDIF
      RETURN
      END

