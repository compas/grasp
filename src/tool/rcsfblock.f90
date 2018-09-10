      PROGRAM RCSFBLOCK
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
      PARAMETER (MAX_CSF = 1000000)
      CHARACTER*300 LINE1, LINE2, LINE3
      CHARACTER*500 string
      CHARACTER(LEN=*), PARAMETER:: input='rcsf.inp', output='rcsf.out'
      CHARACTER*5  blkid(100), jblk
      INTEGER indx(100)
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
      write(*,*) 'RCSFBLOCK'
      write(*,*) 'This program groups CSFs into symmetry blocks'
      write(*,*) 'Inputfile:  rcsf.inp'
      write(*,*) 'Outputfile: rcsf.out'
      write(*,*)

      WRITE (0,*) 'Perform duplicate check and remove them ?'
      CheckDuplicate = getyn()

      OPEN (19, FILE = input, STATUS = 'OLD', FORM = 'FORMATTED')
      OPEN (20, FILE = output, STATUS = 'UNKNOWN', FORM = 'FORMATTED')

      nblock = 0

!***************************************************
! Headers, read and write without modification
!***************************************************

      DO i = 1, 5
         READ (19,'(A)') string
         WRITE (20,'(A)') TRIM(string)
      ENDDO

!***************************************************
! Read in all CSF's and place them in different scratch files
! for different blocks
!***************************************************

   10 CONTINUE
      READ (19, '(A)', END = 99) LINE1
      IF (LINE1(1:2) .EQ. ' *') READ (19, '(A)') LINE1
      READ (19, '(A)') LINE2
      READ (19, '(A)') LINE3

      i = LEN_TRIM(line3)
      jblk = line3(i-4:i)

      iunit = 0
      DO i = 1, nblock
         IF (jblk .EQ. blkid(i)) THEN
            iunit = 20 + i
            ithisblk = i
         ENDIF
      ENDDO

      IF (iunit .EQ. 0) THEN
         !... we have a new block
         nblock = nblock + 1
         ithisblk = nblock
         blkid(nblock) = jblk
         iunit = 20 + nblock
         OPEN (iunit, STATUS = 'SCRATCH', FORM = 'FORMATTED')
      ENDIF

      WRITE (iunit,'(A)') TRIM(line1)
      WRITE (iunit,'(A)') TRIM(line2)
      WRITE (iunit,'(A)') TRIM(line3)

      ncfblk(ithisblk) = ncfblk(ithisblk) + 1
      GOTO 10

   99 CONTINUE

      write(*,*)  nblock, ' blocks were found'
      write(*,*) '      block  J/P            NCSF'

!************************************************************
! Sort the the block order so that it always comes out like
!  0, 1/2, 1, 3/2, 2, 5/2 ...
!************************************************************

      CALL INDEX (nblock, blkid, .TRUE., indx)

!************************************************************
!  Detect and remove possible duplicates, and then write to
!  the rcsl.out file
!************************************************************

      IF (CheckDuplicate) THEN
         DO ii = nblock, 1, -1
            idx = indx(ii)
            PRINT *
            PRINT *, nblock-ii+1, blkid(idx),'   ', ncfblk(idx)
            iunit = 20 + idx
            REWIND (iunit)
            !  Read in all records of the current block
            i = 1
  100       READ (iunit, '(A)', END = 999) a1(i)
            READ (iunit, '(A)') b1(i)
            READ (iunit, '(A)') c1(i)
            i = i + 1
            GOTO 100
  999       CLOSE (iunit)

            ncfg = i - 1

            do k=1,ncfg
              aa1(k)=' '
              bb1(k)=' '
              cc1(k)=' '
              i1=0
               ct=a1(k)
              do i=1,299
               if(ct(i:i).ne.' ') then
                 i1=i1+1
                 ct(i1:i1)=ct(i:i)
               else
                 if(ct(i+1:i+1).ne.' ') then
                   i1=i1+1
                   ct(i1:i1)=ct(i:i)
                 endif
               endif
              enddo
              aa1(k)=ct(1:i1)

              i1=0
               ct=b1(k)
              do i=1,299
               if(ct(i:i).ne.' ') then
                 i1=i1+1
                 ct(i1:i1)=ct(i:i)
               else
                 if(ct(i+1:i+1).ne.' ') then
                   i1=i1+1
                   ct(i1:i1)=ct(i:i)
                 endif
               endif
              enddo
              bb1(k)=ct(1:i1)

              i1=0
               ct=c1(k)
              do i=1,299
               if(ct(i:i).ne.' ') then
                 i1=i1+1
                 ct(i1:i1)=ct(i:i)
               else
                 if(ct(i+1:i+1).ne.' ') then
                   i1=i1+1
                   ct(i1:i1)=ct(i:i)
                 endif
               endif
              enddo
              cc1(k)=ct(1:i1)
            enddo

            !  Check for possible duplicate and remove it
            nn = 0
            DO k = 1, ncfg
               n = 0
               DO i = k+1, ncfg
                  IF (aa1(k) .EQ. aa1(i) .AND. bb1(k) .EQ. bb1(i) .AND. &
                     cc1(k) .EQ. cc1(i)) THEN
                     n = 1
                  ENDIF
               ENDDO
               ! Write only the unique records
               IF (n .EQ. 0) THEN
                  WRITE (20,'(A)') TRIM(a1(k))
                  WRITE (20,'(A)') TRIM(b1(k))
                  WRITE (20,'(A)') TRIM(c1(k))
               else
                  WRITE (90,'(i5)') k
                  WRITE (90,'(A)') TRIM(aa1(k))
                  WRITE (90,'(A)') TRIM(bb1(k))
                  WRITE (90,'(A)') TRIM(cc1(k))
                     nn = nn + 1
               ENDIF
            ENDDO
            PRINT *, '--- ', nn,' duplicates removed from the list'

            IF (ii .NE. 1) WRITE (20,'(A)') ' *'
         ENDDO

      ELSE

!************************************************************
!  No detection for the possible duplicates, directly write to
!  the rcsl.out file
!************************************************************

         DO ii = nblock, 1, -1
            idx = indx(ii)
            PRINT *, nblock-ii+1, blkid(idx),'   ', ncfblk(idx)
            iunit = 20 + idx
            REWIND (iunit)
  200       READ (iunit, '(A)', END = 9999) line1
            WRITE (20,'(A)') TRIM(line1)
            GOTO 200
 9999       CLOSE (iunit)

            IF (ii .NE. 1) WRITE (20,'(A)') ' *'
         ENDDO

      ENDIF
     
      CONTAINS

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                     index
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE index(n,a,ldown,indx)
      IMPLICIT NONE
!
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
      END SUBROUTINE
    
      END  PROGRAM
