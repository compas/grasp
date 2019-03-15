      MODULE jlabl_C
!...Created by Pacific-Sierra Research 77to90  4.3E  06:16:25   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      CHARACTER(LEN=4), DIMENSION(32) :: JLBL, JLBR
      CHARACTER(LEN=4), DIMENSION(2) :: JLBP
!
!   Left-justified strings
!
      DATA JLBL/'0   ', '1/2 ', '1   ', '3/2 ', '2   ', '5/2 ', '3   ', '7/2 ',&
                '4   ', '9/2 ', '5   ', '11/2', '6   ', '13/2', '7   ', '15/2',&
                '8   ', '17/2', '9   ', '19/2', '10  ', '21/2', '11  ', '23/2',&
                '12  ', '25/2', '13  ', '27/2', '14  ', '29/2', '15  ', '31/2'/
!
!   Right-justified strings
!
      DATA JLBR/'   0', ' 1/2', '   1', ' 3/2', '   2', ' 5/2', '   3', ' 7/2',&
                '   4', ' 9/2', '   5', '11/2', '   6', '13/2', '   7', '15/2',&
                '   8', '17/2', '   9', '19/2', '  10', '21/2', '  11', '23/2',&
                '  12', '25/2', '  13', '27/2', '  14', '29/2', '  15', '31/2'/
!
!   Parity signs
!
      DATA JLBP/ ' -  ', ' +  '/
!

      END MODULE jlabl_C
