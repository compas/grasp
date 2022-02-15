program lib9290_quad
    implicit none

    call setup
    call runtest

contains

    subroutine runtest
        use parameter_def
        use grid_C
        use tatb_C
        use quad_I

        real*8 :: real64_kind_
        integer, parameter :: real64 = kind(real64_kind_)
        integer, parameter :: dp = real64

        integer :: i
        real(real64) :: result

        ! Integrate exp(-2x) from 0 to infinity. Expected result is 1/2.
        do i = 1, NNN1
            TA(i) = exp(-2*R(i)) * RP(i)
        end do
        call quad(result)

        print *, "RESULT from QUAD", result
        print *, "expected        ", 0.5_dp
        print *, "difference      ", abs(result - 0.5_dp)

        if(abs(result - 0.5_dp) > 1e-12_dp) then
            stop 1
        end if
    end subroutine runtest

    subroutine setup
        use parameter_def
        use grid_C
        use tatb_C
        use setqic_I
        use radgrd_I

        H = 5.0D-2
        RNT = 2.0D-6
        HP = 0.0D0
        N = NNNP
        MTP = NNNP

        call setqic
        call radgrd
    end subroutine setup

end program lib9290_quad
