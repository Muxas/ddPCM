!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/dd_operators.f90
!! Tests for dd_operators module
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-11

program test_dd_operators
use dd_operators
implicit none

integer :: i, iprint = 1
! Some implementation do not work for alpha=1d-307, so range of test values is
! reduced. This can be fixed by enforcing proper order of multiplications like
! a*b*c into a*(b*c).
real(dp) :: alpha(4)=(/1d0, -1d0, 1d-298, 1d+300/)

contains

end program

