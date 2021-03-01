!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/main.f90
!! Main ddx driver.
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-25

program main
use dd_core
use dd_operators
use dd_solvers
use dd_cosmo
use dd_pcm
implicit none

character(len=255) :: fname
type(dd_data_type) :: dd_data
integer :: info
real(dp), allocatable :: phi_cav(:), gradphi_cav(:, :), psi(:, :), force(:, :)
real(dp) :: esolv, start_time, finish_time
integer :: i, j, isph

! Read input file name
call getarg(1, fname)
write(*, *) "Using provided file ", trim(fname), " as a config file"
call ddfromfile(fname, dd_data, info)
if(info .ne. 0) stop "info != 0"
allocate(phi_cav(dd_data % ncav), gradphi_cav(3, dd_data % ncav), &
    & psi(dd_data % nbasis, dd_data % nsph), force(3, dd_data % nsph))
call cpu_time(start_time)
call mkrhs(dd_data, phi_cav, gradphi_cav, psi)
call cpu_time(finish_time)
write(*, "(A,ES11.4E2,A)") "mkrhs time:", finish_time-start_time, " seconds"
call cpu_time(start_time)
call ddpcm(dd_data, phi_cav, gradphi_cav, psi, esolv, force)
call cpu_time(finish_time)
write(*, "(A,ES11.4E2,A)") "ddpcm time:", finish_time-start_time, " seconds"
write(*, "(A,ES25.16E3)") "ddpcm esolv:", esolv
write(*, *) "Full forces"
do isph = 1, dd_data % nsph
    write(6,'(1x,i5,3ES25.16E3)') isph, force(:,isph)
end do
deallocate(phi_cav, gradphi_cav, psi, force)
call ddfree(dd_data)

end program main

