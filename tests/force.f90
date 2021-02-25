!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/force.f90
!! Test of analytical forces against numerical
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
integer :: info, pmax=30
real(dp), allocatable :: phi_cav(:), gradphi_cav(:, :), psi(:, :), &
    & force(:, :), force_num(:, :)
real(dp) :: esolv1, esolv2, start_time, finish_time, step=0.0001, relerr
integer :: isph, i
real(dp), external :: dnrm2

! Read input file name
call getarg(1, fname)
write(*, *) "Using provided file ", trim(fname), " as a config file"
call ddfromfile(fname, dd_data, info)
if(info .ne. 0) stop "info != 0"
allocate(phi_cav(dd_data % ncav), gradphi_cav(3, dd_data % ncav), &
    & psi(dd_data % nbasis, dd_data % nsph), force(3, dd_data % nsph), &
    & force_num(3, dd_data % nsph))
call mkrhs(dd_data, phi_cav, gradphi_cav, psi)
call ddpcm(dd_data, phi_cav, gradphi_cav, psi, esolv1, force)
do isph = 1, dd_data % nsph
    do i = 1, 3
        dd_data % csph(i, isph) = dd_data % csph(i, isph) + step
        call ddpcm_solve(dd_data, esolv1, phi_cav, gradphi_cav, psi, force)
        dd_data % csph(i, isph) = dd_data % csph(i, isph) - two*step
        call ddpcm_solve(dd_data, esolv2, phi_cav, gradphi_cav, psi, force)
        dd_data % csph(i, isph) = dd_data % csph(i, isph) + step
        force_num(i, isph) = (esolv1-esolv2) / two / step
    end do
end do
relerr = dnrm2(3*dd_data % nsph, force_num-force, 1) / &
    & dnrm2(3*dd_data % nsph, force, 1)

write(6,'(2A60)') 'Analytical forces', 'Numerical forces'
do i = 1, dd_data % nsph
  write(6,'(6E20.10)') force(1,i), force(2,i), force(3,i), force_num(1,i), &
      & force_num(2,i), force_num(3,i)
end do

deallocate(phi_cav, gradphi_cav, psi, force, force_num)
call ddfree(dd_data)

write(*, *) "Rel.error of forces:", relerr
if (relerr .gt. 1d-6) stop 1
contains 

subroutine ddpcm_solve(dd_data, esolv, phi_cav, gradphi_cav, psi, force)
    type(dd_data_type), intent(inout) :: dd_data
    real(dp), intent(out) :: esolv, phi_cav(dd_data % ncav), &
        & gradphi_cav(3, dd_data % ncav), &
        & psi(dd_data % nbasis, dd_data % nsph), force(3, dd_data % nsph)
    type(dd_data_type) :: dd_data2
    call ddinit(dd_data % nsph, dd_data % charge, dd_data % csph(1, :), &
        & dd_data % csph(2, :), dd_data % csph(3, :), dd_data % rsph, &
        & dd_data % model, dd_data % lmax, dd_data % ngrid, 0, &
        & dd_data % fmm, dd_data % pm, dd_data % pl, &
        & dd_data % fmm_precompute, dd_data % iprint, dd_data % se, &
        & dd_data % eta, dd_data % eps, dd_data % kappa, &
        & dd_data % itersolver, dd_data % tol, dd_data % maxiter, &
        & dd_data % ndiis, dd_data % nproc, dd_data2, info)
    call mkrhs(dd_data2, phi_cav, gradphi_cav, psi)
    call ddpcm(dd_data2, phi_cav, gradphi_cav, psi, esolv, force)
    call ddfree(dd_data2)
end subroutine ddpcm_solve

end program main


