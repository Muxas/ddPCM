!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/dd_pcm.f90
!! PCM model in the domain decomposition framework
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-11

!> Core routines and parameters of ddX software
module dd_pcm
use dd_core
use dd_operators
use dd_solvers
implicit none

contains

subroutine ddpcm(dd_data, phi, psi)
    ! Inputs:
    type(dd_data_type), intent(inout)  :: dd_data
    real(dp), intent(in) :: phi(dd_data % ncav), &
        & psi(dd_data % nbasis, dd_data % nsph)
    ! Local variables
    real(dp), allocatable :: g(:, :), rhs(:, :), phieps(:, :), &
        & phiinf(:, :), xs(:, :)
    integer :: istatus, isph, niter
    real(dp) :: start_time, finish_time
    logical :: ok
    double precision, external :: ddot
    ! Accumulate right hand side
    allocate(g(dd_data % ngrid, dd_data % nsph), &
        & rhs(dd_data % nbasis, dd_data % nsph), stat=istatus)
    call wghpot(dd_data, phi, g)
    do isph = 1, dd_data % nsph
        call intrhs(dd_data % iprint, dd_data % ngrid, dd_data % lmax, &
            & dd_data % vwgrid, dd_data % vgrid_nbasis, isph, g(:, isph), &
            & rhs(:, isph))
    end do
    allocate(phiinf(dd_data % nbasis, dd_data % nsph))
    call rinfx(dd_data, rhs, phiinf)
    allocate(phieps(dd_data % nbasis, dd_data % nsph))
    phieps = rhs
    niter = dd_data % maxiter
    call cpu_time(start_time)
    call jacobi_diis(dd_data, dd_data % n, dd_data % iprint, dd_data % ndiis, &
        & 4, dd_data % tol, phiinf, phieps, niter, ok, rx, apply_repsx_prec, &
        & hnorm)
    call cpu_time(finish_time)
    if (dd_data % iprint .ge. 1) then
        write(*, "(A,ES11.4E2,A)") " ddpcm step time:", &
            & finish_time-start_time, " seconds"
        write(*, "(A,I0)") " ddpcm step iterations: ", niter
    endif
    allocate(xs(dd_data % nbasis, dd_data % nsph))
    niter = dd_data % maxiter
    xs = zero
    call cpu_time(start_time)
    call jacobi_diis(dd_data, dd_data % n, dd_data % iprint, dd_data % ndiis, &
        & 4, dd_data % tol, phieps, xs, niter, ok, lx, ldm1x, hnorm)
    call cpu_time(finish_time)
    if (dd_data % iprint .ge. 1) then
        write(*, "(A,ES11.4E2,A)") " ddcosmo step time:", &
            & finish_time-start_time, " seconds"
        write(*, "(A,I0)") " ddcosmo step iterations: ", niter
    endif
    write(*, *) "esolv=", pt5*ddot(dd_data % n, xs, 1, psi, 1)
    deallocate(g, rhs, xs, phieps, phiinf)
end subroutine ddpcm

end module dd_pcm

