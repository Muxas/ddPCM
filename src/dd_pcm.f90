!> @copyright (c) 2020-2020 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/dd_pcm.f90
!! PCM model in the domain decomposition framework
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2020-12-17

!> Core routines and parameters of ddX software
module dd_pcm
use dd_core
use dd_operators
use dd_solvers
implicit none

contains

subroutine ddpcm(dd_data, phi, psi, tol, ndiis, maxiter)
    ! Inputs:
    type(dd_data_type), intent(in)  :: dd_data
    real(dp), intent(in) :: phi(dd_data % ncav), &
        & psi(dd_data % nbasis, dd_data % nsph)
    real(dp), intent(in) :: tol
    integer, intent(in) :: ndiis
    integer, intent(in) :: maxiter
    ! Local variables
    real(dp), allocatable :: g(:, :), rhs(:, :), phieps(:, :), &
        & phiinf(:, :), xs(:, :)
    integer :: istatus, isph, niter
    logical :: ok
    double precision, external :: ddot
    ! Accumulate right hand side
    allocate(g(dd_data % ngrid, dd_data % nsph), &
        & rhs(dd_data % nbasis, dd_data % nsph), stat=istatus)
    call wghpot(dd_data, phi, g)
    do isph = 1, dd_data % nsph
        call intrhs(dd_data % iprint, dd_data % ngrid, dd_data % lmax, &
            & dd_data % vwgrid, isph, g(:, isph), rhs(:, isph))
    end do
    call prtsph("rhs", dd_data % nbasis, dd_data % lmax, dd_data % nsph, 0, rhs)
    allocate(phiinf(dd_data % nbasis, dd_data % nsph))
    call rinfx(dd_data, rhs, phiinf)
    allocate(phieps(dd_data % nbasis, dd_data % nsph))
    phieps = rhs
    call prtsph("phiinf", dd_data % nbasis, dd_data % lmax, dd_data % nsph, 0, phiinf)
    niter = maxiter
    call jacobi_diis(dd_data, dd_data % n, dd_data % iprint, ndiis, 4, tol, &
        & phiinf, phieps, niter, ok, rx, apply_repsx_prec, hnorm)
    return
    allocate(xs(dd_data % nbasis, dd_data % nsph))
    call prtsph('phieps', dd_data % nbasis, dd_data % lmax, dd_data % nsph, 0, phieps)
    niter = maxiter
    xs = zero
    call jacobi_diis(dd_data, dd_data % n, dd_data % iprint, ndiis, 4, tol, &
        & phieps, xs, niter, ok, lx, ldm1x, hnorm)
    call prtsph('xs', dd_data % nbasis, dd_data % lmax, dd_data % nsph, 0, xs)
    write(*, *) "esolv=", pt5*ddot(dd_data % n, xs, 1, psi, 1)
    deallocate(g, rhs, xs, phieps, phiinf)
end subroutine ddpcm

end module dd_pcm

