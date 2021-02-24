!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/dd_pcm.f90
!! PCM model in the domain decomposition framework
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-24

!> Core routines and parameters of ddX software
module dd_pcm
use dd_core
use dd_operators
use dd_solvers
implicit none

contains

!> ddPCM solver
!!
!! Solves the problem within PCM model using a domain decomposition approach.
!!
!! @param[in] dd_data: ddX object with all input information
!! @param[in] phi_cav: Potential at cavity points
!! @param[in] gradphi_cav: Gradient of a potential at cavity points
!! @param[in] psi: TODO
!! @param[out] esolv: Solvation energy
!! @param[out] force: Analytical forces
subroutine ddpcm(dd_data, phi_cav, gradphi_cav, psi, esolv, force)
    ! Inputs:
    type(dd_data_type), intent(inout)  :: dd_data
    real(dp), intent(in) :: phi_cav(dd_data % ncav), &
        & gradphi_cav(3, dd_data % ncav), psi(dd_data % nbasis, dd_data % nsph)
    ! Outputs
    real(dp), intent(out) :: esolv, force(3, dd_data % nsph)
    ! Local variables
    real(dp), allocatable :: g(:, :), rhs(:, :), phieps(:, :), &
        & phiinf(:, :), xs(:, :)
    integer :: istatus, isph, niter
    real(dp) :: start_time, finish_time
    logical :: ok
    double precision, external :: ddot
    external :: dgemm
    ! Unwrap sparsely stored potential at cavity points and multiply by ui
    call wghpot(dd_data, phi_cav, dd_data % tmp_grid)
    ! Integrate against spherical harmonics and Lebedev weights to get Phi
    call dgemm('N', 'N', dd_data % nbasis, dd_data % nsph, dd_data % ngrid, &
        & one, dd_data % vwgrid, dd_data % vgrid_nbasis, dd_data % tmp_grid, &
        & dd_data % ngrid, zero, dd_data % phi, dd_data % nbasis)
    ! Compute Phi_infty
    call rinfx(dd_data, dd_data % phi, dd_data % phiinf)
    ! Set initial guess on Phi_epsilon as Phi
    dd_data % phieps = dd_data % phi
    ! Maximum number of iterations for an iterative solver
    niter = dd_data % maxiter
    ! Solve ddPCM system R_eps Phi_epsilon = Phi_infty
    call cpu_time(start_time)
    call jacobi_diis(dd_data, dd_data % n, dd_data % iprint, dd_data % ndiis, &
        & 4, dd_data % tol, dd_data % phiinf, dd_data % phieps, niter, ok, &
        & rx, apply_repsx_prec, hnorm)
    call cpu_time(finish_time)
    if (dd_data % iprint .ge. 1) then
        write(*, "(A,ES11.4E2,A)") " ddpcm step time:", &
            & finish_time-start_time, " seconds"
        write(*, "(A,I0)") " ddpcm step iterations: ", niter
    endif
    ! Solve ddCOSMO system L X = -Phi_epsilon with a zero initial guess
    niter = dd_data % maxiter
    dd_data % xs = zero
    call cpu_time(start_time)
    call jacobi_diis(dd_data, dd_data % n, dd_data % iprint, dd_data % ndiis, &
        & 4, dd_data % tol, dd_data % phieps, dd_data % xs, niter, ok, lx, &
        & ldm1x, hnorm)
    call cpu_time(finish_time)
    if (dd_data % iprint .ge. 1) then
        write(*, "(A,ES11.4E2,A)") " ddcosmo step time:", &
            & finish_time-start_time, " seconds"
        write(*, "(A,I0)") " ddcosmo step iterations: ", niter
    endif
    ! Solvation energy is computed
    esolv = pt5*ddot(dd_data % n, dd_data % xs, 1, psi, 1)
    ! Get forces if needed
    if (dd_data % force .eq. 1) then
        ! Solve adjoint ddCOSMO system
        niter = dd_data % maxiter
        dd_data % s = zero
        call cpu_time(start_time)
        call jacobi_diis(dd_data, dd_data % n, dd_data % iprint, &
            & dd_data % ndiis, 4, dd_data % tol, psi, dd_data % s, niter, ok, &
            & lstarx, ldm1x, hnorm)
        call cpu_time(finish_time)
        if (dd_data % iprint.ge.1) then
            write(*, "(A,ES11.4E2,A)") " adjoint ddcosmo step time:", &
                & finish_time-start_time, " seconds"
            write(*, "(A,I0)") " adjoint ddcosmo step iterations: ", niter
        endif
        ! Solve adjoint ddPCM system
        niter = dd_data % maxiter
        dd_data % y = zero
        call cpu_time(start_time)
        call jacobi_diis(dd_data, dd_data % n, dd_data % iprint, &
            & dd_data % ndiis, 4, dd_data % tol, dd_data % s, dd_data % y, &
            & niter, ok, rstarx, apply_rstarepsx_prec, hnorm)
        call cpu_time(finish_time)
        if (dd_data % iprint .ge. 1) then
            write(*,"(A,ES11.4E2,A)") " adjoint ddpcm step time:", &
                & finish_time-start_time, " seconds"
            write(*,"(A,I0)") " adjoint ddpcm step iterations: ", &
                & niter
        end if
    end if
end subroutine ddpcm

end module dd_pcm

