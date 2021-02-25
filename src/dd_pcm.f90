!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/dd_pcm.f90
!! PCM model in the domain decomposition framework
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-25

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
        & phiinf(:, :), xs(:, :), vsin(:), vcos(:), vplm(:), basloc(:), &
        & dbsloc(:, :), fx(:, :), ef(:, :)
    integer :: istatus, isph, niter, igrid, icav, inear, inode, jnear, jnode, &
        & jsph
    real(dp) :: start_time, finish_time, tmp1, tmp2, d(3), dnorm
    logical :: ok
    double precision, external :: ddot, dnrm2
    external :: dgemm
    ! Unwrap sparsely stored potential at cavity points and multiply by ui
    call wghpot(dd_data, phi_cav, dd_data % phi_grid, dd_data % tmp_grid)
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
        call dgemm('T', 'N', dd_data % ngrid, dd_data % nsph, &
            & dd_data % nbasis, one, dd_data % vgrid, dd_data % vgrid_nbasis, &
            & dd_data % s, dd_data % nbasis, zero, dd_data % sgrid, &
            & dd_data % ngrid)
        call dgemm('T', 'N', dd_data % ngrid, dd_data % nsph, &
            & dd_data % nbasis, one, dd_data % vgrid, dd_data % vgrid_nbasis, &
            & dd_data % y, dd_data % nbasis, zero, dd_data % ygrid, &
            & dd_data % ngrid)
        dd_data % g = dd_data % phi - dd_data % phieps
        dd_data % q = dd_data % s - fourpi/(dd_data % eps-one)*dd_data % y
        dd_data % qgrid = dd_data % sgrid - &
            & fourpi/(dd_data % eps-one)*dd_data % ygrid
        allocate(vsin(dd_data % lmax+1), vcos(dd_data % lmax+1), &
            & vplm(dd_data % nbasis), basloc(dd_data % nbasis), &
            & dbsloc(3, dd_data % nbasis), fx(3, dd_data % nsph), stat=istatus)
        if (istatus.ne.0) write(6,*) 'ddpcm forces allocation failed'
        call gradr(dd_data, force)
        do isph = 1, dd_data % nsph
            call fdoka(dd_data, isph, dd_data % xs, dd_data % sgrid(:, isph), &
                & basloc, dbsloc, vplm, vcos, vsin, force(:,isph)) 
            call fdokb(dd_data, isph, dd_data % xs, dd_data % sgrid, basloc, &
                & dbsloc, vplm, vcos, vsin, force(:, isph))
            call fdoga(dd_data, isph, dd_data % qgrid, dd_data % phi_grid, &
                & force(:, isph)) 
        end do
        force = -pt5 * force
        icav = 0
        do isph = 1, dd_data % nsph
            do igrid = 1, dd_data % ngrid
                if(dd_data % ui(igrid, isph) .ne. zero) then
                    icav = icav + 1
                    dd_data % zeta(icav) = -pt5 * dd_data % wgrid(igrid) * &
                        & dd_data % ui(igrid, isph) * ddot(dd_data % nbasis, &
                        & dd_data % vgrid(1, igrid), 1, &
                        & dd_data % q(1, isph), 1)
                    force(:, isph) = force(:, isph) + &
                        & dd_data % zeta(icav)*gradphi_cav(:, icav)
                end if
            end do
        end do
        !! Last term where we compute gradients of potential at centers of atoms
        !! spawned by intermediate zeta.
        allocate(ef(3, dd_data % nsph))
        if(dd_data % fmm .eq. 1) then
            !! This step can be substituted by a proper dgemm if zeta
            !! intermediate is converted from cavity points to all grid points
            !! with zero values at internal grid points
            ! P2M step
            icav = 0
            do isph = 1, dd_data % nsph
                inode = dd_data % snode(isph)
                dd_data % tmp_node_m(:, inode) = zero
                do igrid = 1, dd_data % ngrid
                    if(dd_data % ui(igrid, isph) .eq. zero) cycle
                    icav = icav + 1
                    call fmm_p2m(dd_data % cgrid(:, igrid), dd_data % zeta(icav), &
                        & one, dd_data % pm, dd_data % vscales, one, &
                        & dd_data % tmp_node_m(:, inode))
                end do
                dd_data % tmp_node_m(:, inode) = dd_data % tmp_node_m(:, inode) / &
                    & dd_data % rsph(isph)
            end do
            ! M2M, M2L and L2L translations
            if(dd_data % fmm_precompute .eq. 1) then
                call tree_m2m_reflection_use_mat(dd_data, dd_data % tmp_node_m)
                call tree_m2l_reflection_use_mat(dd_data, dd_data % tmp_node_m, &
                    & dd_data % tmp_node_l)
                call tree_l2l_reflection_use_mat(dd_data, dd_data % tmp_node_l)
            else
                call tree_m2m_rotation(dd_data, dd_data % tmp_node_m)
                call tree_m2l_rotation(dd_data, dd_data % tmp_node_m, &
                    & dd_data % tmp_node_l)
                call tree_l2l_rotation(dd_data, dd_data % tmp_node_l)
            end if
            ! Now compute near-field FMM gradients
            ! Cycle over all spheres
            icav = 0
            ef = zero
            do isph = 1, dd_data % nsph
                ! Cycle over all external grid points
                do igrid = 1, dd_data % ngrid
                    if(dd_data % ui(igrid, isph) .eq. zero) cycle
                    icav = icav + 1
                    ! Cycle over all near-field admissible pairs of spheres,
                    ! including pair (isph, isph) which is a self-interaction
                    inode = dd_data % snode(isph)
                    do jnear = dd_data % snear(inode), dd_data % snear(inode+1)-1
                        ! Near-field interactions are possible only between leaf
                        ! nodes, which must contain only a single input sphere
                        jnode = dd_data % near(jnear)
                        jsph = dd_data % order(dd_data % cluster(1, jnode))
                        d = dd_data % csph(:, isph) + &
                            & dd_data % cgrid(:, igrid)*dd_data % rsph(isph) - &
                            & dd_data % csph(:, jsph)
                        dnorm = dnrm2(3, d, 1)
                        ef(:, jsph) = ef(:, jsph) + &
                            & dd_data % zeta(icav)*d/(dnorm**3)
                    end do
                end do
            end do
            ! Take into account far-field FMM gradients
            tmp1 = one / sqrt(3d0) / dd_data % vscales(1)
            do isph = 1, dd_data % nsph
                inode = dd_data % snode(isph)
                tmp2 = tmp1 / dd_data % rsph(isph)
                ef(3, isph) = ef(3, isph) + tmp2*dd_data % tmp_node_l(3, inode)
                ef(1, isph) = ef(1, isph) + tmp2*dd_data % tmp_node_l(4, inode)
                ef(2, isph) = ef(2, isph) + tmp2*dd_data % tmp_node_l(2, inode)
            end do
            do isph = 1, dd_data % nsph
                force(:, isph) = force(:, isph) + ef(:, isph)*dd_data % charge(isph)
            end do
        ! Naive quadratically scaling implementation
        else
            ! This routines actually computes -grad, not grad
            call efld(dd_data % ncav, dd_data % zeta, dd_data % ccav, &
                & dd_data % nsph, dd_data % csph, ef)
            do isph = 1, dd_data % nsph
                force(:, isph) = force(:, isph) - ef(:, isph)*dd_data % charge(isph)
            end do
        end if
        deallocate(ef)
    end if
end subroutine ddpcm

end module dd_pcm

