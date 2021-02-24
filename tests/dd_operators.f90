!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/dd_operators.f90
!! Tests for dd_operators module
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-24

program test_dd_operators
use dd_operators
implicit none

integer :: i, iprint=1, info, ngrid=590
! Some implementation do not work for alpha=1d-307, so range of test values is
! reduced. This can be fixed by enforcing proper order of multiplications like
! a*b*c into a*(b*c).
real(dp) :: alpha(4)=(/1d0, -1d0, 1d-300, 1d+300/)
type(dd_data_type) :: dd_data
integer, parameter :: nsph=10, lmax=7, force=1, itersolver=1, &
    & maxiter=1000, ndiis=25, nproc=1
real(dp), parameter :: se=0d0, eta=0.1d0, eps=78d0, kappa=0d0, tol=1d-7
real(dp) :: gcsph(3, nsph), csph(3, nsph), grsph(nsph), rsph(nsph), &
    & gcharge(nsph), charge(nsph)

gcsph(:, 1) = (/1d0, 1d0, 1d0/)
gcsph(:, 2) = (/2d0, 2d0, 2d0/)
gcsph(:, 3) = (/1d0, 1d0, 3d0/)
gcsph(:, 4) = (/2d0, 2d0, 4d0/)
gcsph(:, 5) = (/1d0, 1d0, 5d0/)
gcsph(:, 6) = (/2d0, 2d0, 6d0/)
gcsph(:, 7) = (/1d0, 1d0, 7d0/)
gcsph(:, 8) = (/2d0, 2d0, 8d0/)
gcsph(:, 9) = (/1d0, 1d0, 9d0/)
gcsph(:, 10) = (/2d0, 2d0, 10d0/)
grsph = (/1d-1, 2d-1, 3d-1, 4d-1, 5d-1, 6d-1, 7d-1, 8d-1, 9d-1, 1d0/)
gcharge = one

do i = 1, size(alpha)
    charge = abs(alpha(i)) * gcharge
    csph = alpha(i) * gcsph
    rsph = abs(alpha(i)) * grsph
    call ddinit(nsph, charge, csph(1, :), csph(2, :), csph(3, :), rsph, 2, &
        lmax, ngrid, force, 0, -1, -1, -1, 0, se, eta, eps, kappa, &
        & itersolver, tol, maxiter, ndiis, nproc, dd_data, info)
    if(info .ne. 0) stop 1
    !call check_dx(dd_data, lmax, lmax, iprint, 1d-4)
    !call check_dx(dd_data, 40, 40, iprint, 1d-15)
    call check_gradr(dd_data, lmax, lmax, iprint, 1d-4)
    call check_gradr(dd_data, lmax+1, lmax+1, iprint, 1d-4)
    call check_gradr(dd_data, 40, 40, iprint, 1d-15)
    call ddfree(dd_data)
end do

contains

subroutine check_dx(dd_data, pm, pl, iprint, threshold)
    ! Inputs
    type(dd_data_type), intent(inout) :: dd_data
    integer, intent(in) :: pm, pl, iprint
    real(dp), intent(in) :: threshold
    ! Local variables
    type(dd_data_type) :: dd_data_fmm, dd_data_fmm2
    integer :: info, irand, iseed(4)=(/0, 0, 0, 1/), do_diag
    integer, parameter :: nrand=10
    real(dp) :: x(dd_data % nbasis, dd_data % nsph, nrand), &
        & y(dd_data % nbasis, dd_data % nsph, nrand), &
        & z(dd_data % nbasis, dd_data % nsph, nrand), &
        & xx(nrand, nrand), yy(nrand, nrand), full_norm, diff_norm, &
        & forces(3, dd_data % nsph), forces2(3, dd_data % nsph)
    real(dp), external :: dnrm2
    ! Init FMM-related dd_data
    call ddinit(dd_data % nsph, dd_data % charge, dd_data % csph(1, :), &
        & dd_data % csph(2, :), dd_data % csph(3, :), dd_data % rsph, &
        & dd_data % model, dd_data % lmax, dd_data % ngrid, dd_data % force, &
        & 1, pm, pl, 0, dd_data % iprint, dd_data % se, dd_data % eta, &
        & dd_data % eps, dd_data % kappa, dd_data % itersolver, &
        & dd_data % tol, dd_data % maxiter, dd_data % ndiis, dd_data % nproc, &
        & dd_data_fmm, info)
    call ddinit(dd_data % nsph, dd_data % charge, dd_data % csph(1, :), &
        & dd_data % csph(2, :), dd_data % csph(3, :), dd_data % rsph, &
        & dd_data % model, dd_data % lmax, dd_data % ngrid, dd_data % force, &
        & 1, pm, pl, 1, dd_data % iprint, dd_data % se, dd_data % eta, &
        & dd_data % eps, dd_data % kappa, dd_data % itersolver, &
        & dd_data % tol, dd_data % maxiter, dd_data % ndiis, dd_data % nproc, &
        & dd_data_fmm2, info)
    ! Dense operator dx is trusted to have no errors, this must be somehow
    ! checked in the future.
    ! Get random x
    call dlarnv(2, iseed, dd_data % n * nrand, x)
    write(*, *) "pm=", pm, "pl=", pl
    do do_diag = 0, 1
        write(*, *) "do_diag=", do_diag
        ! Random check of FMM dx operator against dense dx operator
        do irand = 1, nrand
            call dx(dd_data, do_diag, x(:, :, irand), y(:, :, irand))
        end do
        full_norm = dnrm2(dd_data % n * nrand, y, 1)
        do irand = 1, nrand
        call dx(dd_data_fmm, do_diag, x(:, :, irand), z(:, :, irand))
        end do
        diff_norm = dnrm2(dd_data % n * nrand, y-z, 1)
        write(*, *) "dx_dense vs dx_fmm(no precompute) rel.error=", diff_norm/full_norm
        ! Random check FMM dx operators (with/out precomputed FMM translations)
        do irand = 1, nrand
            call dx(dd_data_fmm2, do_diag, x(:, :, irand), y(:, :, irand))
        end do
        diff_norm = dnrm2(dd_data % n * nrand, y-z, 1)
        write(*, *) "dx_fmm(no precompute) vs dx_fmm(precompute) rel.error=", diff_norm/full_norm
        ! Check dense adjoint operator dstarx
        do irand = 1, nrand
            call dstarx(dd_data, do_diag, x(:, :, irand), y(:, :, irand))
        end do
        call dgemm('T', 'N', nrand, nrand, dd_data % n, one, y, dd_data % n, &
            & y, dd_data % n, zero, xx, nrand)
        full_norm = dnrm2(nrand**2, xx, 1)
        do irand = 1, nrand
            call dx(dd_data, do_diag, y(:, :, irand), z(:, :, irand))
        end do
        call dgemm('T', 'N', nrand, nrand, dd_data % n, one, z, dd_data % n, &
            & x, dd_data % n, zero, yy, nrand)
        diff_norm = dnrm2(nrand**2, xx-yy, 1)
        write(*, *) "dstarx_dense vs dx_dense rel.error=", diff_norm/full_norm
        ! Check FMM adjoint operator dstarx (without precomputed FMM matrices)
        do irand = 1, nrand
            call dstarx(dd_data_fmm, do_diag, x(:, :, irand), y(:, :, irand))
        end do
        call dgemm('T', 'N', nrand, nrand, dd_data % n, one, y, dd_data % n, &
            & y, dd_data % n, zero, xx, nrand)
        full_norm = dnrm2(nrand**2, xx, 1)
        do irand = 1, nrand
            call dx(dd_data_fmm, do_diag, y(:, :, irand), z(:, :, irand))
        end do
        call dgemm('T', 'N', nrand, nrand, dd_data % n, one, z, dd_data % n, &
            & x, dd_data % n, zero, yy, nrand)
        diff_norm = dnrm2(nrand**2, xx-yy, 1)
        write(*, *) "dstarx_fmm vs dx_fmm (no precompute) rel.error=", diff_norm/full_norm
        ! Check FMM adjoint operator dstarx (with/out precomputed FMM matrices)
        do irand = 1, nrand
            call dstarx(dd_data_fmm2, do_diag, x(:, :, irand), z(:, :, irand))
        end do
        diff_norm = dnrm2(dd_data % n * nrand, y-z, 1)
        write(*, *) "dstarx_fmm (precompute) vs dstarx_fmm(no precompute) rel.error=", diff_norm/full_norm
    end do
    ! Free temporary objects
    call ddfree(dd_data_fmm)
    call ddfree(dd_data_fmm2)
end subroutine check_dx

subroutine check_gradr(dd_data, pm, pl, iprint, threshold)
    ! Inputs
    type(dd_data_type), intent(inout) :: dd_data
    integer, intent(in) :: pm, pl, iprint
    real(dp), intent(in) :: threshold
    ! Local variables
    type(dd_data_type) :: dd_data_fmm, dd_data_fmm2
    integer :: info, irand, iseed(4)=(/0, 0, 0, 1/), do_diag
    integer, parameter :: nrand=10
    real(dp) :: x(dd_data % nbasis, dd_data % nsph, nrand), &
        & y(dd_data % nbasis, dd_data % nsph, nrand), &
        & z(dd_data % nbasis, dd_data % nsph, nrand), &
        & xx(nrand, nrand), yy(nrand, nrand), full_norm, diff_norm, &
        & forces(3, dd_data % nsph), forces2(3, dd_data % nsph)
    real(dp), external :: dnrm2
    ! Init FMM-related dd_data
    call ddinit(dd_data % nsph, dd_data % charge, dd_data % csph(1, :), &
        & dd_data % csph(2, :), dd_data % csph(3, :), dd_data % rsph, &
        & dd_data % model, dd_data % lmax, dd_data % ngrid, dd_data % force, &
        & 1, pm, pl, 0, dd_data % iprint, dd_data % se, dd_data % eta, &
        & dd_data % eps, dd_data % kappa, dd_data % itersolver, &
        & dd_data % tol, dd_data % maxiter, dd_data % ndiis, dd_data % nproc, &
        & dd_data_fmm, info)
    call ddinit(dd_data % nsph, dd_data % charge, dd_data % csph(1, :), &
        & dd_data % csph(2, :), dd_data % csph(3, :), dd_data % rsph, &
        & dd_data % model, dd_data % lmax, dd_data % ngrid, dd_data % force, &
        & 1, pm, pl, 1, dd_data % iprint, dd_data % se, dd_data % eta, &
        & dd_data % eps, dd_data % kappa, dd_data % itersolver, &
        & dd_data % tol, dd_data % maxiter, dd_data % ndiis, dd_data % nproc, &
        & dd_data_fmm2, info)
    ! Dense operator dx is trusted to have no errors, this must be somehow
    ! checked in the future.
    ! Get random ygrid and g
    call dlarnv(2, iseed, dd_data % ngrid * dd_data % nsph, dd_data % ygrid)
    call dlarnv(2, iseed, dd_data % n, dd_data % g)
    ! Check gradr
    call gradr_dense(dd_data, forces)
    full_norm = dnrm2(3*dd_data % nsph, forces, 1)
    dd_data_fmm % ygrid = dd_data % ygrid
    dd_data_fmm % g = dd_data % g
    call gradr_fmm(dd_data_fmm, forces2)
    diff_norm = dnrm2(3*dd_data % nsph, forces-forces2, 1)
    write(*, *) "gradr dense vs fmm rel.error=", diff_norm / full_norm
    full_norm = dnrm2(3*dd_data % nsph, forces2, 1)
    dd_data_fmm2 % ygrid = dd_data % ygrid
    dd_data_fmm2 % g = dd_data % g
    call gradr_fmm(dd_data_fmm2, forces)
    diff_norm = dnrm2(3*dd_data % nsph, forces-forces2, 1)
    write(*, *) "gradr fmm (precompute) vs fmm (no precompute) rel.error=", &
        & diff_norm / full_norm
    ! Free temporary objects
    call ddfree(dd_data_fmm)
    call ddfree(dd_data_fmm2)
end subroutine check_gradr

end program

