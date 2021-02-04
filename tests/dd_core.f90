!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/dd_core.f90
!! Tests for dd_core module
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-02-04

program test_dd_core
use dd_core
implicit none

!character(len=255) :: testname, iprint_string
integer :: argc, iprint=1
integer :: p=10, i

! Input points and spheres
integer, parameter :: nsph = 10
real(dp) :: gsrc0(3), gsrc(3), gdst(3, 2), gcsph(3, nsph), grsph(nsph), &
    & gdst_csph(3, 2), gdst_rsph(2)
! Some implementation do not work for alpha=1d-307, so range of test values is
! reduced. This can be fixed by enforcing proper order of multiplications like
! a*b*c into a*(b*c).
real(dp) :: alpha(4)=(/1d0, -1d0, 1d-300, 1d+300/)

real(dp), external :: dnrm2
! Set inputs
gsrc0 = zero
gsrc = (/1d-1, 2d-1, 3d-1/)
! This test shall be aligned along OZ axis
gdst(:, 1) = (/-0.2d0, 0.1d0, 25.2d0/)
gdst_csph(:, 1) = (/0d0, 0d0, 25d0/)
gdst_rsph(1) = pt5
! This test shall NOT be aligned along OZ axis
gdst(:, 2) = (/-16.3d0, 20.2d0, 0.1d0/)
gdst_csph(:, 2) = (/-16d0, 20d0, 0d0/)
gdst_rsph(2) = pt5
gcsph(:, 1) = zero
grsph(1) = one
gcsph(:, 2) = (/0d0, 0d0, -1.1d0/)
grsph(2) = grsph(1) + dnrm2(3, gcsph(:, 2)-gcsph(:, 1), 1)
gcsph(:, 3) = (/0d0, 0d0, 7d-1/)
grsph(3) = grsph(1) + dnrm2(3, gcsph(:, 3)-gcsph(:, 1), 1)
gcsph(:, 4) = -pt5
grsph(4) = grsph(1) + dnrm2(3, gcsph(:, 4)-gcsph(:, 1), 1)
gcsph(:, 5) = pt5
grsph(5) = grsph(1) + dnrm2(3, gcsph(:, 4)-gcsph(:, 1), 1)
gcsph(:, 6) = zero
grsph(6) = three
gcsph(:, 7) = (/1d-1, 2d-1, 1.1d0/)
grsph(7) = grsph(1) + dnrm2(3, gcsph(:, 7)-gcsph(:, 1), 1)
gcsph(:, 8) = (/4d-1, 2d-9, 1.1d0/)
grsph(8) = grsph(1) + dnrm2(3, gcsph(:, 8)-gcsph(:, 1), 1)
gcsph(:, 9) = (/1.1d0, 0d0, 0d0/)
grsph(9) = grsph(1) + dnrm2(3, gcsph(:, 9)-gcsph(:, 1), 1)
gcsph(:, 10) = (/-4d-1, 1.1d0, 0d0/)
grsph(10) = grsph(1) + dnrm2(3, gcsph(:, 10)-gcsph(:, 1), 1)

! Check correctness of info for valid and invalid input parameters of ddinit
call check_ddinit_args()

! Check P2M and M2P operations of the FMM
do i = 1, size(alpha)
    call check_p2m_m2p(0, alpha(i), iprint, 6d-2)
    call check_p2m_m2p(1, alpha(i), iprint, 3d-3)
    call check_p2m_m2p(p, alpha(i), iprint, 120*epsilon(zero))
end do

! Check P2L and L2P operations of the FMM
do i = 1, size(alpha)
    call check_p2l_l2p(0, alpha(i), iprint, 6d-2)
    call check_p2l_l2p(1, alpha(i), iprint, 3d-3)
    call check_p2l_l2p(p, alpha(i), iprint, 120*epsilon(zero))
end do

! Check M2M operations of the FMM
do i = 1, size(alpha)
    call check_m2m(0, alpha(i), iprint, 10d0*epsilon(zero))
    call check_m2m(1, alpha(i), iprint, 10d0*epsilon(zero))
    call check_m2m(p, alpha(i), iprint, 10d0*epsilon(zero))
end do

! Check L2L operations of the FMM
do i = 1, size(alpha)
    call check_l2l(0, alpha(i), iprint, 10d0*epsilon(zero))
    call check_l2l(1, alpha(i), iprint, 10d0*epsilon(zero))
    call check_l2l(p, alpha(i), iprint, 40d0*epsilon(zero))
end do

! Check M2L operations of the FMM
do i = 1, size(alpha)
    call check_m2l(0, 0, alpha(i), iprint, 6d-2)
    call check_m2l(0, 1, alpha(i), iprint, 6d-2)
    call check_m2l(0, p, alpha(i), iprint, 6d-2)
    call check_m2l(1, 0, alpha(i), iprint, 2d-2)
    call check_m2l(1, 1, alpha(i), iprint, 3d-3)
    call check_m2l(1, p, alpha(i), iprint, 3d-3)
    call check_m2l(p, 0, alpha(i), iprint, 2d-2)
    call check_m2l(p, 1, alpha(i), iprint, 2d-4)
    call check_m2l(p, p, alpha(i), iprint, 20d0*epsilon(zero))
end do

! Check recursive inertial tree
do i = 1, size(alpha)
    call check_tree_rib(alpha(i))
end do

contains

subroutine check_ddinit_args()
    ! Example of correct args
    integer :: n=1, model=1, lmax=0, ngrid=0, force=1, fmm=1, pm=0, pl=0, &
        & iprint=0
    real(dp) :: x(10), y(10), z(10), rvdw(10), se=zero, eta=zero, eps=zero, &
        & kappa=zero
    type(dd_data_type) :: dd_data
    integer :: info=0, i, j
    real(dp) :: tmp
    ! Generate coordinates and radii
    rvdw = one
    do i = 1, 10
        x(i) = dble(2*i)
        y(i) = x(i)
        z(i) = x(i)
    end do
    ! Check correct input
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, se, eta, eps, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    ! Check different correct inputs with different n <= 10 (hardcoded value)
    do i = 1, 10
        call ddinit(i, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
            & pl, iprint, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect input n = 0
    i = 0
    call ddinit(i, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -1) stop 1
    call ddfree(dd_data)
    ! Check all possible models with other correct inputs
    do i = 1, 3
        call ddinit(n, x, y, z, rvdw, i, lmax, ngrid, force, fmm, pm, &
            & pl, iprint, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect models
    i = -1
    call ddinit(n, x, y, z, rvdw, i, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -6) stop 1
    call ddfree(dd_data)
    i = 4
    call ddinit(n, x, y, z, rvdw, i, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -6) stop 1
    call ddfree(dd_data)
    ! Check correct lmax
    do i = 1, 6
        call ddinit(n, x, y, z, rvdw, model, i, ngrid, force, fmm, pm, &
            & pl, iprint, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect lmax < 0
    i = -1
    call ddinit(n, x, y, z, rvdw, model, i, ngrid, force, fmm, pm, &
        & pl, iprint, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -7) stop 1
    call ddfree(dd_data)
    ! Check correct ngrid
    do i = 0, 1000, 100
        j = i
        call ddinit(n, x, y, z, rvdw, model, lmax, j, force, fmm, pm, &
            & pl, iprint, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect ngrid < 0
    i = -1
    call ddinit(n, x, y, z, rvdw, model, lmax, i, force, fmm, pm, &
        & pl, iprint, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -8) stop 1
    call ddfree(dd_data)
    ! Check correct force
    do i = 0, 1
        call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, i, fmm, pm, &
            & pl, iprint, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect force
    i = -1
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, i, fmm, pm, &
        & pl, iprint, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -9) stop 1
    call ddfree(dd_data)
    i = 2
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, i, fmm, pm, &
        & pl, iprint, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -9) stop 1
    call ddfree(dd_data)
    ! Check correct fmm
    do i = 0, 1
        call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, i, pm, &
            & pl, iprint, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect fmm
    i = -1
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, i, pm, &
        & pl, iprint, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -10) stop 1
    call ddfree(dd_data)
    i = 2
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, i, pm, &
        & pl, iprint, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -10) stop 1
    call ddfree(dd_data)
    ! Check correct pm (ignored if fmm=0)
    j = 0
    do i = -2, 2
        call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, j, i, &
            & pl, iprint, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check correct pm (fmm=1)
    j = 1
    do i = 0, 20, 5
        call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, j, i, &
            & pl, iprint, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect pm (fmm=1)
    j = 1
    i = -1
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, j, i, &
        & pl, iprint, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -11) stop 1
    call ddfree(dd_data)
    ! Check correct pl (ignored if fmm=0)
    j = 0
    do i = -2, 2
        call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, j, pm, &
            & i, iprint, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check correct pl (fmm=1)
    j = 1
    do i = 0, 20, 5
        call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, j, pm, &
            & i, iprint, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect pl (fmm=1)
    j = 1
    i = -1
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, j, pm, &
        & i, iprint, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -12) stop 1
    call ddfree(dd_data)
    ! Check correct iprint
    do i = 0, 10
        call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
            & pl, i, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect iprint
    i = -1
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, i, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -13) stop 1
    call ddfree(dd_data)
    ! Check correct se
    tmp = -one
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, tmp, eta, eps, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    tmp = zero
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, tmp, eta, eps, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    tmp = one
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, tmp, eta, eps, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    ! Check incorrect se
    tmp = 1.01d0
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, tmp, eta, eps, kappa, dd_data, info)
    if (info .ne. -14) stop 1
    call ddfree(dd_data)
    tmp = -1.01d0
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, tmp, eta, eps, kappa, dd_data, info)
    if (info .ne. -14) stop 1
    call ddfree(dd_data)
    ! Check correct eta
    tmp = zero
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, se, tmp, eps, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    tmp = pt5
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, se, tmp, eps, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    tmp = one
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, se, tmp, eps, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    ! Check incorrect eta
    tmp = 1.01d0
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, se, tmp, eps, kappa, dd_data, info)
    if (info .ne. -15) stop 1
    call ddfree(dd_data)
    tmp = -1d-2
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, se, tmp, eps, kappa, dd_data, info)
    if (info .ne. -15) stop 1
    call ddfree(dd_data)
    ! Check correct eps
    tmp = zero
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, se, eta, tmp, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    tmp = pt5
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, se, eta, tmp, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    tmp = one
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, se, eta, tmp, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    tmp = dble(1000)
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, se, eta, tmp, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    ! Check incorrect eps
    tmp = -1d-2
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, se, eta, tmp, kappa, dd_data, info)
    if (info .ne. -16) stop 1
    call ddfree(dd_data)
    ! Check incorrect kappa
    tmp = -1d-2
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, se, eta, eps, tmp, dd_data, info)
    if (info .ne. -17) stop 1
    call ddfree(dd_data)
end subroutine check_ddinit_args

! Check P2M and M2P for spherical harmonics
! List of explicitly checked functions:
!       fmm_p2m
!       fmm_m2p
! List of implicitly checked functions:
!       ylmscale
!       ylmbas
!       trgev
!       polleg
subroutine check_p2m_m2p(p, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: p, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
    integer :: nbasis, i, j, istatus
    real(dp) :: vscales((p+1)**2)
    real(dp) :: src0(3), src(3), dst(3), csph(3, nsph), rsph(nsph), v0, v(nsph)
    real(dp), dimension((p+1)**2, nsph) :: coef, coef2
    real(dp) :: ztranslate_mat((p+1)*(p+2)*(p+3)/6)
    real(dp) :: transform_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(dp) :: full_mat((p+1)**2, (p+1)**2)
    real(dp) :: z, r1(3, 3), dst1(3), full_norm, diff_norm
    logical :: ok
    real(dp), external :: dnrm2
    !! Scale inputs
    src0 = gsrc0
    src = alpha * gsrc
    dst = alpha * gdst(:, 1)
    csph = alpha * gcsph
    rsph = abs(alpha) * grsph
    ! Preliminaries
    nbasis = (p+1) * (p+1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(p, vscales)
    !! Check P2M parameters q and beta
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0,A,ES12.4E3)") "Check P2M q and beta params for p=", &
            & p, " alpha=", alpha
        write(*, "(A)") "================================"
    end if
    ! Check beta=zero
    coef(:, 1) = zero
    call fmm_p2m(src, one, alpha, p, vscales, zero, coef(:, 1))
    coef2(:, 1) = coef(:, 1)
    call fmm_p2m(src, one, alpha, p, vscales, zero, coef2(:, 1))
    v(1) = dnrm2(nbasis, coef(:, 1)-coef2(:, 1), 1) / &
        & dnrm2(nbasis, coef(:, 1), 1)
    if (iprint .gt. 0) then
        write(*, *) "P2M beta=zero param check: ", v(1) .eq. zero
    end if
    if (v(1) .ne. zero) stop 1
    coef2(:, 1) = coef(:, 1)
    ! Check non-zero beta
    call fmm_p2m(src, one, alpha, p, vscales, one, coef2(:, 1))
    v(1) = dnrm2(nbasis, two*coef(:, 1)-coef2(:, 1), 1) / &
        & dnrm2(nbasis, coef(:, 1), 1) / two
    if (iprint .gt. 0) then
        write(*, *) "P2M q and beta params check:", v(1) .eq. zero
    end if
    if (v(1) .ne. zero) stop 1
    call fmm_p2m(src, -pt5, alpha, p, vscales, -two, coef2(:, 1))
    v(1) = dnrm2(nbasis, 4.5d0*coef(:, 1)+coef2(:, 1), 1) / &
        & dnrm2(nbasis, coef(:, 1), 1) / 4.5d0
    if (iprint .gt. 0) then
        write(*, *) "P2M q and beta params check:", v(1) .eq. zero
    end if
    if (v(1) .ne. zero) stop 1
    !! Check M2P parameters alpha and beta
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,A,I0,A,ES12.4E3)") "Check M2P alpha and beta params ", &
            & "for p=", p, " alpha=", alpha
        write(*, "(A)") "================================"
    end if
    ! Check beta=zero
    v(1) = zero
    call fmm_m2p(dst, alpha, p, vscales, one, coef(:, 1), zero, v(1))
    v0 = v(1)
    call fmm_m2p(dst, alpha, p, vscales, one, coef(:, 1), zero, v(1))
    if (iprint .gt. 0) then
        write(*, *) "M2P beta=zero param check:", v0 .eq. v(1)
    end if
    if (v0 .ne. v(1)) stop 1
    ! Check alpha=zero
    call fmm_m2p(dst, alpha, p, vscales, zero, coef(:, 1), -one, v(1))
    if (iprint .gt. 0) then
        write(*, *) "M2P alpha=zero param check:", v0 .eq. -v(1)
    end if
    if (v0 .ne. -v(1)) stop 1
    ! Check non-zero alpha and beta
    v(1) = v0
    call fmm_m2p(dst, alpha, p, vscales, one, coef(:, 1), one, v(1))
    ok = abs(v(1)/v0-two) .le. 10d0*epsilon(zero)
    if (iprint .gt. 0) then
        write(*, *) "M2P alpha and beta params check:", ok
    end if
    if (.not. ok) stop 1
    call fmm_m2p(dst, alpha, p, vscales, -two, coef(:, 1), pt5, v(1))
    ok = abs(v(1)/v0+one) .le. 10d0*epsilon(zero)
    if (iprint .gt. 0) then
        write(*, *) "M2P alpha and beta parameter check:", ok
    end if
    if (.not. ok) stop 1
    !! Check M2P with a location of P equal to the center of M
    call fmm_m2p(src0, rsph(1), p, vscales, one, coef(:, 1), zero, v0)
    if (iprint .gt. 0) then
        write(*, *) "Check M2P with center of harmonics=particle location"
        write(*, *) "================================="
        write(*, *) "Required result must be 0"
        write(*, *) "Got", v0
    end if
    if (v0 .ne. zero) stop 1
    !! Check P2M+M2P with source particle and multipole harmonics are centered
    !! at the origin
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0,A,ES12.4E3)") "Check P2M and M2P for p=", p, &
            & " alpha=", alpha
        write(*, "(A)") "================================"
        write(*, *) "threshold=" , threshold
    end if
    call fmm_p2m(src0, one, alpha, p, vscales, zero, coef(:, 1))
    call fmm_m2p(dst, alpha, p, vscales, one, coef(:, 1), zero, v(1))
    v0 = one / dnrm2(3, dst, 1)
    v(1) = abs((v(1)-v0) / v0)
    if (iprint .gt. 0) then
        write(*, "(A,A,ES24.16E3)") "|| P2M(0) + M(0)2P - P2P ||  /  ", &
            & "|| P2P || =", v(1)
    end if
    if (v(1) .gt. threshold) stop 1
    call fmm_p2m(src0, pt5, alpha, p, vscales, pt5, coef(:, 1))
    call fmm_m2p(dst, alpha, p, vscales, one, coef(:, 1), zero, v(1))
    v(1) = abs((v(1)-v0) / v0)
    if (iprint .gt. 0) then
        write(*, "(A,A,ES24.16E3)") "|| P2M(0) + M(0)2P - P2P ||  /  ", &
            & "|| P2P || =", v(1)
    end if
    if (v(1) .gt. threshold) stop 1
    !! Check P2M+M2P for predefined spheres
    ! Compute multipole coefficients for given spheres from source particle
    ! Get potential, spawned by each sphere with its multipole expansion
    v0 = one / dnrm2(3, src-dst, 1)
    do i = 1, nsph
        call fmm_p2m(src-csph(:, i), one, rsph(i), p, vscales, zero, &
            & coef(:, i))
        call fmm_m2p(dst-csph(:, i), rsph(i), p, vscales, one, coef(:, i), &
            & zero, v(i))
        ! Finally check p2m+m2p
        v(i) = abs((v(i)-v0) / v0)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,ES24.16E3)") &
                & "|| P2M(", i, ") + M(", i, &
                & ")2P - P2P ||  /  || P2P || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
    end do
end subroutine check_p2m_m2p

! Check transformation of spherical harmonics
! List of explicitly checked functions:
!       coord_reflect_get_mat
!       fmm_sph_transform
!       fmm_sph_transform_get_mat
!       fmm_sph_transform_use_mat
! TODO explicitly:
!       fmm_sph_rotate_oz
!       fmm_sph_transform_oxz
!subroutine check_sph_transform(p, alpha, iprint, threshold)
!    !! Check reflection of coordinates
!    if (iprint .gt. 0) then
!        write(*, *)
!        write(*, "(A)") "Check coord_reflect_get_mat"
!        write(*, "(A)") "================================"
!        write(*, *) "threshold=", threshold
!    end if
!    do i = 1, nsph
!        ! Get reflection matrix from new coordinates (y, z, x) to old (y, z, x)
!        call coord_reflect_get_mat(dst-csph(:, i), z, r1)
!        ! Get where vector is actually reflected
!        dst1 = zero
!        do j = 1, 3
!            dst1(j) = (dst(1)-csph(1, i))*r1(3, j) + &
!                & (dst(2)-csph(2, i))*r1(1, j) + &
!                & (dst(3)-csph(3, i))*r1(2, j)
!        end do
!        ! Convert (y,z,x) to (x,y,z)
!        v(1) = dst1(1)
!        dst1(1) = dst1(3)
!        dst1(3) = dst1(2)
!        dst1(2) = v(1)
!        if (iprint .gt. 0) then
!            write(*, "(A)") "Reflected vector (x,y,z) must be aligned along OZ"
!            write(*, "(3ES25.16E3)") dst1(1), dst1(2), dst1(3)
!        end if
!        ! Check that x and y coordinates are zeros compared to z coordinate
!        v(1) = abs(dst1(1) / dst1(3))
!        v(2) = abs(dst1(2) / dst1(3))
!        if (max(v(1), v(2)) .gt. threshold) stop 1
!        ! Check that dst1(3) is equal to z output of coord_reflect_get_mat
!        if (abs((z-dst1(3)) / z) .gt. threshold) stop 1
!    end do
!    !! Check transformation of spherical harmonics under reflection of system
!    !! of coordinates
!    if (iprint .gt. 0) then
!        write(*, *)
!        write(*, "(A,I0)") "Check fmm_sph_transform for p=", p
!        write(*, "(A)") "================================"
!        write(*, *) "threshold=", threshold
!    end if
!    do i = 1, nsph
!        ! Get reflection of coordinates
!        call coord_reflect_get_mat(dst-csph(:, i), z, r1)
!        ! Apply reflection of spherical harmonics
!        call fmm_sph_transform(p, r1, one, msph(:, i), zero, m2sph(:, i))
!        ! Compute corresponding potential with change of coordinates
!        dst1(1:2) = zero
!        dst1(3) = z
!        call fmm_m2p(dst1, rsph(i), p, vscales, one, m2sph(:, i), zero, v(i))
!        ! Get reference value of potential (without change of coordinates)
!        call fmm_m2p(dst-csph(:, i), rsph(i), p, vscales, one, msph(:, i), &
!            & zero, v0)
!        if (iprint .gt. 0) then
!            write(*, "(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES24.16E3)") &
!                & "|| M(", i, ")2M(", i, "') + M(", i, "')2P - M(", i, &
!                & ")2P ||  /  || M(", i, ")2P || =", abs((v(i)-v0)/v0)
!        end if
!        if (abs((v(i)-v0)/v0) .gt. threshold) stop 1
!        ! Repeat with different alpha and beta
!        call fmm_sph_transform(p, r1, -one, msph(:, i), two, m2sph(:, i))
!        call fmm_m2p(dst1, rsph(i), p, vscales, one, m2sph(:, i), zero, v(i))
!        if (iprint .gt. 0) then
!            write(*, "(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES24.16E3)") &
!                & "|| M(", i, ")2M(", i, "') + M(", i, "')2P - M(", i, &
!                & ")2P ||  /  || M(", i, ")2P || =", abs((v(i)-v0)/v0)
!        end if
!        if (abs((v(i)-v0)/v0) .gt. threshold) stop 1
!    end do
!    !! Check transformation of spherical harmonics under reflection of system
!    !! of coordinates
!    if (iprint .gt. 0) then
!        write(*, *)
!        write(*, "(A,I0)") "Check fmm_sph_transform_get_mat for p=", p
!        write(*, "(A,I0)") "Check fmm_sph_transform_use_mat for p=", p
!        write(*, "(A)") "================================"
!        write(*, *) "threshold=", threshold
!    end if
!    do i = 1, nsph
!        ! Get reflection of coordinates
!        call coord_reflect_get_mat(dst-csph(:, i), z, r1)
!        ! Get matrix of reflection of spherical harmonics
!        call fmm_sph_transform_get_mat(p, r1, transform_mat)
!        ! Apply reflection of spherical harmonics
!        call fmm_sph_transform_use_mat(p, transform_mat, one, msph(:, i), &
!            & zero, m2sph(:, i))
!        ! Compute corresponding potential with change of coordinates
!        dst1(1:2) = zero
!        dst1(3) = z
!        call fmm_m2p(dst1, rsph(i), p, vscales, one, m2sph(:, i), zero, v(i))
!        ! Get reference value of potential (without change of coordinates)
!        call fmm_m2p(dst-csph(:, i), rsph(i), p, vscales, one, msph(:, i), &
!            & zero, v0)
!        if (iprint .gt. 0) then
!            write(*, "(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES24.16E3)") &
!                & "|| M(", i, ")2M(", i, "') + M(", i, "')2P - M(", i, &
!                & ")2P ||  /  || M(", i, ")2P || =", abs((v(i)-v0)/v0)
!        end if
!        if (abs((v(i)-v0)/v0) .gt. threshold) stop 1
!        ! Repeat with different alpha and beta
!        call fmm_sph_transform_use_mat(p, transform_mat, -one, msph(:, i), &
!            & two, m2sph(:, i))
!        call fmm_m2p(dst1, rsph(i), p, vscales, one, m2sph(:, i), zero, v(i))
!        if (iprint .gt. 0) then
!            write(*, "(A,I0,A,I0,A,I0,A,I0,A,I0,A,ES24.16E3)") &
!                & "|| M(", i, ")2M(", i, "') + M(", i, "')2P - M(", i, &
!                & ")2P ||  /  || M(", i, ")2P || =", abs((v(i)-v0)/v0)
!        end if
!        if (abs((v(i)-v0)/v0) .gt. threshold) stop 1
!    end do
!end subroutine check_sph_transform

! Check P2L and L2P for spherical harmonics
! List of explicitly checked functions:
!       fmm_p2l
!       fmm_l2p
subroutine check_p2l_l2p(p, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: p, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
    integer :: nbasis, i, j
    real(dp) :: delta = 10d0 * epsilon(one)
    real(dp) :: vscales((p+1)**2)
    real(dp) :: dst0(3), dst(3), src(3), csph(3, nsph), rsph(nsph), v0, v(nsph)
    real(dp), dimension((p+1)**2, nsph) :: coef, coef2
    logical :: ok
    real(dp), external :: dnrm2
    ! Scale inputs
    dst0 = gsrc0
    dst = alpha * gsrc
    src = alpha * gdst(:, 1)
    csph = alpha * gcsph
    rsph = alpha * grsph
    ! Preliminaries
    nbasis = (p+1) * (p+1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(p, vscales)
    !! Check P2L parameters q and beta
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0,A,ES12.4E3)") "Check P2L q and beta params for p=", &
            & p, " alpha=", alpha
        write(*, "(A)") "================================"
    end if
    ! Check beta=zero
    coef(:, 1) = zero
    call fmm_p2l(src, one, alpha, p, vscales, zero, coef(:, 1))
    coef2(:, 1) = coef(:, 1)
    call fmm_p2l(src, one, alpha, p, vscales, zero, coef2(:, 1))
    v(1) = dnrm2(nbasis, coef(:, 1)-coef2(:, 1), 1) / &
        & dnrm2(nbasis, coef(:, 1), 1)
    ok = v(1) .lt. delta
    if (iprint .gt. 0) then
        write(*, *) "P2L beta=zero param check: ", ok
    end if
    if (.not. ok) stop 1
    ! Check non-zero beta
    coef2(:, 1) = coef(:, 1)
    call fmm_p2l(src, one, alpha, p, vscales, one, coef2(:, 1))
    v(1) = dnrm2(nbasis, two*coef(:, 1)-coef2(:, 1), 1) / &
        & dnrm2(nbasis, coef(:, 1), 1) / two
    ok = v(1) .lt. delta
    if (iprint .gt. 0) then
        write(*, *) "P2L q and beta params check:", ok
    end if
    if (.not. ok) stop 1
    call fmm_p2l(src, -pt5, alpha, p, vscales, -two, coef2(:, 1))
    v(1) = dnrm2(nbasis, 4.5d0*coef(:, 1)+coef2(:, 1), 1) / &
        & dnrm2(nbasis, coef(:, 1), 1) / 4.5d0
    ok = v(1) .lt. delta
    if (iprint .gt. 0) then
        write(*, *) "P2L q and beta params check:", ok
    end if
    if (.not. ok) stop 1
    !! Check L2P parameters alpha and beta
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,A,I0,A,ES12.4E3)") "Check L2P alpha and beta params ", &
            & "for p=", p, " alpha=", alpha
        write(*, "(A)") "================================"
    end if
    ! Check beta=zero
    v(1) = zero
    call fmm_l2p(dst, alpha, p, vscales, one, coef(:, 1), zero, v(1))
    v0 = v(1)
    call fmm_l2p(dst, alpha, p, vscales, one, coef(:, 1), zero, v(1))
    ok = v0 .eq. v(1)
    if (iprint .gt. 0) then
        write(*, *) "L2P beta=zero param check:", ok
    end if
    if (.not. ok) stop 1
    ! Check alpha=zero
    call fmm_l2p(dst, alpha, p, vscales, zero, coef(:, 1), -one, v(1))
    ok = v0 .eq. -v(1)
    if (iprint .gt. 0) then
        write(*, *) "L2P alpha=zero param check:", ok
    end if
    if (.not. ok) stop 1
    ! Check non-zero alpha and beta
    v(1) = v0
    call fmm_l2p(dst, alpha, p, vscales, one, coef(:, 1), one, v(1))
    ok = abs(v(1)/v0-two) .le. delta
    if (iprint .gt. 0) then
        write(*, *) "L2P alpha and beta params check:", ok
    end if
    if (.not. ok) stop 1
    call fmm_l2p(dst, alpha, p, vscales, -two, coef(:, 1), pt5, v(1))
    ok = abs(v(1)/v0+one) .le. delta
    if (iprint .gt. 0) then
        write(*, *) "L2P alpha and beta parameter check:", ok
    end if
    if (.not. ok) stop 1
    !! Check P2L with a location of P equal to the center of L
    if (iprint .gt. 0) then
        write(*, *) "Check P2L with center of harmonics=particle location"
        write(*, *) "================================="
        write(*, *) "Required norm of resulting vector must be 0"
    end if
    call fmm_p2l(dst0, one, rsph(1), p, vscales, zero, coef(:, 1))
    v0 = dnrm2((p+1)**2, coef(:, 1), 1)
    ok = v0 .eq. zero
    if (iprint .gt. 0) then
        write(*, *) "|| diff || =", v0
    end if
    if (.not. ok) stop 1
    coef(:, 1) = one
    call fmm_p2l(dst0, one, rsph(1), p, vscales, -one, coef(:, 1))
    v0 = dnrm2((p+1)**2, coef(:, 1)+one, 1)
    ok = v0 .eq. zero
    if (iprint .gt. 0) then
        write(*, *) "|| diff || =", v0
    end if
    if (.not. ok) stop 1
    !! Check P2L+L2P with source particle at the origin and local harmonics
    !! and target particle are at the same point
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0,A,ES12.4E3)") "Check P2L and L2P for p=", p, &
            & " alpha=", alpha
        write(*, "(A)") "================================"
        write(*, *) "threshold=" , threshold
    end if
    call fmm_p2l(src, one, alpha, p, vscales, zero, coef(:, 1))
    call fmm_l2p(dst0, alpha, p, vscales, one, coef(:, 1), zero, v(1))
    v0 = one / dnrm2(3, src, 1)
    v(1) = abs((v(1)-v0) / v0)
    if (iprint .gt. 0) then
        write(*, "(A,A,ES24.16E3)") "|| P2L(0) + L(0)2P - P2P ||  /  ", &
            & "|| P2P || =", v(1)
    end if
    if (v(1) .gt. threshold) stop 1
    call fmm_p2l(src, one, alpha, p, vscales, zero, coef(:, 1))
    call fmm_l2p(dst0, alpha, p, vscales, one, coef(:, 1), zero, v(1))
    v(1) = abs((v(1)-v0) / v0)
    if (iprint .gt. 0) then
        write(*, "(A,A,ES24.16E3)") "|| P2L(0) + L(0)2P - P2P ||  /  ", &
            & "|| P2P || =", v(1)
    end if
    if (v(1) .gt. threshold) stop 1
    !! Check P2L+L2P for predefined spheres
    ! Get potential, spawned by each sphere with its multipole expansion
    v0 = one / dnrm2(3, src-dst, 1)
    ! Compute multipole coefficients for given spheres from source particle
    do i = 1, nsph
        call fmm_p2l(src-csph(:, i), one, rsph(i), p, vscales, zero, &
            & coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales, one, coef(:, i), &
            & zero, v(i))
        ! Finally check p2m+m2p
        v(i) = abs((v(i)-v0) / v0)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,ES24.16E3)") &
                & "|| P2L(", i, ") + L(", i, &
                & ")2P - P2P ||  /  || P2P || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
    end do
end subroutine check_p2l_l2p

! Check M2M for spherical harmonics
! List of explicitly checked functions:
!       fmm_m2m_ztranslate
!       fmm_m2m_ztranslate_adj
!       fmm_m2m_ztranslate_get_mat
!       fmm_m2m_ztranslate_use_mat
!       fmm_m2m_ztranslate_use_mat_adj
!       fmm_m2m_reflection
!       fmm_m2m_reflection2
!       fmm_m2m_rotation
!       fmm_m2m_reflection_adj
!       fmm_m2m_reflection2_adj
!       fmm_m2m_rotation_adj
!       fmm_m2m_reflection_get_mat
!       fmm_m2m_reflection_use_mat
!       fmm_m2m_reflection_use_mat_adj
subroutine check_m2m(p, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: p, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
    integer :: nbasis, i, j
    real(dp) :: vscales((p+1)**2), vfact(2*p+1)
    real(dp) :: src0(3), src(3), dst(3), csph(3, nsph), rsph(nsph), v0, v(nsph)
    real(dp), dimension((p+1)**2, nsph) :: coef, coef2
    real(dp) :: ztranslate_mat((p+1)*(p+2)*(p+3)/6)
    real(dp) :: transform_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(dp) :: full_mat((p+1)**2, (p+1)**2)
    real(dp) :: z, r1(3, 3), dst1(3), full_norm, diff_norm
    logical :: ok
    real(dp), external :: dnrm2
    !! Scale inputs
    src0 = gsrc0
    src = alpha * gsrc
    dst = alpha * gdst(:, 1)
    csph = alpha * gcsph
    rsph = alpha * grsph
    !! Preliminaries
    nbasis = (p+1) * (p+1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(p, vscales)
    vfact(1) = one
    do i = 2, 2*p+1
        vfact(i) = vfact(i-1) * sqrt(dble(i-1))
    end do
    !! Get multipole coefficients of all spheres by P2M
    do i = 1, nsph
        call fmm_p2m(src-csph(:, i), one, rsph(i), p, vscales, zero, &
            & coef(:, i))
    end do
    !! Check M2M OZ translation by fmm_m2m_ztranslate. Spheres 2 and 3 are
    !! explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2m_ztranslate"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || P2M(1) + M(1)2M(i) - P2M(i) || / ", &
            & "|| P2M(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        call fmm_m2m_ztranslate(-csph(3, i), rsph(1), rsph(i), p, vscales, &
            & vfact, one, coef(:, 1), zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_ztranslate(-csph(3, i), rsph(1), rsph(i), p, vscales, &
            & vfact, -one, coef(:, 1), two, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    ! Check M2M with zero OZ translation. Sphere 6 is explicitly set for this
    call fmm_m2m_ztranslate(zero, rsph(1), rsph(6), p, vscales, vfact, one, &
        & coef(:, 1), zero, coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
    end if
    if (.not. ok) stop 1
    call fmm_m2m_ztranslate(zero, rsph(1), rsph(6), p, vscales, vfact, -one, &
        & coef(:, 1), two, coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check M2M OZ translation by precomputed matrices. Spheres 2 and 3 are
    !! explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2m_ztranslate_get_mat"
        write(*, "(A)") "Check fmm_m2m_ztranslate_use_mat"
        write(*, "(A)") "Check fmm_m2m_scale"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || P2M(1) + M(1)2M(i) - P2M(i) || / ", &
            & "|| P2M(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        call fmm_m2m_ztranslate_get_mat(-csph(3, i), rsph(1), rsph(i), p, &
            & vscales, vfact, ztranslate_mat)
        call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, one, coef(:, 1), &
            & zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, -one, coef(:, 1), &
            & two, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    ! M2M with 0 shift is simply scaling
    call fmm_m2m_scale(rsph(1), rsph(6), p, one, coef(:, 1), zero, coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
    end if
    if (.not. ok) stop 1
    call fmm_m2m_scale(rsph(1), rsph(6), p, -one, coef(:, 1), two, &
        & coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check adjoint M2M OZ translation by precomputed matrices of direct M2M
    !! OZ translation. Spheres 2 and 3 are explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2m_ztranslate_use_mat_adj"
        write(*, "(A)") "Check fmm_m2m_scale_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || M(1)2M(i) - [adj M(1)2M(i)]^T || ", &
            & "/ || M(1)2M(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        ! Generate entire matrix of a direct M2M OZ translation
        call fmm_m2m_ztranslate_get_mat(-csph(3, i), rsph(1), rsph(i), p, &
            & vscales, vfact, ztranslate_mat)
        do j = 1, (p+1)**2
            coef2(:, i) = zero
            coef2(j, i) = one
            call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, one, &
                & coef2(:, i), zero, full_mat(:, j))
        end do
        full_norm = dnrm2((p+1)**4, full_mat, 1)
        ! Subtract transpose of an adjoint M2M OZ translation matrix
        do j = 1, (p+1)**2
            coef2(:, i) = zero
            coef2(j, i) = one
            call fmm_m2m_ztranslate_use_mat_adj(p, ztranslate_mat, -two, &
                & coef2(:, i), two, full_mat(j, :))
        end do
        diff_norm = dnrm2((p+1)**4, full_mat, 1)
        v(i) = diff_norm / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_ztranslate_use_mat_adj(p, ztranslate_mat, one, &
            & coef(:, 1), zero, full_mat(:, 1))
        full_norm = dnrm2((p+1)**2, full_mat(:, 1), 1)
        call fmm_m2m_ztranslate_use_mat_adj(p, ztranslate_mat, -two, &
            & coef(:, 1), two, full_mat(:, 1))
        diff_norm = dnrm2((p+1)**2, full_mat(:, 1), 1)
        v(i) = diff_norm / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    ! Check fmm_m2m_scale_adj. Sphere 6 is intended for this. Result is the
    ! same as with direct fmm_m2m_scale.
    call fmm_m2m_scale_adj(rsph(6), rsph(1), p, one, coef(:, 1), zero, &
        & coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
    end if
    if (.not. ok) stop 1
    call fmm_m2m_scale_adj(rsph(6), rsph(1), p, -one, coef(:, 1), two, &
        & coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check adjoint M2M OZ translation without precomputed matrices. Spheres 2
    !! and 3 are explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2m_ztranslate_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || [adj M(1)2M(i)] - [adj M(1)2M(i)]", &
            & " || / || [adj M(1)2M(i)] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        call fmm_m2m_ztranslate_get_mat(-csph(3, i), rsph(1), rsph(i), p, &
            & vscales, vfact, ztranslate_mat)
        call fmm_m2m_ztranslate_adj(csph(3, i), rsph(i), rsph(1), p, &
            & vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        full_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        call fmm_m2m_ztranslate_use_mat_adj(p, ztranslate_mat, -one, &
            & coef(:, 1), one, coef2(:, i))
        diff_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        v(i) = diff_norm / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_ztranslate_use_mat_adj(p, ztranslate_mat, one, &
            & coef(:, 1), zero, coef2(:, i))
        full_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        call fmm_m2m_ztranslate_adj(csph(3, i), rsph(i), rsph(1), p, &
            & vscales, vfact, -two, coef(:, 1), two, coef2(:, i))
        diff_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        v(i) = diff_norm / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    !! Check that attempt to use M2M matrix with 0 shift raises NaNs
    call fmm_m2m_ztranslate_get_mat(zero, rsph(1), rsph(6), p, vscales, &
        & vfact, ztranslate_mat)
    call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, one, coef(:, 1), zero, &
        & coef2(:, 6))
    ! Adjoint with zero translation is the same as direct with zero translation
    call fmm_m2m_ztranslate_adj(zero, rsph(6), rsph(1), p, vscales, vfact, &
        & one, coef(:, 1), zero, coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
    end if
    if (.not. ok) stop 1
    call fmm_m2m_ztranslate_adj(zero, rsph(6), rsph(1), p, vscales, vfact, &
        & -one, coef(:, 1), two, coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    ok = v(6) .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v(6)
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check M2M by reflection
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2m_reflection"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || P2M(1) + M(1)2M(i) - P2M(i) || / ", &
            & "|| P2M(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, nsph
        call fmm_m2m_reflection(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_reflection(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, -one, coef(:, 1), two, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    if(iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check M2M by reflection with help of matrices
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_m2m_reflection_get_mat"
        write(*, "(A,I0)") "Check fmm_m2m_reflection_use_mat"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || P2M(1) + M(1)2M(i) - P2M(i) || / ", &
            & "|| P2M(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, nsph
        call fmm_m2m_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, vfact, transform_mat, ztranslate_mat)
        call fmm_m2m_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, transform_mat, ztranslate_mat, one, coef(:, 1), &
            & zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, transform_mat, ztranslate_mat, -one, coef(:, 1), &
            & two, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    if(iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check adjoint M2M by reflection with help of matrices
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_m2m_reflection_use_mat_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || M(1)2M(i) - [adj M(1)2M(i)]^T || ", &
            & "/ || M(1)2M(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, nsph
        ! Generate entire matrix of a direct M2M translation
        call fmm_m2m_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, vfact, transform_mat, ztranslate_mat)
        do j = 1, (p+1)**2
            coef2(:, i) = zero
            coef2(j, i) = one
            call fmm_m2m_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
                & rsph(i), p, transform_mat, ztranslate_mat, one, &
                & coef2(:, i), zero, full_mat(:, j))
        end do
        full_norm = dnrm2((p+1)**4, full_mat, 1)
        ! Subtract transpose of an adjoint M2M translation matrix
        do j = 1, (p+1)**2
            coef2(:, i) = zero
            coef2(j, i) = one
            call fmm_m2m_reflection_use_mat_adj(csph(:, i)-csph(:, 1), &
                & rsph(i), rsph(1), p, transform_mat, ztranslate_mat, -two, &
                & coef2(:, i), two, full_mat(j, :))
        end do
        diff_norm = dnrm2((p+1)**4, full_mat, 1)
        v(i) = diff_norm / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_reflection_use_mat_adj(csph(:, i)-csph(:, 1), rsph(i), &
            & rsph(1), p, transform_mat, ztranslate_mat, one, coef(:, i), &
            & zero, full_mat(:, 1))
        full_norm = dnrm2((p+1)**2, full_mat(:, 1), 1)
        call fmm_m2m_reflection_use_mat_adj(csph(:, i)-csph(:, 1), rsph(i), &
            & rsph(1), p, transform_mat, ztranslate_mat, -two, coef(:, i), &
            & two, full_mat(:, 1))
        diff_norm = dnrm2((p+1)**2, full_mat(:, 1), 1)
        v(i) = diff_norm / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    if(iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check adjoint M2M by reflection
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_m2m_reflection_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || [adj M(1)2M(i)] - [adj M(1)2M(i)]", &
            & " || / || [adj M(1)2M(i)] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, nsph
        ! Check with help of fmm_m2m_reflection_use_mat that is checked already
        call fmm_m2m_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, vfact, transform_mat, ztranslate_mat)
        ! Check with beta=zero
        call fmm_m2m_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        ! Subtract reference value
        call fmm_m2m_reflection_use_mat_adj(csph(:, i)-csph(:, 1), &
            & rsph(i), rsph(1), p, transform_mat, ztranslate_mat, -one, &
            & coef(:, 1), one, coef2(:, i))
        v(i) = dnrm2((p+1)**2, coef2(:, i), 1) / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        ! Get reference value
        call fmm_m2m_reflection_use_mat_adj(csph(:, i)-csph(:, 1), &
            & rsph(i), rsph(1), p, transform_mat, ztranslate_mat, one, &
            & coef(:, 1), zero, coef2(:, i))
        full_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        ! Subtract value to be checked
        call fmm_m2m_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -two, coef(:, 1), two, coef2(:, i))
        v(i) = dnrm2((p+1)**2, coef2(:, i), 1) / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    if(iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check M2M by reflection (another implementation)
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2m_reflection2"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || P2M(1) + M(1)2M(i) - P2M(i) || / ", &
            & "|| P2M(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, nsph
        call fmm_m2m_reflection2(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_reflection2(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, -one, coef(:, 1), two, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    if(iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check adjoint M2M by reflection (another implementation)
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_m2m_reflection2_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || [adj M(1)2M(i)] - [adj M(1)2M(i)]", &
            & " || / || [adj M(1)2M(i)] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, nsph
        ! Check with beta=zero
        call fmm_m2m_reflection2_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        ! Subtract reference value
        call fmm_m2m_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -one, coef(:, 1), one, coef2(:, i))
        v(i) = dnrm2((p+1)**2, coef2(:, i), 1) / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        ! Get reference value
        call fmm_m2m_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        full_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        ! Subtract value to be checked
        call fmm_m2m_reflection2_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -two, coef(:, 1), two, coef2(:, i))
        v(i) = dnrm2((p+1)**2, coef2(:, i), 1) / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    if(iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check M2M by rotation
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2m_rotation"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || P2M(1) + M(1)2M(i) - P2M(i) || / ", &
            & "|| P2M(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, nsph
        call fmm_m2m_rotation(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        call fmm_m2m_rotation(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    if(iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check adjoint M2M by rotation
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_m2m_rotation_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || [adj M(1)2M(i)] - [adj M(1)2M(i)]", &
            & " || / || [adj M(1)2M(i)] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, nsph
        ! Check with beta=zero
        call fmm_m2m_rotation_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        ! Subtract reference value
        call fmm_m2m_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -one, coef(:, 1), one, coef2(:, i))
        v(i) = dnrm2((p+1)**2, coef2(:, i), 1) / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
        ! Get reference value
        call fmm_m2m_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        full_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        ! Subtract value to be checked
        call fmm_m2m_rotation_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -two, coef(:, 1), two, coef2(:, i))
        v(i) = dnrm2((p+1)**2, coef2(:, i), 1) / full_norm
        ok = v(i) .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v(i)
        end if
        if (.not. ok) stop 1
    end do
    if(iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
end subroutine check_m2m

! Check L2L for spherical harmonics
! List of explicitly checked functions:
!       fmm_l2l_ztranslate
!       fmm_l2l_ztranslate_adj
!       fmm_l2l_ztranslate_get_mat
!       fmm_l2l_ztranslate_use_mat
!       fmm_l2l_ztranslate_use_mat_adj
!       fmm_l2l_reflection
!       fmm_l2l_reflection2
!       fmm_l2l_rotation
!       fmm_l2l_reflection_adj
!       fmm_l2l_reflection2_adj
!       fmm_l2l_rotation_adj
!       fmm_l2l_reflection_get_mat
!       fmm_l2l_reflection_use_mat
!       fmm_l2l_reflection_use_mat_adj
subroutine check_l2l(p, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: p, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
    integer :: nbasis, i, j
    real(dp) :: vscales((p+1)**2), vfact(2*p+1)
    real(dp) :: dst0(3), dst(3), src(3), csph(3, nsph), rsph(nsph), v0, v1
    real(dp), dimension((p+1)**2, nsph) :: coef
    real(dp) :: ztranslate_mat((p+1)*(p+2)*(p+3)/6)
    real(dp) :: transform_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(dp) :: full_mat((p+1)**2, (p+1)**2)
    real(dp) :: z, r1(3, 3), dst1(3), full_norm, diff_norm
    logical :: ok
    real(dp), external :: dnrm2
    !! Scale inputs
    dst0 = gsrc0
    dst = alpha * gsrc
    src = alpha * gdst(:, 1)
    csph = alpha * gcsph
    rsph = abs(alpha) * grsph
    !! Preliminaries
    nbasis = (p+1) * (p+1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(p, vscales)
    vfact(1) = one
    do i = 2, 2*p+1
        vfact(i) = vfact(i-1) * sqrt(dble(i-1))
    end do
    !! Get local coefficients of the main sphere (source of L2L) and
    !! corresponding potential
    call fmm_p2l(src-csph(:, 1), one, rsph(1), p, vscales, zero, coef(:, 1))
    call fmm_l2p(dst-csph(:, 1), rsph(1), p, vscales, one, coef(:, 1), zero, &
        & v0)
    !! Check L2L OZ translation by fmm_l2l_ztranslate. Spheres 2 and 3 are
    !! explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_l2l_ztranslate"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || L(1)2L(i) + L(i)2P - L(1)2P || / ", &
            & "|| L(1)2P ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        call fmm_l2l_ztranslate(-csph(3, i), rsph(1), rsph(i), p, vscales, &
            & vfact, one, coef(:, 1), zero, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales, one, coef(:, i), &
            & zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
        call fmm_l2l_ztranslate(-csph(3, i), rsph(1), rsph(i), p, vscales, &
            & vfact, -one, coef(:, 1), two, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales, one, coef(:, i), &
            & zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
    end do
    ! Check L2L with zero OZ translation. Sphere 6 is explicitly set for this
    call fmm_l2l_ztranslate(zero, rsph(1), rsph(6), p, vscales, vfact, one, &
        & coef(:, 1), zero, coef(:, 6))
    call fmm_l2p(dst-csph(:, 6), rsph(6), p, vscales, one, coef(:, 6), zero, &
        & v1)
    v1 = abs(v1/v0 - one)
    ok = v1 .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
    end if
    if (.not. ok) stop 1
    call fmm_l2l_ztranslate(zero, rsph(1), rsph(6), p, vscales, vfact, -one, &
        & coef(:, 1), two, coef(:, 6))
    call fmm_l2p(dst-csph(:, 6), rsph(6), p, vscales, one, coef(:, 6), zero, &
        & v1)
    v1 = abs(v1/v0 - one)
    ok = v1 .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check L2L OZ translation by precomputed matrices. Spheres 2 and 3 are
    !! explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_l2l_ztranslate_get_mat"
        write(*, "(A)") "Check fmm_l2l_ztranslate_use_mat"
        write(*, "(A)") "Check fmm_l2l_scale"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || L(1)2L(i) + L(i)2P - L(1)2P || / ", &
            & "|| L(1)2P ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        call fmm_l2l_ztranslate_get_mat(-csph(3, i), rsph(1), rsph(i), p, &
            & vscales, vfact, ztranslate_mat)
        call fmm_l2l_ztranslate_use_mat(p, ztranslate_mat, one, coef(:, 1), &
            & zero, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales, one, coef(:, i), &
            & zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        call fmm_l2l_ztranslate_use_mat(p, ztranslate_mat, -one, coef(:, 1), &
            & two, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales, one, coef(:, i), &
            & zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    ! L2L with 0 shift is simply scaling
    call fmm_l2l_scale(rsph(1), rsph(6), p, one, coef(:, 1), zero, coef(:, 6))
    call fmm_l2p(dst-csph(:, 6), rsph(6), p, vscales, one, coef(:, 6), &
        & zero, v1)
    v1 = abs(v1/v0 - one)
    ok = v1 .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
    end if
    if (.not. ok) stop 1
    call fmm_l2l_scale(rsph(1), rsph(6), p, -one, coef(:, 1), two, &
        & coef(:, 6))
    call fmm_l2p(dst-csph(:, 6), rsph(6), p, vscales, one, coef(:, 6), &
        & zero, v1)
    v1 = abs(v1/v0 - one)
    ok = v1 .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check adjoint L2L OZ translation by precomputed matrices of direct L2L
    !! OZ translation. Spheres 2 and 3 are explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_l2l_ztranslate_use_mat_adj"
        write(*, "(A)") "Check fmm_l2l_scale_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || L(1)2L(i) - [adj L(1)2L(i)]^T || ", &
            & "/ || L(1)2L(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        ! Generate entire matrix of a direct L2L OZ translation
        call fmm_l2l_ztranslate_get_mat(-csph(3, i), rsph(1), rsph(i), p, &
            & vscales, vfact, ztranslate_mat)
        do j = 1, (p+1)**2
            coef(:, i) = zero
            coef(j, i) = one
            call fmm_l2l_ztranslate_use_mat(p, ztranslate_mat, one, &
                & coef(:, i), zero, full_mat(:, j))
        end do
        full_norm = dnrm2((p+1)**4, full_mat, 1)
        ! Subtract transpose of an adjoint L2L OZ translation matrix
        do j = 1, (p+1)**2
            coef(:, i) = zero
            coef(j, i) = one
            call fmm_l2l_ztranslate_use_mat_adj(p, ztranslate_mat, -two, &
                & coef(:, i), two, full_mat(j, :))
        end do
        diff_norm = dnrm2((p+1)**4, full_mat, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        call fmm_l2l_ztranslate_use_mat_adj(p, ztranslate_mat, one, &
            & coef(:, 1), zero, full_mat(:, 1))
        full_norm = dnrm2((p+1)**2, full_mat(:, 1), 1)
        call fmm_l2l_ztranslate_use_mat_adj(p, ztranslate_mat, -two, &
            & coef(:, 1), two, full_mat(:, 1))
        diff_norm = dnrm2((p+1)**2, full_mat(:, 1), 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    ! Check fmm_l2l_scale_adj. Sphere 6 is intended for this. Result is the
    ! same as with direct fmm_l2l_scale.
    call fmm_l2l_scale_adj(rsph(6), rsph(1), p, one, coef(:, 1), zero, &
        & coef(:, 6))
    full_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    call fmm_l2l_scale(rsph(1), rsph(6), p, -one, coef(:, 1), one, coef(:, 6))
    diff_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    v1 = diff_norm / full_norm
    ok = v1 .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
    end if
    call fmm_l2l_scale(rsph(1), rsph(6), p, one, coef(:, 1), zero, coef(:, 6))
    full_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    call fmm_l2l_scale_adj(rsph(6), rsph(1), p, one, coef(:, 1), -one, &
        & coef(:, 6))
    diff_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    v1 = diff_norm / full_norm
    ok = v1 .le. threshold
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check adjoint L2L OZ translation without precomputed matrices. Spheres 2
    !! and 3 are explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_l2l_ztranslate_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || [adj L(1)2L(i)] - [adj L(1)2L(i)]", &
            & " || / || [adj L(1)2L(i)] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 2, 3
        call fmm_l2l_ztranslate_get_mat(-csph(3, i), rsph(1), rsph(i), p, &
            & vscales, vfact, ztranslate_mat)
        call fmm_l2l_ztranslate_adj(csph(3, i), rsph(i), rsph(1), p, &
            & vscales, vfact, one, coef(:, 1), zero, coef(:, i))
        full_norm = dnrm2((p+1)**2, coef(:, i), 1)
        call fmm_l2l_ztranslate_use_mat_adj(p, ztranslate_mat, -two, &
            & coef(:, 1), two, coef(:, i))
        diff_norm = dnrm2((p+1)**2, coef(:, i), 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        call fmm_l2l_ztranslate_use_mat_adj(p, ztranslate_mat, one, &
            & coef(:, 1), zero, coef(:, i))
        full_norm = dnrm2((p+1)**2, coef(:, i), 1)
        call fmm_l2l_ztranslate_adj(csph(3, i), rsph(i), rsph(1), p, &
            & vscales, vfact, -two, coef(:, 1), two, coef(:, i))
        diff_norm = dnrm2((p+1)**2, coef(:, i), 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    !! Check that attempt to use L2L matrix with 0 shift raises NaNs
    call fmm_l2l_ztranslate_get_mat(zero, rsph(1), rsph(6), p, vscales, &
        & vfact, ztranslate_mat)
    call fmm_l2l_ztranslate_use_mat(p, ztranslate_mat, one, coef(:, 1), zero, &
        & coef(:, 6))
    ! Adjoint with zero translation is the same as direct with zero translation
    call fmm_l2l_ztranslate_adj(zero, rsph(6), rsph(1), p, vscales, vfact, &
        & one, coef(:, 1), zero, coef(:, 6))
    full_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    call fmm_l2l_ztranslate(zero, rsph(1), rsph(6), p, vscales, vfact, &
        & one, coef(:, 1), -one, coef(:, 6))
    diff_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    v1 = diff_norm / full_norm
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
    end if
    if (.not. ok) stop 1
    call fmm_l2l_ztranslate(zero, rsph(1), rsph(6), p, vscales, vfact, &
        & one, coef(:, 1), zero, coef(:, 6))
    full_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    call fmm_l2l_ztranslate_adj(zero, rsph(6), rsph(1), p, vscales, vfact, &
        & one, coef(:, 1), -one, coef(:, 6))
    diff_norm = dnrm2((p+1)**2, coef(:, 6), 1)
    v1 = diff_norm / full_norm
    if (iprint .gt. 0) then
        write(*, "(I3,A,L3,A,ES23.16E3)") 6, " |", ok, " | ", v1
        write(*, "(A,/)") repeat("=", 40)
    end if
    if (.not. ok) stop 1
    !! Check L2L by reflection
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_l2l_reflection"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || L(1)2L(i) + L(i)2P - L(1)2P || / ", &
            & "|| L(1)2P ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    ! ignore i=1 since coef(:, 1) would be overwritten while read -- this is
    ! an undefined behaviour
    do i = 2, nsph
        call fmm_l2l_reflection(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, one, coef(:, 1), zero, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales, one, coef(:, i), &
            & zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        call fmm_l2l_reflection(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, -one, coef(:, 1), two, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales, one, coef(:, i), &
            & zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check L2L by reflection with help of matrices
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_l2l_reflection_get_mat"
        write(*, "(A,I0)") "Check fmm_l2l_reflection_use_mat"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || L(1)2L(i) + L(i)2P - L(1)2P || / ", &
            & "|| L(1)2P ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    ! ignore i=1 since coef(:, 1) would be overwritten while read -- this is
    ! an undefined behaviour
    do i = 2, nsph
        call fmm_l2l_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, vfact, transform_mat, ztranslate_mat)
        call fmm_l2l_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, transform_mat, ztranslate_mat, one, coef(:, 1), &
            & zero, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales, one, coef(:, i), &
            & zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        call fmm_l2l_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, transform_mat, ztranslate_mat, -one, coef(:, 1), &
            & two, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales, one, coef(:, i), &
            & zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check adjoint L2L by reflection with help of matrices
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_l2l_reflection_use_mat_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || L(1)2L(i) - [adj L(1)2L(i)]^T || ", &
            & "/ || L(1)2L(i) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    ! ignore i=1 since coef(:, 1) would be overwritten while read -- this is
    ! an undefined behaviour
    do i = 2, nsph
        ! Generate entire matrix of a direct L2L translation
        call fmm_l2l_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, vfact, transform_mat, ztranslate_mat)
        do j = 1, (p+1)**2
            coef(:, i) = zero
            coef(j, i) = one
            call fmm_l2l_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
                & rsph(i), p, transform_mat, ztranslate_mat, one, coef(:, i), &
                & zero, full_mat(:, j))
        end do
        full_norm = dnrm2((p+1)**4, full_mat, 1)
        ! Subtract transpose of an adjoint L2L translation matrix
        do j = 1, (p+1)**2
            coef(:, i) = zero
            coef(j, i) = one
            call fmm_l2l_reflection_use_mat_adj(csph(:, i)-csph(:, 1), &
                & rsph(i), rsph(1), p, transform_mat, ztranslate_mat, -two, &
                & coef(:, i), two, full_mat(j, :))
        end do
        diff_norm = dnrm2((p+1)**4, full_mat, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        call fmm_l2l_reflection_use_mat_adj(csph(:, i)-csph(:, 1), rsph(i), &
            & rsph(1), p, transform_mat, ztranslate_mat, one, coef(:, 1), &
            & zero, full_mat(:, 1))
        full_norm = dnrm2((p+1)**2, full_mat(:, 1), 1)
        call fmm_l2l_reflection_use_mat_adj(csph(:, i)-csph(:, 1), rsph(i), &
            & rsph(1), p, transform_mat, ztranslate_mat, -two, coef(:, 1), &
            & two, full_mat(:, 1))
        diff_norm = dnrm2((p+1)**2, full_mat(:, 1), 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check adjoint L2L by reflection
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_l2l_reflection_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || [adj L(1)2L(i)] - [adj L(1)2L(i)]", &
            & " || / || [adj L(1)2L(i)] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    ! ignore i=1 since coef(:, 1) would be overwritten while read -- this is
    ! an undefined behaviour
    do i = 2, nsph
        ! Check with help of fmm_l2l_reflection_use_mat that is checked already
        call fmm_l2l_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, vfact, transform_mat, ztranslate_mat)
        ! Get reference value
        call fmm_l2l_reflection_use_mat_adj(csph(:, i)-csph(:, 1), &
            & rsph(i), rsph(1), p, transform_mat, ztranslate_mat, one, &
            & coef(:, 1), zero, coef(:, i))
        full_norm = dnrm2((p+1)**2, coef(:, i), 1)
        ! Subtract value to be checked
        call fmm_l2l_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -two, coef(:, 1), two, coef(:, i))
        v1 = dnrm2((p+1)**2, coef(:, i), 1) / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        ! Check with beta=zero
        call fmm_l2l_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef(:, i))
        ! Subtract reference value
        call fmm_l2l_reflection_use_mat_adj(csph(:, i)-csph(:, 1), &
            & rsph(i), rsph(1), p, transform_mat, ztranslate_mat, -one, &
            & coef(:, 1), one, coef(:, i))
        v1 = dnrm2((p+1)**2, coef(:, i), 1) / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check L2L by reflection (another implementation)
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_l2l_reflection2"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || L(1)2L(i) + L(i)2P - L(1)2P || / ", &
            & "|| L(1)2P ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    ! ignore i=1 since coef(:, 1) would be overwritten while read -- this is
    ! an undefined behaviour
    do i = 2, nsph
        call fmm_l2l_reflection2(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, one, coef(:, 1), zero, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales, one, coef(:, i), &
            & zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        call fmm_l2l_reflection2(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, -one, coef(:, 1), two, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales, one, coef(:, i), &
            & zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check adjoint L2L by reflection (another implementation)
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_l2l_reflection2_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || [adj L(1)2L(i)] - [adj L(1)2L(i)]", &
            & " || / || [adj L(1)2L(i)] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    ! ignore i=1 since coef(:, 1) would be overwritten while read -- this is
    ! an undefined behaviour
    do i = 2, nsph
        ! Get reference value
        call fmm_l2l_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef(:, i))
        full_norm = dnrm2((p+1)**2, coef(:, i), 1)
        ! Subtract value to be checked
        call fmm_l2l_reflection2_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -two, coef(:, 1), two, coef(:, i))
        v1 = dnrm2((p+1)**2, coef(:, i), 1) / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        ! Check with beta=zero
        call fmm_l2l_reflection2_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef(:, i))
        ! Subtract reference value
        call fmm_l2l_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -one, coef(:, 1), one, coef(:, i))
        v1 = dnrm2((p+1)**2, coef(:, i), 1) / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check L2L by rotation
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_l2l_rotation"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || L(1)2L(i) + L(i)2P - L(1)2P || / ", &
            & "|| L(1)2P ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    ! ignore i=1 since coef(:, 1) would be overwritten while read -- this is
    ! an undefined behaviour
    do i = 2, nsph
        call fmm_l2l_rotation(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, one, coef(:, 1), zero, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales, one, coef(:, i), &
            & zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        call fmm_l2l_rotation(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, -one, coef(:, 1), two, coef(:, i))
        call fmm_l2p(dst-csph(:, i), rsph(i), p, vscales, one, coef(:, i), &
            & zero, v1)
        v1 = abs(v1/v0 - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check adjoint L2L by rotation
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A,I0)") "Check fmm_l2l_rotation_adj"
        write(*, "(4x,A,I0)") "p = ", p
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || [adj L(1)2L(i)] - [adj L(1)2L(i)]", &
            & " || / || [adj L(1)2L(i)] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    ! ignore i=1 since coef(:, 1) would be overwritten while read -- this is
    ! an undefined behaviour
    do i = 2, nsph
        ! Get reference value
        call fmm_l2l_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef(:, i))
        full_norm = dnrm2((p+1)**2, coef(:, i), 1)
        ! Subtract value to be checked
        call fmm_l2l_rotation_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -two, coef(:, 1), two, coef(:, i))
        v1 = dnrm2((p+1)**2, coef(:, i), 1) / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        ! Check with beta=zero
        call fmm_l2l_rotation_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef(:, i))
        ! Subtract reference value
        call fmm_l2l_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -one, coef(:, 1), one, coef(:, i))
        v1 = dnrm2((p+1)**2, coef(:, i), 1) / full_norm
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
end subroutine check_l2l

subroutine check_m2l(pm, pl, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: pm, pl, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
    integer :: nbasism, nbasisl, i, j, k
    real(dp) :: vscales((pm+pl+1)**2), vfact(2*(pm+pl)+1), dst_csph(3, 2), &
        & dst_rsph(2)
    real(dp) :: src0(3), src(3), dst(3, 2), csph(3, nsph), rsph(nsph), v0(2), &
        & v1
    real(dp) :: coefm((pm+1)**2, nsph), coef2m((pm+1)**2), coefl((pl+1)**2)
    real(dp) :: ztranslate_mat((min(pm,pl)+1)*(min(pm,pl)+2)* &
        & (3*max(pm,pl)+3-min(pm,pl))/6)
    real(dp) :: transform_mat((max(pm,pl)+1)*(2*max(pm,pl)+1)* &
        & (2*max(pm,pl)+3)/3)
    real(dp) :: full_mat((pl+1)**2, (pm+1)**2)
    real(dp) :: z, r1(3, 3), dst1(3), full_norm, diff_norm
    logical :: ok
    real(dp), external :: dnrm2
    !! Scale inputs
    src0 = gsrc0
    src = alpha * gsrc
    dst = alpha * gdst
    csph = alpha * gcsph
    rsph = abs(alpha) * grsph
    dst_csph = alpha * gdst_csph
    dst_rsph = abs(alpha) * gdst_rsph
    !! Preliminaries
    nbasism = (pm+1)**2
    nbasisl = (pl+1)**2
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(pm+pl, vscales)
    vfact(1) = one
    do i = 2, 2*(pm+pl)+1
        vfact(i) = vfact(i-1) * sqrt(dble(i-1))
    end do
    !! Get multipole coefficients of all spheres by P2M
    do i = 1, nsph
        call fmm_p2m(src-csph(:, i), one, rsph(i), pm, vscales, zero, &
            & coefm(:, i))
    end do
    !! Get reference value of potentials
    v0(1) = one / dnrm2(3, src-dst(:, 1), 1)
    v0(2) = one / dnrm2(3, src-dst(:, 2), 1)
    !! Check M2L OZ translation by fmm_m2l_ztranslate. Spheres 1, 2 and 3 are
    !! explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2l_ztranslate"
        write(*, "(4x,A,I0)") "pm = ", pm
        write(*, "(4x,A,I0)") "pl = ", pl
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || P2M(i) + M(i)2L(1) + L(1)2P - P2P", &
            & " || / || P2P ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, 3
        call fmm_m2l_ztranslate(csph(3, i)-dst_csph(3, 1), rsph(i), &
            & dst_rsph(1), pm, pl, vscales, vfact, one, coefm(:, i), zero, &
            & coefl)
        call fmm_l2p(dst(:, 1)-dst_csph(:, 1), dst_rsph(1), pl, vscales, one, &
            & coefl, zero, v1)
        v1 = abs(v1/v0(1) - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
        call fmm_m2l_ztranslate(csph(3, i)-dst_csph(3, 1), rsph(i), &
            & dst_rsph(1), pm, pl, vscales, vfact, -one, coefm(:, i), two, &
            & coefl)
        call fmm_l2p(dst(:, 1)-dst_csph(:, 1), dst_rsph(1), pl, vscales, one, &
            & coefl, zero, v1)
        v1 = abs(v1/v0(1) - one)
        ok = v1 .le. threshold
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check M2L OZ translation by precomputed matrices. Spheres 1, 2 and 3 are
    !! explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2l_ztranslate_get_mat"
        write(*, "(A)") "Check fmm_m2l_ztranslate_use_mat"
        write(*, "(4x,A,I0)") "pm = ", pm
        write(*, "(4x,A,I0)") "pl = ", pl
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || M(i)2L(1) - M(i)2L(1) || / ", &
            & "|| M(i)2L(1) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, 3
        call fmm_m2l_ztranslate_get_mat(csph(3, i)-dst_csph(3, 1), rsph(i), &
            & dst_rsph(1), pm, pl, vscales, vfact, ztranslate_mat)
        call fmm_m2l_ztranslate_use_mat(pm, pl, ztranslate_mat, one, &
            & coefm(:, i), zero, coefl)
        full_norm = dnrm2((pl+1)**2, coefl, 1)
        call fmm_m2l_ztranslate(csph(3, i)-dst_csph(3, 1), rsph(i), &
            & dst_rsph(1), pm, pl, vscales, vfact, -one, coefm(:, i), one, &
            & coefl)
        diff_norm = dnrm2((pl+1)**2, coefl, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. 10d0*epsilon(one)
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
        call fmm_m2l_ztranslate(csph(3, i)-dst_csph(3, 1), rsph(i), &
            & dst_rsph(1), pm, pl, vscales, vfact, one, coefm(:, i), zero, &
            & coefl)
        full_norm = dnrm2((pl+1)**2, coefl, 1)
        call fmm_m2l_ztranslate_use_mat(pm, pl, ztranslate_mat, -two, &
            & coefm(:, i), two, coefl)
        diff_norm = dnrm2((pl+1)**2, coefl, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. 10d0*epsilon(one)
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Since centers of M and L harmonics are not allowed to be in the same
    !! place, there is no fmm_m2l_scale function.
    !! Check adjoint M2L OZ translation by precomputed matrices of direct L2L
    !! OZ translation. Spheres 1, 2 and 3 are explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2l_ztranslate_use_mat_adj"
        write(*, "(4x,A,I0)") "pm = ", pm
        write(*, "(4x,A,I0)") "pl = ", pl
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || M(i)2L(1) - [adj M(i)2L(1)]^T || ", &
            & "/ || M(i)2L(1) ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, 3
        ! Generate entire matrix of a direct M2L OZ translation
        call fmm_m2l_ztranslate_get_mat(csph(3, i)-dst_csph(3, 1), rsph(i), &
            & dst_rsph(1), pm, pl, vscales, vfact, ztranslate_mat)
        do j = 1, (pm+1)**2
            coef2m(:) = zero
            coef2m(j) = one
            call fmm_m2l_ztranslate_use_mat(pm, pl, ztranslate_mat, one, &
                & coef2m, zero, coefl)
            full_mat(:, j) = coefl
        end do
        full_norm = dnrm2(((pm+1)*(pl+1))**2, full_mat, 1)
        ! Subtract transpose of an adjoint M2L OZ translation matrix
        do j = 1, (pl+1)**2
            coefl(:) = zero
            coefl(j) = one
            call fmm_m2l_ztranslate_use_mat_adj(pl, pm, ztranslate_mat, -one, &
                & coefl, one, full_mat(j, :))
        end do
        diff_norm = dnrm2(((pm+1)*(pl+1))**2, full_mat, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. 10d0*epsilon(one)
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
        call fmm_m2l_ztranslate_use_mat(pm, pl, ztranslate_mat, one, &
            & coefm(:, i), zero, coefl)
        call fmm_m2l_ztranslate_use_mat_adj(pl, pm, ztranslate_mat, one, &
            & coefl, zero, coef2m)
        full_norm = dnrm2((pm+1)**2, coef2m, 1)
        call fmm_m2l_ztranslate_use_mat_adj(pl, pm, ztranslate_mat, -two, &
            & coefl, two, coef2m)
        diff_norm = dnrm2((pm+1)**2, coef2m, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. 10d0*epsilon(one)
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check adjoint M2L OZ translation without precomputed matrices. Spheres
    !! 1, 2 and 3 are explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, "(/,A)") repeat("=", 40)
        write(*, "(A)") "Check fmm_m2l_ztranslate_adj"
        write(*, "(4x,A,I0)") "pm = ", pm
        write(*, "(4x,A,I0)") "pl = ", pl
        write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
        write(*, "(4x,A,A)") "err(i) = || [adj M(i)2L(1)] - [adj M(i)2L(1)]", &
            & " || / || [adj M(i)2L(1)] ||"
        write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
        write(*, "(A)") repeat("=", 40)
        write(*, "(A)") "  i | ok | err(i)"
        write(*, "(A)") repeat("=", 40)
    end if
    do i = 1, 3
        call fmm_m2l_ztranslate_get_mat(csph(3, i)-dst_csph(3, 1), rsph(i), &
            & dst_rsph(1), pm, pl, vscales, vfact, ztranslate_mat)
        call fmm_m2l_ztranslate_use_mat(pm, pl, ztranslate_mat, one, &
            & coefm(:, i), zero, coefl)
        call fmm_m2l_ztranslate_adj(dst_csph(3, 1)-csph(3, i), dst_rsph(1), &
            & rsph(i), pl, pm, vscales, vfact, one, coefl, zero, coef2m)
        full_norm = dnrm2((pm+1)**2, coef2m, 1)
        call fmm_m2l_ztranslate_use_mat_adj(pl, pm, ztranslate_mat, -two, &
            & coefl, two, coef2m)
        diff_norm = dnrm2((pm+1)**2, coef2m, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. 10d0*epsilon(one)
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
        call fmm_m2l_ztranslate_use_mat_adj(pl, pm, ztranslate_mat, one, &
            & coefl, zero, coef2m)
        full_norm = dnrm2((pm+1)**2, coef2m, 1)
        call fmm_m2l_ztranslate_adj(dst_csph(3, 1)-csph(3, i), dst_rsph(1), &
            & rsph(i), pl, pm, vscales, vfact, -two, coefl, two, coef2m)
        diff_norm = dnrm2((pm+1)**2, coef2m, 1)
        v1 = diff_norm / full_norm
        ok = v1 .le. 10d0*epsilon(one)
        if (iprint .gt. 0) then
            write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
        end if
        if (.not. ok) stop 1
    end do
    if (iprint .gt. 0) then
        write(*, "(A,/)") repeat("=", 40)
    end if
    !! Check M2L by reflection
    do k = 1, 2
        if (iprint .gt. 0) then
            write(*, "(/,A)") repeat("=", 40)
            write(*, "(A)") "Check fmm_m2l_reflection"
            write(*, "(4x,A,I0)") "pm = ", pm
            write(*, "(4x,A,I0)") "pl = ", pl
            write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
            write(*, "(4x,A,I0,A,I0,A)") "err(i) = || P2M(i) + M(i)2L(", k, &
                & ") + L(", k, ")2P - P2P || / || P2P ||"
            write(*, "(4x,A,ES23.16E3)") "threshold = ", threshold
            write(*, "(A)") repeat("=", 40)
            write(*, "(A)") "  i | ok | err(i)"
            write(*, "(A)") repeat("=", 40)
        end if
        do i = 1, nsph
            call fmm_m2l_reflection(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, vfact, one, coefm(:, i), &
                & zero, coefl)
            call fmm_l2p(dst(:, k)-dst_csph(:, k), dst_rsph(k), pl, vscales, &
                & one, coefl, zero, v1)
            v1 = abs(v1/v0(k) - one)
            ok = v1 .le. threshold
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
            call fmm_m2l_reflection(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, vfact, -one, coefm(:, i), &
                & two, coefl)
            call fmm_l2p(dst(:, k)-dst_csph(:, k), dst_rsph(k), pl, vscales, &
                & one, coefl, zero, v1)
            v1 = abs(v1/v0(k) - one)
            ok = v1 .le. threshold
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
        end do
        if (iprint .gt. 0) then
            write(*, "(A,/)") repeat("=", 40)
        end if
    end do
    !! Check M2L by reflection with help of matrices
    do k = 1, 2
        if (iprint .gt. 0) then
            write(*, "(/,A)") repeat("=", 40)
            write(*, "(A)") "Check fmm_m2l_reflection_get_mat"
            write(*, "(A)") "Check fmm_m2l_reflection_use_mat"
            write(*, "(4x,A,I0)") "pm = ", pm
            write(*, "(4x,A,I0)") "pl = ", pl
            write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
            write(*, "(4x,3(A,I0),A)") "err(i) = || M(i)2L(", k, &
                & ") - M(i)2L(", k, ") || / || M(i)2L(", k, ") ||"
            write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
            write(*, "(A)") repeat("=", 40)
            write(*, "(A)") "  i | ok | err(i)"
            write(*, "(A)") repeat("=", 40)
        end if
        do i = 1, nsph
            call fmm_m2l_reflection_get_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(k), pm, pl, vscales, vfact, &
                & transform_mat, ztranslate_mat)
            call fmm_m2l_reflection_use_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(k), pm, pl, transform_mat, &
                & ztranslate_mat, one, coefm(:, i), zero, coefl)
            full_norm = dnrm2((pl+1)**2, coefl, 1)
            call fmm_m2l_reflection(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, vfact, -one, coefm(:, i), &
                & one, coefl)
            diff_norm = dnrm2((pl+1)**2, coefl, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
            call fmm_m2l_reflection(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, vfact, one, coefm(:, i), &
                & zero, coefl)
            full_norm = dnrm2((pl+1)**2, coefl, 1)
            call fmm_m2l_reflection_use_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(K), pm, pl, transform_mat, &
                & ztranslate_mat, -two, coefm(:, i), two, coefl)
            diff_norm = dnrm2((pl+1)**2, coefl, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
        end do
        if (iprint .gt. 0) then
            write(*, "(A,/)") repeat("=", 40)
        end if
    end do
    !! Check adjoint M2L by reflection with help of matrices
    do k = 1, 2
        if (iprint .gt. 0) then
            write(*, "(/,A)") repeat("=", 40)
            write(*, "(A)") "Check fmm_m2l_reflection_use_mat_adj"
            write(*, "(4x,A,I0)") "pm = ", pm
            write(*, "(4x,A,I0)") "pl = ", pl
            write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
            write(*, "(4x,3(A,I0),A)") "err(i) = || M(i)2L(", k, &
                & ") - [adj M(i)2L(", k, ")]^T || / || M(i)2L(", k, ") ||"
            write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
            write(*, "(A)") repeat("=", 40)
            write(*, "(A)") "  i | ok | err(i)"
            write(*, "(A)") repeat("=", 40)
        end if
        do i = 1, nsph
            ! Generate entire matrix of a direct M2L translation
            call fmm_m2l_reflection_get_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(k), pm, pl, vscales, vfact, &
                & transform_mat, ztranslate_mat)
            do j = 1, (pm+1)**2
                coef2m = zero
                coef2m(j) = one
                call fmm_m2l_reflection_use_mat(csph(:, i)-dst_csph(:, k), &
                    & rsph(i), dst_rsph(k), pm, pl, transform_mat, &
                    & ztranslate_mat, one, coef2m, zero, full_mat(:, j))
            end do
            full_norm = dnrm2(((pm+1)*(pl+1))**2, full_mat, 1)
            ! Subtract transpose of an adjoint L2L translation matrix
            do j = 1, (pl+1)**2
                coefl = zero
                coefl(j) = one
                call fmm_m2l_reflection_use_mat_adj( &
                    & dst_csph(:, k)-csph(:, i), dst_rsph(k), rsph(i), pl, &
                    & pm, transform_mat, ztranslate_mat, -two, coefl, two, &
                    & full_mat(j, :))
            end do
            diff_norm = dnrm2(((pm+1)*(pl+1))**2, full_mat, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
            call fmm_m2l_reflection_use_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(K), pm, pl, transform_mat, &
                & ztranslate_mat, one, coefm(:, i), zero, coefl)
            call fmm_m2l_reflection_use_mat_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, transform_mat, &
                & ztranslate_mat, one, coefl, zero, coef2m)
            full_norm = dnrm2((pm+1)**2, coef2m, 1)
            call fmm_m2l_reflection_use_mat_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, transform_mat, &
                & ztranslate_mat, -two, coefl, two, coef2m)
            diff_norm = dnrm2((pm+1)**2, coef2m, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
        end do
        if (iprint .gt. 0) then
            write(*, "(A,/)") repeat("=", 40)
        end if
    end do
    !! Check adjoint M2L by reflection
    do k = 1, 2
        if (iprint .gt. 0) then
            write(*, "(/,A)") repeat("=", 40)
            write(*, "(A)") "Check fmm_m2l_reflection_adj"
            write(*, "(4x,A,I0)") "pm = ", pm
            write(*, "(4x,A,I0)") "pl = ", pl
            write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
            write(*, "(4x,3(A,I0),A)") "err(i) = || [adj M(i)2L(", k, &
                & ")] - [adj M(i)2L(", k, ")] || / || [adj M(i)2L(", k, ")] ||"
            write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
            write(*, "(A)") repeat("=", 40)
            write(*, "(A)") "  i | ok | err(i)"
            write(*, "(A)") repeat("=", 40)
        end if
        do i = 1, nsph
            ! Check with help of fmm_m2l_reflection_use_mat_adj
            call fmm_m2l_reflection_get_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(k), pm, pl, vscales, vfact, &
                & transform_mat, ztranslate_mat)
            call fmm_m2l_reflection_use_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(k), pm, pl, transform_mat, &
                & ztranslate_mat, one, coefm(:, i), zero, coefl)
            ! Get reference value
            call fmm_m2l_reflection_use_mat_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, transform_mat, &
                & ztranslate_mat, one, coefl, zero, coef2m)
            full_norm = dnrm2((pm+1)**2, coef2m, 1)
            ! Subtract value to be checked
            call fmm_m2l_reflection_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, vscales, vfact, -two, coefl, &
                & two, coef2m)
            diff_norm = dnrm2((pm+1)**2, coef2m, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
            ! Check with beta=zero
            call fmm_m2l_reflection_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, vscales, vfact, one, coefl, &
                & zero, coef2m)
            full_norm = dnrm2((pm+1)**2, coef2m, 1)
            ! Subtract reference value
            call fmm_m2l_reflection_use_mat_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, transform_mat, &
                & ztranslate_mat, -one, coefl, one, coef2m)
            diff_norm = dnrm2((pm+1)**2, coef2m, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
        end do
        if (iprint .gt. 0) then
            write(*, "(A,/)") repeat("=", 40)
        end if
    end do
    !! Check M2L by reflection2 (another implementation)
    do k = 1, 2
        if (iprint .gt. 0) then
            write(*, "(/,A)") repeat("=", 40)
            write(*, "(A)") "Check fmm_m2l_reflection2"
            write(*, "(4x,A,I0)") "pm = ", pm
            write(*, "(4x,A,I0)") "pl = ", pl
            write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
            write(*, "(4x,3(A,I0),A)") "err(i) = || M(i)2L(", k, &
                & ") - M(i)2L(", k, ") || / || M(i)2L(", k, ") ||"
            write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
            write(*, "(A)") repeat("=", 40)
            write(*, "(A)") "  i | ok | err(i)"
            write(*, "(A)") repeat("=", 40)
        end if
        do i = 1, nsph
            call fmm_m2l_reflection2(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, vfact, one, coefm(:, i), &
                & zero, coefl)
            full_norm = dnrm2((pl+1)**2, coefl, 1)
            call fmm_m2l_reflection(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, vfact, -one, coefm(:, i), &
                & one, coefl)
            diff_norm = dnrm2((pl+1)**2, coefl, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
            call fmm_m2l_reflection(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, vfact, one, coefm(:, i), &
                & zero, coefl)
            full_norm = dnrm2((pl+1)**2, coefl, 1)
            call fmm_m2l_reflection2(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, vfact, -two, coefm(:, i), &
                & two, coefl)
            diff_norm = dnrm2((pl+1)**2, coefl, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
        end do
        if (iprint .gt. 0) then
            write(*, "(A,/)") repeat("=", 40)
        end if
    end do
    !! Check adjoint M2L by reflection2
    do k = 1, 2
        if (iprint .gt. 0) then
            write(*, "(/,A)") repeat("=", 40)
            write(*, "(A)") "Check fmm_m2l_reflection2_adj"
            write(*, "(4x,A,I0)") "pm = ", pm
            write(*, "(4x,A,I0)") "pl = ", pl
            write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
            write(*, "(4x,3(A,I0),A)") "err(i) = || [adj M(i)2L(", k, &
                & ")] - [adj M(i)2L(", k, ")] || / || [adj M(i)2L(", k, ")] ||"
            write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
            write(*, "(A)") repeat("=", 40)
            write(*, "(A)") "  i | ok | err(i)"
            write(*, "(A)") repeat("=", 40)
        end if
        do i = 1, nsph
            ! Check with help of fmm_m2l_reflection_use_mat_adj
            call fmm_m2l_reflection_get_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(k), pm, pl, vscales, vfact, &
                & transform_mat, ztranslate_mat)
            call fmm_m2l_reflection_use_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(k), pm, pl, transform_mat, &
                & ztranslate_mat, one, coefm(:, i), zero, coefl)
            ! Get reference value
            call fmm_m2l_reflection_use_mat_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, transform_mat, &
                & ztranslate_mat, one, coefl, zero, coef2m)
            full_norm = dnrm2((pm+1)**2, coef2m, 1)
            ! Subtract value to be checked
            call fmm_m2l_reflection2_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, vscales, vfact, -two, coefl, &
                & two, coef2m)
            diff_norm = dnrm2((pm+1)**2, coef2m, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
            ! Check with beta=zero
            call fmm_m2l_reflection2_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, vscales, vfact, one, coefl, &
                & zero, coef2m)
            full_norm = dnrm2((pm+1)**2, coef2m, 1)
            ! Subtract reference value
            call fmm_m2l_reflection_use_mat_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, transform_mat, &
                & ztranslate_mat, -one, coefl, one, coef2m)
            diff_norm = dnrm2((pm+1)**2, coef2m, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
        end do
        if (iprint .gt. 0) then
            write(*, "(A,/)") repeat("=", 40)
        end if
    end do
    !! Check M2L by rotation
    do k = 1, 2
        if (iprint .gt. 0) then
            write(*, "(/,A)") repeat("=", 40)
            write(*, "(A)") "Check fmm_m2l_rotation"
            write(*, "(4x,A,I0)") "pm = ", pm
            write(*, "(4x,A,I0)") "pl = ", pl
            write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
            write(*, "(4x,3(A,I0),A)") "err(i) = || M(i)2L(", k, &
                & ") - M(i)2L(", k, ") || / || M(i)2L(", k, ") ||"
            write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
            write(*, "(A)") repeat("=", 40)
            write(*, "(A)") "  i | ok | err(i)"
            write(*, "(A)") repeat("=", 40)
        end if
        do i = 1, nsph
            call fmm_m2l_rotation(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, vfact, one, coefm(:, i), &
                & zero, coefl)
            full_norm = dnrm2((pl+1)**2, coefl, 1)
            call fmm_m2l_reflection(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, vfact, -one, coefm(:, i), &
                & one, coefl)
            diff_norm = dnrm2((pl+1)**2, coefl, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
            call fmm_m2l_reflection(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, vfact, one, coefm(:, i), &
                & zero, coefl)
            full_norm = dnrm2((pl+1)**2, coefl, 1)
            call fmm_m2l_rotation(csph(:, i)-dst_csph(:, k), rsph(i), &
                & dst_rsph(k), pm, pl, vscales, vfact, -two, coefm(:, i), &
                & two, coefl)
            diff_norm = dnrm2((pl+1)**2, coefl, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
        end do
        if (iprint .gt. 0) then
            write(*, "(A,/)") repeat("=", 40)
        end if
    end do
    !! Check adjoint M2L by rotation
    do k = 1, 2
        if (iprint .gt. 0) then
            write(*, "(/,A)") repeat("=", 40)
            write(*, "(A)") "Check fmm_m2l_rotation_adj"
            write(*, "(4x,A,I0)") "pm = ", pm
            write(*, "(4x,A,I0)") "pl = ", pl
            write(*, "(4x,A,ES23.16E3)") "alpha = ", alpha
            write(*, "(4x,3(A,I0),A)") "err(i) = || [adj M(i)2L(", k, &
                & ")] - [adj M(i)2L(", k, ")] || / || [adj M(i)2L(", k, ")] ||"
            write(*, "(4x,A,ES23.16E3)") "threshold = ", 10d0*epsilon(one)
            write(*, "(A)") repeat("=", 40)
            write(*, "(A)") "  i | ok | err(i)"
            write(*, "(A)") repeat("=", 40)
        end if
        do i = 1, nsph
            ! Check with help of fmm_m2l_reflection_use_mat_adj
            call fmm_m2l_reflection_get_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(k), pm, pl, vscales, vfact, &
                & transform_mat, ztranslate_mat)
            call fmm_m2l_reflection_use_mat(csph(:, i)-dst_csph(:, k), &
                & rsph(i), dst_rsph(k), pm, pl, transform_mat, &
                & ztranslate_mat, one, coefm(:, i), zero, coefl)
            ! Get reference value
            call fmm_m2l_reflection_use_mat_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, transform_mat, &
                & ztranslate_mat, one, coefl, zero, coef2m)
            full_norm = dnrm2((pm+1)**2, coef2m, 1)
            ! Subtract value to be checked
            call fmm_m2l_rotation_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, vscales, vfact, -two, coefl, &
                & two, coef2m)
            diff_norm = dnrm2((pm+1)**2, coef2m, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
            ! Check with beta=zero
            call fmm_m2l_rotation_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, vscales, vfact, one, coefl, &
                & zero, coef2m)
            full_norm = dnrm2((pm+1)**2, coef2m, 1)
            ! Subtract reference value
            call fmm_m2l_reflection_use_mat_adj(dst_csph(:, k)-csph(:, i), &
                & dst_rsph(k), rsph(i), pl, pm, transform_mat, &
                & ztranslate_mat, -one, coefl, one, coef2m)
            diff_norm = dnrm2((pm+1)**2, coef2m, 1)
            v1 = diff_norm / full_norm
            ok = v1 .le. 10d0*epsilon(one)
            if (iprint .gt. 0) then
                write(*, "(I3,A,L3,A,ES23.16E3)") i, " |", ok, " | ", v1
            end if
            if (.not. ok) stop 1
        end do
        if (iprint .gt. 0) then
            write(*, "(A,/)") repeat("=", 40)
        end if
    end do
end subroutine

subroutine check_tree_rib(alpha)
    real(dp), intent(in) :: alpha
    integer, parameter :: nsph = 10
    real(dp) :: csph(3, nsph), rsph(nsph), csph2(3, nsph), rsph2(nsph), &
        & cnode(3, 2*nsph-1), rnode(2*nsph-1)
    integer :: order(nsph), i, reorder(nsph), cluster(2, 2*nsph-1), &
        & children(2, 2*nsph-1), parent(2*nsph-1), snode(nsph)
    ! Scale inputs
    csph(:, 1) = alpha * (/1d0, 1d0, 1d0/)
    csph(:, 2) = alpha * (/2d0, 2d0, 2d0/)
    csph(:, 3) = alpha * (/1d0, 1d0, 3d0/)
    csph(:, 4) = alpha * (/2d0, 2d0, 4d0/)
    csph(:, 5) = alpha * (/1d0, 1d0, 5d0/)
    csph(:, 6) = alpha * (/2d0, 2d0, 6d0/)
    csph(:, 7) = alpha * (/1d0, 1d0, 7d0/)
    csph(:, 8) = alpha * (/2d0, 2d0, 8d0/)
    csph(:, 9) = alpha * (/1d0, 1d0, 9d0/)
    csph(:, 10) = alpha * (/2d0, 2d0, 10d0/)
    rsph = abs(alpha) * (/1d-1, 2d-1, 3d-1, 4d-1, 5d-1, 6d-1, 7d-1, 8d-1, &
        & 9d-1, 1d0/)
    ! Reorder inputs
    order = (/4, 2, 8, 5, 7, 1, 9, 6, 10, 3/)
    do i = 1, nsph
        csph2(:, i) = csph(:, order(i))
        rsph2(i) = rsph(order(i))
    end do
    ! Build a recursive inertial binary tree
    call tree_rib_build(nsph, csph2, rsph2, reorder, cluster, children, &
        & parent, cnode, rnode, snode)
end subroutine check_tree_rib

end program test_dd_core

