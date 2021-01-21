!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/dd_core.f90
!! Tests for dd_core module
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2020-12-17

program test_dd_core
use dd_core
implicit none

!character(len=255) :: testname, iprint_string
integer :: argc, iprint=1
integer :: p=10, i

! Input points and spheres
integer ,parameter :: nsph = 10
real(dp) :: gsrc0(3), gsrc(3), gdst(3), gcsph(3, nsph), grsph(nsph)
! Some implementation do not work for alpha=1d-307, so range of test values is
! reduced. This can be fixed by enforcing proper order of multiplications like
! a*b*c into a*(b*c).
real(dp) :: alpha(3)=(/1d0, 1d-300, 1d+300/)

real(dp), external :: dnrm2
! Set inputs
gsrc0 = zero
gsrc = (/1d-1, 2d-1, 3d-1/)
gdst = (/14d0, 12.1d0, 11.2d0/)
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
!call check_ddinit_args()

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

! Check P2M and M2P operations of the FMM
do i = 1, size(alpha)
    call check_m2m(0, alpha(i), iprint, 10d0*epsilon(zero))
    call check_m2m(1, alpha(i), iprint, 10d0*epsilon(zero))
    call check_m2m(p, alpha(i), iprint, 10d0*epsilon(zero))
end do


contains

subroutine check_ddinit_args()
    ! Example of correct args
    integer :: n=0, model=1, lmax=0, ngrid=0, force=1, fmm=1, pm=0, pl=0, &
        & iprint=0, nngmax=1
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
        & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    ! Check different correct inputs with different n <= 10 (hardcoded value)
    do i = 1, 10
        call ddinit(i, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
            & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect input n = -1
    i = -1
    call ddinit(i, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -1) stop 1
    call ddfree(dd_data)
    ! Check all possible models with other correct inputs
    do i = 1, 3
        call ddinit(n, x, y, z, rvdw, i, lmax, ngrid, force, fmm, pm, &
            & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect models
    i = -1
    call ddinit(n, x, y, z, rvdw, i, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -6) stop 1
    call ddfree(dd_data)
    i = 4
    call ddinit(n, x, y, z, rvdw, i, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -6) stop 1
    call ddfree(dd_data)
    ! Check correct lmax
    do i = 1, 6
        call ddinit(n, x, y, z, rvdw, model, i, ngrid, force, fmm, pm, &
            & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect lmax < 0
    i = -1
    call ddinit(n, x, y, z, rvdw, model, i, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -7) stop 1
    call ddfree(dd_data)
    ! Check correct ngrid
    do i = 0, 1000, 100
        j = i
        call ddinit(n, x, y, z, rvdw, model, lmax, j, force, fmm, pm, &
            & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect ngrid < 0
    i = -1
    call ddinit(n, x, y, z, rvdw, model, lmax, i, force, fmm, pm, &
        & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -8) stop 1
    call ddfree(dd_data)
    ! Check correct force
    do i = 0, 1
        call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, i, fmm, pm, &
            & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect force
    i = -1
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, i, fmm, pm, &
        & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -9) stop 1
    call ddfree(dd_data)
    i = 2
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, i, fmm, pm, &
        & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -9) stop 1
    call ddfree(dd_data)
    ! Check correct fmm
    do i = 0, 1
        call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, i, pm, &
            & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect fmm
    i = -1
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, i, pm, &
        & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -10) stop 1
    call ddfree(dd_data)
    i = 2
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, i, pm, &
        & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -10) stop 1
    call ddfree(dd_data)
    ! Check correct pm (ignored if fmm=0)
    j = 0
    do i = -2, 2
        call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, j, i, &
            & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check correct pm (fmm=1)
    j = 1
    do i = 0, 20, 5
        call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, j, i, &
            & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect pm (fmm=1)
    j = 1
    i = -1
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, j, i, &
        & pl, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -11) stop 1
    call ddfree(dd_data)
    ! Check correct pl (ignored if fmm=0)
    j = 0
    do i = -2, 2
        call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, j, pm, &
            & i, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check correct pl (fmm=1)
    j = 1
    do i = 0, 20, 5
        call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, j, pm, &
            & i, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect pl (fmm=1)
    j = 1
    i = -1
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, j, pm, &
        & i, iprint, nngmax, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -12) stop 1
    call ddfree(dd_data)
    ! Check correct iprint
    do i = 0, 10
        call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
            & pl, i, nngmax, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect iprint
    i = -1
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, i, nngmax, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -13) stop 1
    call ddfree(dd_data)
    ! Check correct nngmax
    do i = 1, 10
        call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
            & pl, iprint, i, se, eta, eps, kappa, dd_data, info)
        if (info .ne. 0) stop 1
        call ddfree(dd_data)
    end do
    ! Check incorrect nngmax
    i = 0
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, i, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -14) stop 1
    call ddfree(dd_data)
    i = -1
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, i, se, eta, eps, kappa, dd_data, info)
    if (info .ne. -14) stop 1
    call ddfree(dd_data)
    ! Check correct se
    tmp = -one
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, tmp, eta, eps, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    tmp = zero
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, tmp, eta, eps, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    tmp = one
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, tmp, eta, eps, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    ! Check incorrect se
    tmp = 1.01d0
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, tmp, eta, eps, kappa, dd_data, info)
    if (info .ne. -15) stop 1
    call ddfree(dd_data)
    tmp = -1.01d0
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, tmp, eta, eps, kappa, dd_data, info)
    if (info .ne. -15) stop 1
    call ddfree(dd_data)
    ! Check correct eta
    tmp = zero
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, se, tmp, eps, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    tmp = pt5
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, se, tmp, eps, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    tmp = one
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, se, tmp, eps, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    ! Check incorrect eta
    tmp = 1.01d0
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, se, tmp, eps, kappa, dd_data, info)
    if (info .ne. -16) stop 1
    call ddfree(dd_data)
    tmp = -1d-2
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, se, tmp, eps, kappa, dd_data, info)
    if (info .ne. -16) stop 1
    call ddfree(dd_data)
    ! Check correct eps
    tmp = zero
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, se, eta, tmp, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    tmp = pt5
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, se, eta, tmp, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    tmp = one
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, se, eta, tmp, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    tmp = dble(1000)
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, se, eta, tmp, kappa, dd_data, info)
    if (info .ne. 0) stop 1
    call ddfree(dd_data)
    ! Check incorrect eps
    tmp = -1d-2
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, se, eta, tmp, kappa, dd_data, info)
    if (info .ne. -17) stop 1
    call ddfree(dd_data)
    ! Check incorrect kappa
    tmp = -1d-2
    call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, &
        & pl, iprint, nngmax, se, eta, eps, tmp, dd_data, info)
    if (info .ne. -18) stop 1
    call ddfree(dd_data)
end subroutine check_ddinit_args

! Check P2M, M2M and M2P for spherical harmonics
! List of explicitly checked functions:
!       fmm_p2m
!       fmm_m2p
!       coord_reflect_get_mat
!       fmm_sph_transform
!       fmm_sph_transform_get_mat
!       fmm_sph_transform_use_mat
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
! List of implicitly checked functions:
!       ylmscale
!       ylmbas
!       trgev
!       polleg
! TODO explicitly:
!       fmm_sph_rotate_oz
!       fmm_sph_transform_oxz
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
    dst = alpha * gdst
    csph = alpha * gcsph
    rsph = alpha * grsph
    ! Preliminaries
    nbasis = (p+1) * (p+1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(p, vscales)
    !! Check P2M parameters q and beta
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0,A,ES11.4E3)") "Check P2M q and beta params for p=", &
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
        write(*, "(A,A,I0,A,ES11.4E3)") "Check M2P alpha and beta params ", &
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
        write(*, "(A,I0,A,ES11.4E3)") "Check P2M and M2P for p=", p, &
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
    src = alpha * gdst
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
        write(*, "(A,I0,A,ES11.4E3)") "Check P2L q and beta params for p=", &
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
        write(*, "(A,A,I0,A,ES11.4E3)") "Check L2P alpha and beta params ", &
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
        write(*, "(A,I0,A,ES11.4E3)") "Check P2L and L2P for p=", p, &
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
    dst = alpha * gdst
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
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_ztranslate for p=", p
        write(*, "(A)") "================================"
        write(*, *) "threshold=", threshold
    end if
    do i = 2, 3
        call fmm_m2m_ztranslate(-csph(3, i), rsph(1), rsph(i), p, vscales, &
            & vfact, one, coef(:, 1), zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") "|| P2M(1) + M(1)2M(", &
                & i, ") - P2M(", i, ") ||  /  || P2M(", i ,") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
        call fmm_m2m_ztranslate(-csph(3, i), rsph(1), rsph(i), p, vscales, &
            & vfact, -one, coef(:, 1), two, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") "|| P2M(1) + M(1)2M(", &
                & i, ") - P2M(", i, ") ||  /  || P2M(", i ,") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
    end do
    !! Check M2M with zero OZ translation. Sphere 6 is explicitly set for this
    call fmm_m2m_ztranslate(zero, rsph(1), rsph(6), p, vscales, vfact, one, &
        & coef(:, 1), zero, coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    if (iprint .gt. 0) then
        write(*, "(A,ES24.16E3)") &
            & "|| P2M(1) + M(1)2M(6) - P2M(6) ||  /  || P2M(6) || =", v(6)
    end if
    if (v(6) .gt. threshold) stop 1
    call fmm_m2m_ztranslate(zero, rsph(1), rsph(6), p, vscales, vfact, -one, &
        & coef(:, 1), two, coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    if (iprint .gt. 0) then
        write(*, "(A,ES24.16E3)") &
            & "|| P2M(1) + M(1)2M(6) - P2M(6) ||  /  || P2M(6) || =", v(6)
    end if
    if (v(6) .gt. threshold) stop 1
    !! Check M2M OZ translation by precomputed matrices. Spheres 2 and 3 are
    !! explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_ztranslate_get_mat for p=", p
        write(*, "(A,I0)") "Check fmm_m2m_ztranslate_use_mat for p=", p
        write(*, "(A,I0)") "Check fmm_m2m_scale for p=", p
        write(*, "(A)") "================================"
        write(*, *) "threshold=", threshold
    end if
    do i = 2, 3
        call fmm_m2m_ztranslate_get_mat(-csph(3, i), rsph(1), rsph(i), p, &
            & vscales, vfact, ztranslate_mat)
        call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, one, coef(:, 1), &
            & zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") "|| P2M(1) + M(1)2M(", &
                & i, ") - P2M(", i, ") ||  /  || P2M(", i ,") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
        call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, -one, coef(:, 1), &
            & two, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") "|| P2M(1) + M(1)2M(", &
                & i, ") - P2M(", i, ") ||  /  || P2M(", i ,") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
    end do
    ! M2M with 0 shift is simply scaling
    call fmm_m2m_scale(rsph(1), rsph(6), p, one, coef(:, 1), zero, coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    if (iprint .gt. 0) then
        write(*, "(A,ES24.16E3)") &
            & "|| P2M(1) + M(1)2M(6) - P2M(6) ||  /  || P2M(6) || =", v(i)
    end if
    call fmm_m2m_scale(rsph(1), rsph(6), p, -one, coef(:, 1), two, &
        & coef2(:, 6))
    v(6) = dnrm2(nbasis, coef2(:, 6)-coef(:, 6), 1) / &
        & dnrm2(nbasis, coef(:, 6), 1)
    if (iprint .gt. 0) then
        write(*, "(A,ES24.16E3)") &
            & "|| P2M(1) + M(1)2M(6) - P2M(6) ||  /  || P2M(6) || =", v(i)
    end if
    !! Check adjoint M2M OZ translation by precomputed matrices of direct M2M
    !! OZ translation. Spheres 2 and 3 are explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_ztranslate_use_mat_adj for p=", p
        write(*, "(A)") "================================"
        write(*, *) "threshold=", threshold
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
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,A,I0,A,ES24.16E3)") &
                & "|| M(1)2M(", i, ") - [adjoint M(1)2M(", i, ")]^T ||  /  ", &
                & "|| M(1)2M(", i, ") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
    end do
    !! Check adjoint M2M OZ translation without precomputed matrices. Spheres 2
    !! and 3 are explicitly aligned along OZ axis.
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_ztranslate_adj for p=", p
        write(*, "(A)") "================================"
        write(*, *) "threshold=", threshold
    end if
    do i = 2, 3
        call fmm_m2m_ztranslate_get_mat(-csph(3, i), rsph(1), rsph(i), p, &
            & vscales, vfact, ztranslate_mat)
        call fmm_m2m_ztranslate_use_mat_adj(p, ztranslate_mat, one, &
            & coef(:, 1), zero, coef2(:, i))
        full_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        call fmm_m2m_ztranslate_adj(csph(3, i), rsph(i), rsph(1), p, &
            & vscales, vfact, -two, coef(:, 1), two, coef2(:, i))
        diff_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        v(i) = diff_norm / full_norm
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,A,I0,A,ES24.16E3)") &
                & "|| M(1)2M(", i, ") - [adjoint M(1)2M(", i, ")]^T ||  /  ", &
                & "|| M(1)2M(", i, ") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
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
    if (iprint .gt. 0) then
        write(*, "(A,A,ES24.16E3)") &
            & "|| M(1)2M(6) - [adjoint M(1)2M(6)]^T ||  /  ", &
            & "|| M(1)2M(6) || =", v(6)
    end if
    if (v(6) .gt. threshold) stop 1
    !! Check M2M by reflection
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_reflection for p=", p
        write(*, "(A)") "================================"
        write(*, *) "threshold=", threshold
    end if
    do i = 1, nsph
        call fmm_m2m_reflection(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") &
                & "|| P2M(1) + M(1)2M(", i, ") - P2M(", i, &
                & ") ||  /  || P2M(", i ,") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
        call fmm_m2m_reflection(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, -one, coef(:, 1), two, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") &
                & "|| P2M(1) + M(1)2M(", i, ") - P2M(", i, &
                & ") ||  /  || P2M(", i ,") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
    end do
    !! Check M2M by reflection with help of matrices
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_reflection_get_mat for p=", p
        write(*, "(A,I0)") "Check fmm_m2m_reflection_use_mat for p=", p
        write(*, "(A)") "================================"
        write(*, *) "threshold=", threshold
    end if
    do i = 1, nsph
        call fmm_m2m_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, vfact, transform_mat, ztranslate_mat)
        call fmm_m2m_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, transform_mat, ztranslate_mat, one, coef(:, 1), &
            & zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") &
                & "|| P2M(1) + M(1)2M(", i, ") - P2M(", i, &
                & ") ||  /  || P2M(", i ,") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
        call fmm_m2m_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, transform_mat, ztranslate_mat, -one, coef(:, 1), &
            & two, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") &
                & "|| P2M(1) + M(1)2M(", i, ") - P2M(", i, &
                & ") ||  /  || P2M(", i ,") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
    end do
    !! Check adjoint M2M by reflection with help of matrices
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_reflection_use_mat_adj for p=", p
        write(*, "(A)") "================================"
        write(*, *) "threshold=", threshold
    end if
    do i = 1, nsph
        ! Generate entire matrix of a direct M2M translation
        call fmm_m2m_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, vfact, transform_mat, ztranslate_mat)
        do j = 1, (p+1)**2
            coef2(:, i) = zero
            coef2(j, i) = one
            call fmm_m2m_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
                & rsph(i), p, transform_mat, ztranslate_mat, one, coef2(:, i), &
                & zero, full_mat(:, j))
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
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,A,I0,A,ES24.16E3)") &
                & "|| M(1)2M(", i, ") - [adjoint M(1)2M(", i, ")]^T ||  /  ", &
                & "|| M(1)2M(", i, ") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
    end do
    !! Check adjoint M2M by reflection
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_reflection_adj for p=", p
        write(*, "(A)") "================================"
        write(*, *) "threshold=", threshold
    end if
    do i = 1, nsph
        ! Check with help of fmm_m2m_reflection_use_mat that is checked already
        call fmm_m2m_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, vfact, transform_mat, ztranslate_mat)
        ! Get reference value
        call fmm_m2m_reflection_use_mat_adj(csph(:, i)-csph(:, 1), &
            & rsph(i), rsph(1), p, transform_mat, ztranslate_mat, one, &
            & coef(:, 1), zero, coef2(:, i))
        full_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        ! Subtract value to be checked
        call fmm_m2m_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -two, coef(:, 1), two, coef2(:, i))
        v(i) = dnrm2((p+1)**2, coef2(:, i), 1) / full_norm
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,A,I0,A,ES24.16E3)") &
                & "|| M(1)2M(", i, ") - [adjoint M(1)2M(", i, ")]^T ||  /  ", &
                & "|| M(1)2M(", i, ") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
        ! Check with beta=zero
        call fmm_m2m_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        ! Subtract reference value
        call fmm_m2m_reflection_use_mat_adj(csph(:, i)-csph(:, 1), &
            & rsph(i), rsph(1), p, transform_mat, ztranslate_mat, -one, &
            & coef(:, 1), one, coef2(:, i))
        v(i) = dnrm2((p+1)**2, coef2(:, i), 1) / full_norm
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,A,I0,A,ES24.16E3)") &
                & "|| M(1)2M(", i, ") - [adjoint M(1)2M(", i, ")]^T ||  /  ", &
                & "|| M(1)2M(", i, ") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
    end do
    !! Check M2M by reflection (another implementation)
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_reflection2 for p=", p
        write(*, "(A)") "================================"
        write(*, *) "threshold=", threshold
    end if
    do i = 1, nsph
        call fmm_m2m_reflection2(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") &
                & "|| P2M(1) + M(1)2M(", i, ") - P2M(", i, &
                & ") ||  /  || P2M(", i ,") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
        call fmm_m2m_reflection2(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, -one, coef(:, 1), two, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") &
                & "|| P2M(1) + M(1)2M(", i, ") - P2M(", i, &
                & ") ||  /  || P2M(", i ,") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
    end do
    !! Check adjoint M2M by reflection (another implementation)
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_reflection2_adj for p=", p
        write(*, "(A)") "================================"
        write(*, *) "threshold=", threshold
    end if
    do i = 1, nsph
        ! Get reference value
        call fmm_m2m_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        full_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        ! Subtract value to be checked
        call fmm_m2m_reflection2_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -two, coef(:, 1), two, coef2(:, i))
        v(i) = dnrm2((p+1)**2, coef2(:, i), 1) / full_norm
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,A,I0,A,ES24.16E3)") &
                & "|| M(1)2M(", i, ") - [adjoint M(1)2M(", i, ")]^T ||  /  ", &
                & "|| M(1)2M(", i, ") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
        ! Check with beta=zero
        call fmm_m2m_reflection2_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        ! Subtract reference value
        call fmm_m2m_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -one, coef(:, 1), one, coef2(:, i))
        v(i) = dnrm2((p+1)**2, coef2(:, i), 1) / full_norm
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,A,I0,A,ES24.16E3)") &
                & "|| M(1)2M(", i, ") - [adjoint M(1)2M(", i, ")]^T ||  /  ", &
                & "|| M(1)2M(", i, ") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
    end do
    !! Check M2M by rotation
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_rotation for p=", p
        write(*, "(A)") "================================"
        write(*, *) "threshold=", threshold
    end if
    do i = 1, nsph
        call fmm_m2m_rotation(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") &
                & "|| P2M(1) + M(1)2M(", i, ") - P2M(", i, &
                & ") ||  /  || P2M(", i ,") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
        call fmm_m2m_rotation(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        v(i) = dnrm2(nbasis, coef2(:, i)-coef(:, i), 1) / &
            & dnrm2(nbasis, coef(:, i), 1)
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") &
                & "|| P2M(1) + M(1)2M(", i, ") - P2M(", i, &
                & ") ||  /  || P2M(", i ,") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
    end do
    !! Check adjoint M2M by rotation
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_rotation_adj for p=", p
        write(*, "(A)") "================================"
        write(*, *) "threshold=", threshold
    end if
    do i = 1, nsph
        ! Get reference value
        call fmm_m2m_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        full_norm = dnrm2((p+1)**2, coef2(:, i), 1)
        ! Subtract value to be checked
        call fmm_m2m_rotation_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -two, coef(:, 1), two, coef2(:, i))
        v(i) = dnrm2((p+1)**2, coef2(:, i), 1) / full_norm
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,A,I0,A,ES24.16E3)") &
                & "|| M(1)2M(", i, ") - [adjoint M(1)2M(", i, ")]^T ||  /  ", &
                & "|| M(1)2M(", i, ") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
        ! Check with beta=zero
        call fmm_m2m_rotation_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, one, coef(:, 1), zero, coef2(:, i))
        ! Subtract reference value
        call fmm_m2m_reflection_adj(csph(:, i)-csph(:, 1), rsph(i), rsph(1), &
            & p, vscales, vfact, -one, coef(:, 1), one, coef2(:, i))
        v(i) = dnrm2((p+1)**2, coef2(:, i), 1) / full_norm
        if (iprint .gt. 0) then
            write(*, "(A,I0,A,I0,A,A,I0,A,ES24.16E3)") &
                & "|| M(1)2M(", i, ") - [adjoint M(1)2M(", i, ")]^T ||  /  ", &
                & "|| M(1)2M(", i, ") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
    end do
end subroutine check_m2m

end program test_dd_core

