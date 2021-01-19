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
integer :: p=100, i
! Some implementation do not work for alpha=1d-307, so range of test values is
! reduced. This can be fixed by enforcing proper order of multiplications like
! a*b*c into a*(b*c).
real(dp) :: alpha(3)=(/1d0, 1d-300, 1d+300/)

! Check arguments to the executable itself
!argc = iargc()
!if ((argc .eq. 0) .or. (argc .gt. 2)) stop "argc"
!call getarg(1, testname)
!if (argc .eq. 2) then
!    call getarg(2, iprint_string)
!    if (iprint_string(1:1) .eq. '0') then
!        iprint = 0
!    else if (iprint_string(1:1) .eq. '1') then
!        iprint = 1
!    else
!        stop "iprint"
!    end if
!else
!    iprint = 0
!end if

! Check correctness of info for valid and invalid input parameters of ddinit
!if (testname .eq. "ddinit") then
    call check_ddinit_args()
!end if

! Check P2M, M2P and M2M operations of the FMM
!if (testname .eq. "p2m+m2p+m2m") then
    do i = 1, size(alpha)
        call check_p2m_m2p_m2m(0, alpha(i), iprint, 5d-1)
        call check_p2m_m2p_m2m(1, alpha(i), iprint, 3d-1)
        call check_p2m_m2p_m2m(p, alpha(i), iprint, 20*epsilon(zero))
    end do
!end if


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
!       fmm_m2m_ztranslate_get_mat
!       fmm_m2m_ztranslate_use_mat
!       fmm_m2m_reflection
!       fmm_m2m_reflection2
!       fmm_m2m_rotation
! List of implicitly checked functions:
!       ylmscale
!       ylmbas
!       trgev
!       polleg
!       fmm_sph_rotate_oz
!       fmm_sph_transform_oxz
! TODO:
!       fmm_m2m_reflection_get_mat
!       fmm_m2m_reflection_use_mat
subroutine check_p2m_m2p_m2m(p, alpha, iprint, threshold)
    integer, intent(in) :: p, iprint
    real(dp), intent(in) :: alpha, threshold
    integer :: nbasis, nseed, i, j, istatus
    integer, parameter :: nsph = 10
    integer, allocatable :: seed(:)
    real(dp) :: vscales((p+1)**2)
    real(dp) :: src0(3), src(3), dst(3), csph(3, nsph), rsph(nsph), v0, v(nsph)
    real(dp), dimension((p+1)**2, nsph) :: msph, m2sph
    real(dp) :: ztranslate_mat((p+1)*(p+2)*(p+3)/6)
    real(dp) :: transform_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(dp) :: z, r1(3, 3), dst1(3)
    real(dp), external :: dnrm2
    ! Scale inputs
    src0 = zero
    src = alpha * (/1d-1, 2d-1, 3d-1/)
    dst = alpha * (/1d0, 1.1d0, 2.2d0/)
    csph(:, 1) = zero
    rsph(1) = alpha
    csph(:, 2) = alpha * (/0d0, 0d0, -1.1d0/)
    rsph(2) = rsph(1) + dnrm2(3, csph(:, 2)-csph(:, 1), 1)
    csph(:, 3) = alpha * (/0d0, 0d0, 7d-1/)
    rsph(3) = rsph(1) + dnrm2(3, csph(:, 3)-csph(:, 1), 1)
    csph(:, 4) = -alpha / two
    rsph(4) = rsph(1) + dnrm2(3, csph(:, 4)-csph(:, 1), 1)
    csph(:, 5) = alpha / two
    rsph(5) = rsph(1) + dnrm2(3, csph(:, 4)-csph(:, 1), 1)
    csph(:, 6) = zero
    rsph(6) = 3 * alpha
    csph(:, 7) = alpha * (/1d-1, 2d-1, 1.1d0/)
    rsph(7) = rsph(1) + dnrm2(3, csph(:, 7)-csph(:, 1), 1)
    csph(:, 8) = alpha * (/4d-1, 2d-9, 1.1d0/)
    rsph(8) = rsph(1) + dnrm2(3, csph(:, 8)-csph(:, 1), 1)
    csph(:, 9) = alpha * (/1.1d0, 0d0, 0d0/)
    rsph(9) = rsph(1) + dnrm2(3, csph(:, 9)-csph(:, 1), 1)
    csph(:, 10) = alpha * (/-4d-1, 1.1d0, 0d0/)
    rsph(10) = rsph(1) + dnrm2(3, csph(:, 10)-csph(:, 1), 1)
    ! Preliminaries
    nbasis = (p+1) * (p+1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(p, vscales)
    ! Init random number generator with a fixed seed
    call random_seed(size=nseed)
    allocate(seed(nseed), stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) "Allocation for random seed failed!"
        stop 1
    end if
    do i = 1, nseed
        seed(i) = 2*i + 1
    end do
    call random_seed(put=seed)
    deallocate(seed, stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) "Deallocation of random seed failed!"
        stop 1
    end if
    !! Check P2M parameters q and beta
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0,A,ES11.4E3)") "Check P2M q and beta params for p=", &
            & p, " alpha=", alpha
        write(*, "(A)") "================================"
    end if
    ! Check beta=zero
    msph(:, 1) = zero
    call fmm_p2m(src, one, alpha, p, vscales, zero, msph(:, 1))
    m2sph(:, 1) = msph(:, 1)
    call fmm_p2m(src, one, alpha, p, vscales, zero, m2sph(:, 1))
    v(1) = dnrm2(nbasis, msph(:, 1)-m2sph(:, 1), 1) / &
        & dnrm2(nbasis, msph(:, 1), 1)
    if (iprint .gt. 0) then
        write(*, *) "P2M beta=zero param check: ", v(1)
    end if
    if (v(1) .gt. threshold) stop 1
    m2sph(:, 1) = msph(:, 1)
    call fmm_p2m(src, one, alpha, p, vscales, one, m2sph(:, 1))
    v(1) = dnrm2(nbasis, two*msph(:, 1)-m2sph(:, 1), 1) / &
        & dnrm2(nbasis, msph(:, 1), 1) / two
    if (iprint .gt. 0) then
        write(*, *) "P2M q and beta params check:", v(1)
    end if
    if (v(1) .gt. threshold) stop 1
    call fmm_p2m(src, -pt5, alpha, p, vscales, -two, m2sph(:, 1))
    v(1) = dnrm2(nbasis, 4.5d0*msph(:, 1)+m2sph(:, 1), 1) / &
        & dnrm2(nbasis, msph(:, 1), 1) / 4.5d0
    if (iprint .gt. 0) then
        write(*, *) "P2M q and beta params check:", v(1)
    end if
    if (v(1) .gt. threshold) stop 1
    !! Check M2P parameters alpha and beta
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,A,I0,A,ES11.4E3)") "Check M2P alpha and beta params ", &
            & "for p=", p, " alpha=", alpha
        write(*, "(A)") "================================"
    end if
    v(1) = zero
    call fmm_m2p(dst, alpha, p, vscales, one, msph(:, 1), zero, v(1))
    v0 = v(1)
    call fmm_m2p(dst, alpha, p, vscales, one, msph(:, 1), zero, v(1))
    if (iprint .gt. 0) then
        write(*, *) "M2P beta=zero param check:", abs(v(1)/v0 - one)
    end if
    if (abs(v(1)/v0 - one) .gt. threshold) stop 1
    call fmm_m2p(dst, alpha, p, vscales, zero, msph(:, 1), -one, v(1))
    if (iprint .gt. 0) then
        write(*, *) "M2P alpha=zero param check:", abs(v(1)/v0 + one)
    end if
    if (abs(v(1)/v0 + one) .gt. threshold) stop 1
    v(1) = v0
    call fmm_m2p(dst, alpha, p, vscales, one, msph(:, 1), one, v(1))
    if (iprint .gt. 0) then
        write(*, *) "M2P alpha and beta params check:", abs(v(1)/v0 - two)
    end if
    if (abs(v(1)/v0 - two) .gt. threshold) stop 1
    call fmm_m2p(dst, alpha, p, vscales, -two, msph(:, 1), pt5, v(1))
    if (iprint .gt. 0) then
        write(*, *) "M2P alpha and beta parameter check:", abs(v(1)/v0 + one)
    end if
    if (abs(v(1)/v0 + one) .gt. threshold) stop 1
    !! Check M2P with a location of P equal to the center of M
    call fmm_m2p(src0, rsph(1), p, vscales, one, msph(:, 1), zero, v0)
    if (iprint .gt. 0) then
        write(*, *) "Check M2P with center of harmonics=particle location"
        write(*, *) "================================="
        write(*, *) "Required result must be 0"
        write(*, *) "Got", v0
    end if
    if (v0 .ne. zero) stop 1
    !! Check P2M+M2P with source particle and multipole harmonics centered at
    !! the origin
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0,A,ES11.4E3)") "Check P2M and M2P for p=", p, &
            & " alpha=", alpha
        write(*, "(A)") "================================"
    end if
    call fmm_p2m(src0, one, alpha, p, vscales, zero, msph(:, 1))
    call fmm_m2p(dst, alpha, p, vscales, one, msph(:, 1), zero, v(1))
    v0 = one / dnrm2(3, dst, 1)
    v(1) = abs((v(1)-v0) / v0)
    if (iprint .gt. 0) then
        write(*, *) "v=", v0, "threshold=", threshold
        write(*, "(A,A,ES24.16E3)") "|| P2M(0) + M(0)2P - P2P ||  /  ", &
            & "|| P2P || =", v(1)
    end if
    if (v(1) .gt. threshold) stop 1
    !! Check P2M+M2P for predefined spheres
    ! Compute multipole coefficients for given spheres from source particle
    do i = 1, nsph
        call fmm_p2m(src-csph(:, i), one, rsph(i), p, vscales, zero, &
            & msph(:, i))
    end do
    ! Get potential, spawned by each sphere with its multipole expansion
    v0 = one / dnrm2(3, src-dst, 1)
    do i = 1, nsph
        call fmm_m2p(dst-csph(:, i), rsph(i), p, vscales, one, msph(:, i), &
            & zero, v(i))
    end do
    ! Finally check p2m+m2p
    v = abs((v-v0) / v0)
    if (iprint .gt. 0) then
        write(*, *) "v=", v0, "threshold=", threshold
        do i = 1, nsph
            write(*, "(A,I0,A,I0,A,ES24.16E3)") "|| P2M(", i, ") + M(", i, &
                & ")2P - P2P ||  /  || P2P || =", v(i)
        end do
    end if
    if (maxval(v) .gt. threshold) stop 1
    !! Check orthogonal transformation of system of coordinates
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A)") "Check coord_reflect_get_mat"
        write(*, "(A,I0)") "Check fmm_sph_transform for p=", p
        write(*, "(A,I0)") "Check fmm_sph_transform_get_mat for p=", p
        write(*, "(A,I0)") "Check fmm_sph_transform_use_mat for p=", p
        write(*, "(A)") "================================"
    end if
    do i = 1, nsph
        call coord_reflect_get_mat(dst-csph(:, i), z, r1)
        if (iprint .gt. 0) then
            write(*, "(A)") "Reflection from new (y,z,x) to new (y,z,x)"
            write(*, "(3ES25.16E3)") r1(1, 1), r1(1, 2), r1(1, 3)
            write(*, "(3ES25.16E3)") r1(2, 1), r1(2, 2), r1(2, 3)
            write(*, "(3ES25.16E3)") r1(3, 1), r1(3, 2), r1(3, 3)
        end if
        dst1 = zero
        do j = 1, 3
            dst1(j) = (dst(1)-csph(1, i))*r1(3, j) + &
                & (dst(2)-csph(2, i))*r1(1, j) + &
                & (dst(3)-csph(3, i))*r1(2, j)
        end do
        ! Convert (y,z,x) to (x,y,z)
        v(1) = dst1(1)
        dst1(1) = dst1(3)
        dst1(3) = dst1(2)
        dst1(2) = v(1)
        if (iprint .gt. 0) then
            write(*, "(A)") "Reflected vector (x,y,z) must be aligned along OZ"
            write(*, "(3ES25.16E3)") dst1(1), dst1(2), dst1(3)
        end if
        ! Check that x and y coordinates are zeros compared to z coordinate
        v(1) = abs(dst1(1) / dst1(3))
        v(2) = abs(dst1(2) / dst1(3))
        if (max(v(1), v(2)) .gt. threshold) stop 1
        ! Check that dst1(3) is equal to z output of coord_reflect_get_mat
        if (abs((z-dst1(3))/z) .gt. threshold) stop 1
        ! Check that reflected spherical harmonics result in the same potential
        call fmm_sph_transform(p, r1, -pt5, msph(:, i), zero, m2sph(:, i))
        call fmm_m2p(dst-csph(:, i), rsph(i), p, vscales, one, msph(:, i), &
            & zero, v0)
        dst1(1) = zero
        dst1(2) = zero
        call fmm_m2p(dst1, rsph(i), p, vscales, -two, m2sph(:, i), zero, v(i))
        if (iprint .gt. 0) then
            write(*, "(A,A,ES24.16E3,A,ES23.16E3,A,I0)") "fmm_sph_transform ", &
                & "relative error is ", abs((v(i)-v0)/v0), " threshold=", &
                & threshold
        end if
        if (abs((v(i)-v0)/v0) .gt. threshold) stop 1
        call fmm_sph_transform(p, r1, pt5, msph(:, i), two, m2sph(:, i))
        dst1(1) = zero
        dst1(2) = zero
        call fmm_m2p(dst1, rsph(i), p, vscales, -two, m2sph(:, i), zero, v(i))
        if (iprint .gt. 0) then
            write(*, "(A,A,ES24.16E3,A,ES23.16E3,A,I0)") "fmm_sph_transform ", &
                & "relative error is ", abs((v(i)-v0)/v0), " threshold=", &
                & threshold
        end if
        if (abs((v(i)-v0)/v0) .gt. threshold) stop 1
        ! Check transformation of harmonics by corresponding matrix
        call fmm_sph_transform_get_mat(p, r1, transform_mat)
        call fmm_sph_transform_use_mat(p, transform_mat, one, msph(:, i), &
            & zero, m2sph(:, i))
        call fmm_m2p(dst1, rsph(i), p, vscales, one, m2sph(:, i), zero, v(i))
        if (iprint .gt. 0) then
            write(*, "(A,A,ES24.16E3,A,ES23.16E3)") &
                & "fmm_sph_transform_get_mat+fmm_sph_transform_use_mat ", &
                & "relative error is ", abs((v(i)-v0)/v0), " threshold=", &
                & threshold
        end if
        if (abs((v(i)-v0)/v0) .gt. threshold) stop 1
        call fmm_sph_transform_use_mat(p, transform_mat, -one, msph(:, i), &
            & two, m2sph(:, i))
        call fmm_m2p(dst1, rsph(i), p, vscales, one, m2sph(:, i), zero, v(i))
        if (iprint .gt. 0) then
            write(*, "(A,A,ES24.16E3,A,ES23.16E3)") &
                & "fmm_sph_transform_get_mat+fmm_sph_transform_use_mat ", &
                & "relative error is ", abs((v(i)-v0)/v0), " threshold=", &
                & threshold
        end if
        if (abs((v(i)-v0)/v0) .gt. threshold) stop 1
    end do
    !! Check M2M by OZ translation Spheres 2 and 3 are explicitly aligned along
    !! OZ axis
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_ztranslate for p=", p
        write(*, "(A,I0)") "Check fmm_m2m_ztranslate_get_mat for p=", p
        write(*, "(A,I0)") "Check fmm_m2m_ztranslate_use_mat for p=", p
        write(*, "(A)") "================================"
    end if
    do i = 2, 3
        call fmm_m2m_ztranslate(-csph(3, i), rsph(1), rsph(i), p, vscales, &
            & one, msph(:, 1), zero, m2sph(:, i))
        v(i) = dnrm2(nbasis, m2sph(:, i)-msph(:, i), 1) / &
            & dnrm2(nbasis, msph(:, i), 1)
        if (iprint .gt. 0) then
            write(*, *) "threshold=", threshold
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") "|| P2M(1) + M(1)2M(", &
                & i, ") - P2M(", i, ") ||  /  || P2M(", i ,") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
        call fmm_m2m_ztranslate(-csph(3, i), rsph(1), rsph(i), p, vscales, &
            & -one, msph(:, 1), two, m2sph(:, i))
        v(i) = dnrm2(nbasis, m2sph(:, i)-msph(:, i), 1) / &
            & dnrm2(nbasis, msph(:, i), 1)
        if (iprint .gt. 0) then
            write(*, *) "threshold=", threshold
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") "|| P2M(1) + M(1)2M(", &
                & i, ") - P2M(", i, ") ||  /  || P2M(", i ,") || =", v(i)
        end if
        if (v(i) .gt. threshold) stop 1
        ! Check transformation of harmonics by corresponding matrix
        call fmm_m2m_ztranslate_get_mat(-csph(3, i), rsph(1), rsph(i), p, &
            & vscales, ztranslate_mat)
        call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, one, msph(:, 1), zero, &
            & m2sph(:, i))
        v(i) = dnrm2(nbasis, m2sph(:, i)-msph(:, i), 1) / &
            & dnrm2(nbasis, msph(:, i), 1)
        if (iprint .gt. 0) then
            write(*, "(A,A,ES24.16E3,A,ES23.16E3)") &
                & "fmm_m2m_ztranslate_get_mat+fmm_m2m_ztranslate_use_mat ", &
                & "relative error is ", v(i), " threshold=", threshold
        end if
        if (v(i) .gt. threshold) stop 1
        call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, -one, msph(:, 1), two, &
            & m2sph(:, i))
        v(i) = dnrm2(nbasis, m2sph(:, i)-msph(:, i), 1) / &
            & dnrm2(nbasis, msph(:, i), 1)
        if (iprint .gt. 0) then
            write(*, "(A,A,ES24.16E3,A,ES23.16E3)") &
                & "fmm_m2m_ztranslate_get_mat+fmm_m2m_ztranslate_use_mat ", &
                & "relative error is ", v(i), " threshold=", threshold
        end if
        if (v(i) .gt. threshold) stop 1
    end do
    !! Check that attempt to use M2M matrix with 0 shift raises NaNs
    call fmm_m2m_ztranslate_get_mat(zero, rsph(1), rsph(6), p, vscales, &
        & ztranslate_mat)
    call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, one, msph(:, 1), zero, &
        & m2sph(:, 6))
    !! Check M2M with zero OZ translation. Sphere 6 is explicitly set for this
    call fmm_m2m_ztranslate(zero, rsph(1), rsph(6), p, vscales, one, &
        & msph(:, 1), zero, m2sph(:, 6))
    v(6) = dnrm2(nbasis, m2sph(:, 6)-msph(:, 6), 1) / &
        & dnrm2(nbasis, msph(:, 6), 1)
    if (iprint .gt. 0) then
        write(*, *) "threshold=", threshold
        write(*, "(A,A,ES24.16E3)") "|| P2M(1) + M(1)2M(6) - P2M(6) ||  /  ", &
            & "|| P2M(6) || =", v(6)
    end if
    if (v(6) .gt. threshold) stop 1
    ! M2M with 0 shift is simply scaling
    call fmm_m2m_scale(rsph(1), rsph(6), p, one, msph(:, 1), zero, m2sph(:, 6))
    v(6) = dnrm2(nbasis, m2sph(:, 6)-msph(:, 6), 1) / &
        & dnrm2(nbasis, msph(:, 6), 1)
    if (iprint .gt. 0) then
        write(*, "(A,ES24.16E3,A,ES23.16E3)") &
            & "fmm_m2m_scale relative error is ", v(6), " threshold=", zero
    end if
    call fmm_m2m_scale(rsph(1), rsph(6), p, -one, msph(:, 1), two, m2sph(:, 6))
    v(6) = dnrm2(nbasis, m2sph(:, 6)-msph(:, 6), 1) / &
        & dnrm2(nbasis, msph(:, 6), 1)
    if (iprint .gt. 0) then
        write(*, "(A,ES24.16E3,A,ES23.16E3)") &
            & "fmm_m2m_scale relative error is ", v(6), " threshold=", zero
    end if
    !! Check M2M by reflection
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_reflection for p=", p
        write(*, "(A,I0)") "Check fmm_m2m_reflection_get_mat for p=", p
        write(*, "(A,I0)") "Check fmm_m2m_reflection_use_mat for p=", p
        write(*, "(A)") "================================"
    end if
    do i = 1, nsph
        call fmm_m2m_reflection(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, one, msph(:, 1), zero, m2sph(:, i))
        v(i) = dnrm2(nbasis, m2sph(:, i)-msph(:, i), 1) / &
            & dnrm2(nbasis, msph(:, i), 1)
    end do
    if (iprint .gt. 0) then
        write(*, *) "threshold=", threshold
        do i = 1, nsph
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") "|| P2M(1) + M(1)2M(", &
                & i, ") - P2M(", i, ") ||  /  || P2M(", i ,") || =", v(i)
        end do
    end if
    if (maxval(v) .gt. threshold) stop 1
    do i = 1, nsph
        call fmm_m2m_reflection_get_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, vscales, transform_mat, ztranslate_mat)
        call fmm_m2m_reflection_use_mat(csph(:, 1)-csph(:, i), rsph(1), &
            & rsph(i), p, transform_mat, ztranslate_mat, one, msph(:, 1), &
            & zero, m2sph(:, i))
        v(i) = dnrm2(nbasis, m2sph(:, i)-msph(:, i), 1) / &
            & dnrm2(nbasis, msph(:, i), 1)
    end do
    if (iprint .gt. 0) then
        write(*, *) "threshold=", threshold
        do i = 1, nsph
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") "|| P2M(1) + M(1)2M(", &
                & i, ") - P2M(", i, ") ||  /  || P2M(", i ,") || =", v(i)
        end do
    end if
    if (maxval(v) .gt. threshold) stop 1
    !! Check M2M by reflection (another implementation)
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_reflection2 for p=", p
        write(*, "(A)") "================================"
    end if
    do i = 1, nsph
        call fmm_m2m_reflection2(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, one, msph(:, 1), zero, m2sph(:, i))
        v(i) = dnrm2(nbasis, m2sph(:, i)-msph(:, i), 1) / &
            & dnrm2(nbasis, msph(:, i), 1)
    end do
    if (iprint .gt. 0) then
        write(*, *) "threshold=", threshold
        do i = 1, nsph
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") "|| P2M(1) + M(1)2M(", &
                & i, ") - P2M(", i, ") ||  /  || P2M(", i ,") || =", v(i)
        end do
    end if
    if (maxval(v) .gt. threshold) stop 1
    !! Check M2M by rotation
    if (iprint .gt. 0) then
        write(*, *)
        write(*, "(A,I0)") "Check fmm_m2m_rotation for p=", p
        write(*, "(A)") "================================"
    end if
    do i = 1, nsph
        call fmm_m2m_rotation(csph(:, 1)-csph(:, i), rsph(1), rsph(i), p, &
            & vscales, one, msph(:, 1), zero, m2sph(:, i))
        v(i) = dnrm2(nbasis, m2sph(:, i)-msph(:, i), 1) / &
            & dnrm2(nbasis, msph(:, i), 1)
    end do
    if (iprint .gt. 0) then
        write(*, *) "threshold=", threshold
        do i = 1, nsph
            write(*, "(A,I0,A,I0,A,I0,A,ES24.16E3)") "|| P2M(1) + M(1)2M(", &
                & i, ") - P2M(", i, ") ||  /  || P2M(", i ,") || =", v(i)
        end do
    end if
    if (maxval(v) .gt. threshold) stop 1
end subroutine check_p2m_m2p_m2m

end program test_dd_core

