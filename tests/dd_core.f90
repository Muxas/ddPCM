!> @copyright (c) 2020-2020 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/dd_core.f90
!! Tests for dd_Core module
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2020-12-17

program test_dd_core
use dd_core
implicit none
integer :: p, i
real(kind=rp) :: alpha(4)=(/1d0, 1d-307, 1d-308, 1d+307/)

do i = 1, 4
    do p = 0, 40
        call check_p2m_m2p_m2m_baseline(p, alpha(i))
    end do
end do


contains

! Check P2M, M2M and baseline M2P for spherical harmonics
subroutine check_p2m_m2p_m2m_baseline(p, alpha)
    integer, intent(in) :: p
    real(kind=rp), intent(in) :: alpha
    integer :: nbasis
    real(kind=rp) :: vscales((p+1)**2)
    real(kind=rp), dimension(3) :: src, dst, sph1, sph2, sph3, sph4, sph5
    real(kind=rp) :: sph1_r, sph2_r, sph3_r, sph4_r, sph5_r
    real(kind=rp), dimension((p+1)**2) :: sph1_m, sph2_m, sph3_m, sph4_m, &
        & sph5_m
    real(kind=rp) :: v, v2, v3((p+1)**2)
    real(kind=rp), external :: dnrm2
    ! Scale inputs
    src = alpha * (/5d-1, 4d-1, 3d-1/)
    dst = alpha * (/1d0, 1.1d0, 2.2d0/)
    sph1 = 0d0
    sph1_r = alpha
    sph2 = -alpha
    sph3 = alpha * (/0d0, 0d0, -1d0/)
    sph4 = 0d0
    sph4_r = 3 * alpha
    sph5 = alpha * (/1d-9, 1d-9, -1d0/)
    ! Preliminaries
    write(*, *)
    write(*, "(A,I0,A,ES11.4E3)") "Check P2M and M2P for p=", p, " alpha=", &
        & alpha
    write(*, "(A)") "================================"
    nbasis = (p+1) * (p+1)
    ! Compute sphere radiuses such that they contain each other
    sph2_r = sph1_r + dnrm2(3, sph2-sph1, 1)
    sph3_r = sph1_r + dnrm2(3, sph3-sph1, 1)
    sph5_r = sph1_r + dnrm2(3, sph5-sph1, 1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(p, vscales)
    ! Compute multipole coefficients for given spheres from source particle
    call fmm_p2m(src-sph1, sph1_r, p, vscales, sph1_m)
    call fmm_p2m(src-sph2, sph2_r, p, vscales, sph2_m)
    call fmm_p2m(src-sph3, sph3_r, p, vscales, sph3_m)
    call fmm_p2m(src-sph4, sph4_r, p, vscales, sph4_m)
    call fmm_p2m(src-sph5, sph5_r, p, vscales, sph5_m)
    ! Get potential, spawned by each sphere with its multipole expansion
    v = 1.0 / dnrm2(3, src-dst, 1)
    write(*, *) "v=", v
    call fmm_m2p(dst-sph1, sph1_r, p, vscales, sph1_m, v2)
    write(*,*) "|| P2M(1) + M(1)2P - P2P ||  /  || P2P || =", abs((v2-v)/v)
    call fmm_m2p(dst-sph2, sph2_r, p, vscales, sph2_m, v2)
    write(*,*) "|| P2M(2) + M(2)2P - P2P ||  /  || P2P || =", abs((v2-v)/v)
    call fmm_m2p(dst-sph3, sph3_r, p, vscales, sph3_m, v2)
    write(*,*) "|| P2M(3) + M(3)2P - P2P ||  /  || P2P || =", abs((v2-v)/v)
    call fmm_m2p(dst-sph4, sph4_r, p, vscales, sph4_m, v2)
    write(*,*) "|| P2M(4) + M(4)2P - P2P ||  /  || P2P || =", abs((v2-v)/v)
    call fmm_m2p(dst-sph5, sph5_r, p, vscales, sph5_m, v2)
    write(*,*) "|| P2M(5) + M(5)2P - P2P ||  /  || P2P || =", abs((v2-v)/v)
    ! Check M2M operation
    !write(*, *)
    !write(*, "(A,I0)") "Check baseline M2M for p=", p
    !write(*, "(A)") "================================"
    !v3 = 0
    !call fmm_m2m_baseline(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_m, &
    !    & v3)
    !v = dnrm2(nbasis, sph2_m, 1)
    !v2 = dnrm2(nbasis, sph2_m-v3, 1)
    !write(*,*) "|| P2M(1) + M(1)2M(2) - P2M(2) ||  /  || P2M(2) || =", v2/v
    !v3 = 0
    !call fmm_m2m_baseline(sph1-sph3, sph1_r, sph3_r, p, vscales, sph1_m, &
    !    & v3)
    !v = dnrm2(nbasis, sph3_m, 1)
    !v2 = dnrm2(nbasis, sph3_m-v3, 1)
    !write(*,*) "|| P2M(1) + M(1)2M(3) - P2M(3) ||  /  || P2M(3) || =", v2/v
    !v3 = 0
    !call fmm_m2m_baseline(sph1-sph4, sph1_r, sph4_r, p, vscales, sph1_m, &
    !    & v3)
    !v = dnrm2(nbasis, sph4_m, 1)
    !v2 = dnrm2(nbasis, sph4_m-v3, 1)
    !write(*,*) "|| P2M(1) + M(1)2M(4) - P2M(4) ||  /  || P2M(4) || =", v2/v
end subroutine check_p2m_m2p_m2m_baseline

end program test_dd_core

