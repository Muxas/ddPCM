program test_pcm_fmm
!use ddcosmo
implicit none
integer :: p1=6, p2=10, p3=15, ngrid=6

call check_p2m_m2p_m2m_baseline(p1)
call check_p2m_m2p_m2m_baseline(p2)
call check_p2m_m2p_m2m_baseline(p3)
call check_m2l_l2p_l2l_baseline(p1, p2)
call check_m2l_l2p_l2l_baseline(p2, p1)
call check_m2l_l2p_l2l_baseline(p1, p3)
call check_m2l_l2p_l2l_baseline(p3, p1)
call check_m2l_l2p_l2l_baseline(p2, p3)
call check_m2l_l2p_l2l_baseline(p3, p2)
call check_m2m_improved(p1)
call check_m2m_improved(p2)
call check_m2m_improved(p3)
call check_l2l_improved(p1)
call check_l2l_improved(p2)
call check_l2l_improved(p3)
call check_m2l_improved(p1, p2)
call check_m2l_improved(p2, p1)
call check_m2l_improved(p1, p3)
call check_m2l_improved(p3, p1)
call check_m2l_improved(p2, p3)
call check_m2l_improved(p3, p2)

end program test_pcm_fmm

! Check P2M, M2M and baseline M2P for spherical harmonics
subroutine check_p2m_m2p_m2m_baseline(p)
    use pcm_fmm
    implicit none
    integer, intent(in) :: p
    integer :: nbasis
    integer, parameter :: ngrid=1
    real(kind=8) :: vscales((p+1)*(p+1)), w(ngrid), grid(3, ngrid)
    real(kind=8) :: vgrid((p+1)*(p+1), ngrid)
    real(kind=8) :: src(3)=(/0.5, 0.4, 0.3/)
    real(kind=8) :: dst(3)=(/1.0, 1.1, 2.2/)
    real(kind=8) :: sph1(3)=(/0.0, 0.0, 0.0/), sph1_r=1.0
    real(kind=8) :: sph2(3)=(/-1.0, -1.0, -1.0/), sph2_r
    real(kind=8) :: sph3(3)=(/0.0, 0.0, -1.0/), sph3_r
    real(kind=8) :: sph4(3)=(/0.0, 0.0, 0.0/), sph4_r=3
    real(kind=8) :: sph1_m((p+1)*(p+1))
    real(kind=8) :: sph2_m((p+1)*(p+1))
    real(kind=8) :: sph3_m((p+1)*(p+1))
    real(kind=8) :: sph4_m((p+1)*(p+1))
    real(kind=8) :: v, v2, v3((p+1)*(p+1))
    real(kind=8), external :: dnrm2
    ! Preliminaries
    write(*, *)
    write(*, "(A,I0)") "Check P2M and M2P for p=", p
    write(*, "(A)") "================================"
    nbasis = (p+1) * (p+1)
    ! Compute sphere radiuses such that they contain each other
    sph2_r = sph1_r + dnrm2(3, sph2-sph1, 1)
    sph3_r = sph1_r + dnrm2(3, sph3-sph1, 1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics and number of Lebedev grid points
    call init_globals(0, p, vscales, ngrid, w, grid, vgrid)
    ! Compute multipole coefficients for given spheres from source particle
    call fmm_p2m(src-sph1, sph1_r, p, vscales, sph1_m)
    call fmm_p2m(src-sph2, sph2_r, p, vscales, sph2_m)
    call fmm_p2m(src-sph3, sph3_r, p, vscales, sph3_m)
    call fmm_p2m(src-sph4, sph4_r, p, vscales, sph4_m)
    ! Get potential, spawned by each sphere with its multipole expansion
    v = 1.0 / dnrm2(3, src-dst, 1)
    call fmm_m2p(dst-sph1, sph1_r, p, vscales, sph1_m, v2)
    write(*,*) "|| P2M(1) + M(1)2P - P2P ||  /  || P2P || =", abs((v2-v)/v)
    call fmm_m2p(dst-sph2, sph2_r, p, vscales, sph2_m, v2)
    write(*,*) "|| P2M(2) + M(2)2P - P2P ||  /  || P2P || =", abs((v2-v)/v)
    call fmm_m2p(dst-sph3, sph3_r, p, vscales, sph3_m, v2)
    write(*,*) "|| P2M(3) + M(3)2P - P2P ||  /  || P2P || =", abs((v2-v)/v)
    call fmm_m2p(dst-sph4, sph4_r, p, vscales, sph4_m, v2)
    write(*,*) "|| P2M(4) + M(4)2P - P2P ||  /  || P2P || =", abs((v2-v)/v)
    ! Check M2M operation
    write(*, *)
    write(*, "(A,I0)") "Check baseline M2M for p=", p
    write(*, "(A)") "================================"
    v3 = 0
    call fmm_m2m_baseline(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_m, &
        & v3)
    v = dnrm2(nbasis, sph2_m, 1)
    v2 = dnrm2(nbasis, sph2_m-v3, 1)
    write(*,*) "|| P2M(1) + M(1)2M(2) - P2M(2) ||  /  || P2M(2) || =", v2/v
    v3 = 0
    call fmm_m2m_baseline(sph1-sph3, sph1_r, sph3_r, p, vscales, sph1_m, &
        & v3)
    v = dnrm2(nbasis, sph3_m, 1)
    v2 = dnrm2(nbasis, sph3_m-v3, 1)
    write(*,*) "|| P2M(1) + M(1)2M(3) - P2M(3) ||  /  || P2M(3) || =", v2/v
    v3 = 0
    call fmm_m2m_baseline(sph1-sph4, sph1_r, sph4_r, p, vscales, sph1_m, &
        & v3)
    v = dnrm2(nbasis, sph4_m, 1)
    v2 = dnrm2(nbasis, sph4_m-v3, 1)
    write(*,*) "|| P2M(1) + M(1)2M(4) - P2M(4) ||  /  || P2M(4) || =", v2/v
end subroutine check_p2m_m2p_m2m_baseline

! Check L2P and baseline L2L and M2L for spherical harmonics
subroutine check_m2l_l2p_l2l_baseline(pm, pl)
    use pcm_fmm
    implicit none
    integer, intent(in) :: pm, pl
    integer :: nbasism, nbasisl
    integer, parameter :: ngrid=1
    real(kind=8) :: vscales((pm+pl+1)*(pm+pl+1)), w(ngrid), grid(3, ngrid)
    real(kind=8) :: vgrid((pl+1)*(pl+1), ngrid)
    real(kind=8) :: src(3)=(/0.5, 0.4, 0.3/)
    real(kind=8) :: dst(3)=(/-1.0, -1.1, -2.2/)
    ! sph1 contains multipole expansion
    real(kind=8) :: sph1(3)=(/0.0, 0.0, 0.0/), sph1_r=2.0
    ! sph2, sph3, sph4 and sph5 contain local expansion
    real(kind=8) :: sph2(3)=(/-1.0, -1.0, -2.0/), sph2_r=3
    real(kind=8) :: sph3(3)=(/0.0, 0.0, -2.2/), sph3_r=1.6
    real(kind=8) :: sph4(3)=(/-1.0, -1.0, -2.2/), sph4_r=0.6
    real(kind=8) :: sph5(3)=(/-1.0, -1.0, -2.0/), sph5_r=0.6
    real(kind=8) :: sph1_m((pm+1)*(pm+1))
    real(kind=8) :: sph2_l((pl+1)*(pl+1))
    real(kind=8) :: sph3_l((pl+1)*(pl+1))
    real(kind=8) :: sph4_l((pl+1)*(pl+1))
    real(kind=8) :: sph5_l((pl+1)*(pl+1))
    real(kind=8) :: v, v2, v3((pl+1)*(pl+1))
    real(kind=8), external :: dnrm2
    ! Preliminaries
    write(*,*) "Check L2P and baseline M2L"
    write(*, *)
    write(*, "(A,I0,A,I0)") "Check L2P and baseline M2L for pm=", pm, &
        & " and pl=", pl
    write(*, "(A)") "================================"
    nbasism = (pm+1) * (pm+1)
    nbasisl = (pl+1) * (pl+1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics and number of Lebedev grid points
    call init_globals(pm, pl, vscales, ngrid, w, grid, vgrid)
    ! Compute multipole coefficients for sph1 from given source particle
    call fmm_p2m(src-sph1, sph1_r, pm, vscales, sph1_m)
    ! Check chain M2L+L2P
    v = 1.0 / dnrm2(3, src-dst, 1)
    sph2_l = 0
    call fmm_m2l_baseline(sph1-sph2, sph1_r, sph2_r, pm, pl, vscales, sph1_m, &
        & sph2_l)
    call fmm_l2p(dst-sph2, sph2_r, pl, vscales, sph2_l, v2)
    write(*,*) "|| P2M(1) + M(1)2L(2) + L(2)2P - P2P ||  /  || P2P || =", &
        & abs((v2-v)/v)
    sph3_l = 0
    call fmm_m2l_baseline(sph1-sph3, sph1_r, sph3_r, pm, pl, vscales, sph1_m, &
        & sph3_l)
    call fmm_l2p(dst-sph3, sph3_r, pl, vscales, sph3_l, v2)
    write(*,*) "|| P2M(1) + M(1)2L(3) + L(3)2P - P2P ||  /  || P2P || =", &
        & abs((v2-v)/v)
    sph4_l = 0
    call fmm_m2l_baseline(sph1-sph4, sph1_r, sph4_r, pm, pl, vscales, sph1_m, &
        & sph4_l)
    call fmm_l2p(dst-sph4, sph4_r, pl, vscales, sph4_l, v2)
    write(*,*) "|| P2M(1) + M(1)2L(4) + L(4)2P - P2P ||  /  || P2P || =", &
        & abs((v2-v)/v)
    sph5_l = 0
    call fmm_m2l_baseline(sph1-sph5, sph1_r, sph5_r, pm, pl, vscales, sph1_m, &
        & sph5_l)
    call fmm_l2p(dst-sph5, sph5_r, pl, vscales, sph5_l, v2)
    write(*,*) "|| P2M(1) + M(1)2L(5) + L(5)2P - P2P ||  /  || P2P || =", &
        & abs((v2-v)/v)
    ! Check L2L
    write(*, *)
    write(*, "(A,I0)") "Check baseline L2L for p=", pl
    write(*, "(A)") "================================"
    call fmm_l2p(dst-sph2, sph2_r, pl, vscales, sph2_l, v)
    v3 = 0
    call fmm_l2l_baseline(sph2-sph3, sph2_r, sph3_r, pl, vscales, sph2_l, &
        & v3)
    call fmm_l2p(dst-sph3, sph3_r, pl, vscales, v3, v2)
    write(*,*) "|| L(2)2L(3) + L(3)2P - L(2)2P ||  /  || L(2)2P || =", &
        & abs((v2-v)/v)
    v3 = 0
    call fmm_l2l_baseline(sph2-sph4, sph2_r, sph4_r, pl, vscales, sph2_l, &
        & v3)
    call fmm_l2p(dst-sph4, sph4_r, pl, vscales, v3, v2)
    write(*,*) "|| L(2)2L(4) + L(4)2P - L(2)2P ||  /  || L(2)2P || =", &
        & abs((v2-v)/v)
    v3 = 0
    call fmm_l2l_baseline(sph2-sph5, sph2_r, sph5_r, pl, vscales, sph2_l, &
        & v3)
    call fmm_l2p(dst-sph5, sph5_r, pl, vscales, v3, v2)
    write(*,*) "|| L(2)2L(5) + L(5)2P - L(2)2P ||  /  || L(2)2P || =", &
        & abs((v2-v)/v)
end subroutine check_m2l_l2p_l2l_baseline

subroutine check_m2m_improved(p)
    use pcm_fmm
    implicit none
    integer, intent(in) :: p
    integer, parameter :: ngrid=1
    integer :: nbasis, i
    real(kind=8) :: vscales((p+1)*(p+1)), w(ngrid), grid(3, ngrid)
    real(kind=8) :: vgrid((p+1)*(p+1), ngrid)
    real(kind=8) :: sph1(3)=(/0.0, 0.0, 0.0/), sph1_r=1
    real(kind=8) :: sph2(3)=(/1.0, 4.0, -10.0/), sph2_r
    real(kind=8) :: sph3(3)=(/0.0, 0.0, 4.0/), sph3_r
    real(kind=8) :: sph4(3)=(/0.0, 0.0, 0.0/), sph4_r=2
    real(kind=8), external :: dnrm2
    real(kind=8) :: sph1_m((p+1)*(p+1))
    real(kind=8) :: sph2_m((p+1)*(p+1))
    real(kind=8) :: sph3_m((p+1)*(p+1))
    real(kind=8) :: sph4_m((p+1)*(p+1))
    real(kind=8) :: v((p+1)*(p+1))
    real(kind=8) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(kind=8) :: ztrans_mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8) :: start, finish
    ! Preliminaries
    nbasis = (p+1) * (p+1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics and number of Lebedev grid points
    call init_globals(0, p, vscales, ngrid, w, grid, vgrid)
    ! Compute radiuses
    sph2_r = sph1_r + dnrm2(3, sph2-sph1, 1)
    sph3_r = sph1_r + dnrm2(3, sph3-sph1, 1)
    ! Init values of spherical harmonics
    sph1_m = 1
    ! Get coefficients by baseline M2M
    sph2_m = 0
    call fmm_m2m_baseline(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_m, &
        & sph2_m)
    sph3_m = 0
    call fmm_m2m_baseline(sph1-sph3, sph1_r, sph3_r, p, vscales, sph1_m, &
        & sph3_m)
    sph4_m = 0
    call fmm_m2m_baseline(sph1-sph4, sph1_r, sph4_r, p, vscales, sph1_m, &
        & sph4_m)
    ! Check fmm_m2m_reflection
    write(*, *)
    write(*, "(A,I0)") "Check fmm_m2m_reflection vs fmm_m2m_baseline for p=", p
    write(*, "(A)") "======================================================="
    v = 0
    call fmm_m2m_reflection(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2M(2) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph2_m, 1) / dnrm2(nbasis, sph2_m, 1)
    v = 0
    call fmm_m2m_reflection(sph1-sph3, sph1_r, sph3_r, p, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2M(3) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph3_m, 1) / dnrm2(nbasis, sph3_m, 1)
    v = 0
    call fmm_m2m_reflection(sph1-sph4, sph1_r, sph4_r, p, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2M(4) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph4_m, 1) / dnrm2(nbasis, sph4_m, 1)
    ! Check fmm_m2m_reflection3
    write(*, *)
    write(*, "(A,A,I0)") "Check fmm_m2m_reflection3 vs fmm_m2m_baseline ", &
        & "for p=", p
    write(*, "(A)") "======================================================="
    v = 0
    call fmm_m2m_reflection3(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2M(2) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph2_m, 1) / dnrm2(nbasis, sph2_m, 1)
    v = 0
    call fmm_m2m_reflection3(sph1-sph3, sph1_r, sph3_r, p, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2M(3) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph3_m, 1) / dnrm2(nbasis, sph3_m, 1)
    v = 0
    call fmm_m2m_reflection3(sph1-sph4, sph1_r, sph4_r, p, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2M(4) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph4_m, 1) / dnrm2(nbasis, sph4_m, 1)
    ! Check fmm_m2m_fast
    write(*, *)
    write(*, "(A,I0)") "Check fmm_m2m_fast vs fmm_m2m_baseline for p=", p
    write(*, "(A)") "======================================================="
    v = 0
    call fmm_m2m_fast(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2M(2) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph2_m, 1) / dnrm2(nbasis, sph2_m, 1)
    v = 0
    call fmm_m2m_fast(sph1-sph3, sph1_r, sph3_r, p, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2M(3) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph3_m, 1) / dnrm2(nbasis, sph3_m, 1)
    v = 0
    call fmm_m2m_fast(sph1-sph4, sph1_r, sph4_r, p, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2M(4) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph4_m, 1) / dnrm2(nbasis, sph4_m, 1)
    ! Check fmm_m2m_test_mat
    write(*, *)
    write(*, "(A,I0)") "Check fmm_m2m_test_mat vs fmm_m2m_baseline for p=", p
    write(*, "(A)") "======================================================="
    v = 0
    call fmm_m2m_test_mat(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2M(2) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph2_m, 1) / dnrm2(nbasis, sph2_m, 1)
    v = 0
    call fmm_m2m_test_mat(sph1-sph3, sph1_r, sph3_r, p, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2M(3) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph3_m, 1) / dnrm2(nbasis, sph3_m, 1)
    v = 0
    call fmm_m2m_test_mat(sph1-sph4, sph1_r, sph4_r, p, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2M(4) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph4_m, 1) / dnrm2(nbasis, sph4_m, 1)
    ! Check fmm_m2m_get_mat + fmm_m2m_use_mat
    write(*, *)
    write(*, "(A,A,I0)") "Check fmm_m2m_get_mat+fmm_m2m_use_mat vs ", &
        & "fmm_m2m_baseline for p=", p
    write(*, "(A)") "======================================================="
    v = 0
    call fmm_m2m_get_mat(sph1-sph2, sph1_r, sph2_r, p, vscales, reflect_mat, &
        & ztrans_mat)
    call fmm_m2m_use_mat(sph1-sph2, sph1_r, sph2_r, p, reflect_mat, &
        & ztrans_mat, sph1_m, v)
    write(*,*) "|| M(1)2M(2) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph2_m, 1) / dnrm2(nbasis, sph2_m, 1)
    v = 0
    call fmm_m2m_get_mat(sph1-sph3, sph1_r, sph3_r, p, vscales, reflect_mat, &
        & ztrans_mat)
    call fmm_m2m_use_mat(sph1-sph3, sph1_r, sph3_r, p, reflect_mat, &
        & ztrans_mat, sph1_m, v)
    write(*,*) "|| M(1)2M(3) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph3_m, 1) / dnrm2(nbasis, sph3_m, 1)
    v = 0
    call fmm_m2m_get_mat(sph1-sph4, sph1_r, sph4_r, p, vscales, reflect_mat, &
        & ztrans_mat)
    call fmm_m2m_use_mat(sph1-sph4, sph1_r, sph4_r, p, reflect_mat, &
        & ztrans_mat, sph1_m, v)
    write(*,*) "|| M(1)2M(4) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph4_m, 1) / dnrm2(nbasis, sph4_m, 1)
    ! Now check performance
    write(*, *)
    write(*, "(A,I0)") "Check performance of M2M for p=", p
    write(*, "(A)") "======================================================="
    call cpu_time(start)
    do i = 1, 1000
        call fmm_m2m_baseline(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_m, &
            & sph2_m)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_m2m_baseline:", (finish-start)/1000, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 1000
        call fmm_m2m_reflection(sph1-sph2, sph1_r, sph2_r, p, vscales, &
            & sph1_m, sph2_m)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_m2m_reflection:", (finish-start)/1000, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 1000
        call fmm_m2m_reflection3(sph1-sph2, sph1_r, sph2_r, p, vscales, &
            & sph1_m, sph2_m)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_m2m_reflection3:", (finish-start)/1000, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 1000
        call fmm_m2m_fast(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_m, &
            & sph2_m)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_m2m_fast:", (finish-start)/1000, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 1000
        call fmm_m2m_test_mat(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_m, &
            & sph2_m)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_m2m_test_mat:", (finish-start)/1000, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 1000
        call fmm_m2m_get_mat(sph1-sph2, sph1_r, sph2_r, p, vscales, &
            & reflect_mat, ztrans_mat)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_m2m_get_mat:", (finish-start)/1000, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 1000
        call fmm_m2m_use_mat(sph1-sph2, sph1_r, sph2_r, p, reflect_mat, &
            & ztrans_mat, sph1_m, v)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_m2m_use_mat:", (finish-start)/1000, &
        & " seconds"
end subroutine check_m2m_improved

subroutine check_l2l_improved(p)
    use pcm_fmm
    implicit none
    integer, intent(in) :: p
    integer, parameter :: ngrid=1
    integer :: nbasis, i
    real(kind=8) :: vscales((p+1)*(p+1)), w(ngrid), grid(3, ngrid)
    real(kind=8) :: vgrid((p+1)*(p+1), ngrid)
    real(kind=8) :: sph1(3)=(/0.0, 0.0, 0.0/), sph1_r=1
    real(kind=8) :: sph2(3)=(/1.0, 4.0, -10.0/), sph2_r
    real(kind=8) :: sph3(3)=(/0.0, 0.0, 4.0/), sph3_r
    real(kind=8) :: sph4(3)=(/0.0, 0.0, 0.0/), sph4_r=2
    real(kind=8), external :: dnrm2
    real(kind=8) :: sph1_l((p+1)*(p+1))
    real(kind=8) :: sph2_l((p+1)*(p+1))
    real(kind=8) :: sph3_l((p+1)*(p+1))
    real(kind=8) :: sph4_l((p+1)*(p+1))
    real(kind=8) :: v((p+1)*(p+1))
    real(kind=8) :: reflect_mat((p+1)*(2*p+1)*(2*p+3)/3)
    real(kind=8) :: ztrans_mat((p+1)*(p+2)*(p+3)/6)
    real(kind=8) :: start, finish
    ! Preliminaries
    nbasis = (p+1) * (p+1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics and number of Lebedev grid points
    call init_globals(0, p, vscales, ngrid, w, grid, vgrid)
    ! Compute radiuses
    sph2_r = sph1_r + dnrm2(3, sph2-sph1, 1)
    sph3_r = sph1_r + dnrm2(3, sph3-sph1, 1)
    ! Init values of spherical harmonics
    sph1_l = 1
    ! Get coefficients by baseline L2L
    sph2_l = 0
    call fmm_l2l_baseline(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_l, &
        & sph2_l)
    sph3_l = 0
    call fmm_l2l_baseline(sph1-sph3, sph1_r, sph3_r, p, vscales, sph1_l, &
        & sph3_l)
    sph4_l = 0
    call fmm_l2l_baseline(sph1-sph4, sph1_r, sph4_r, p, vscales, sph1_l, &
        & sph4_l)
    ! Check fmm_l2l_fast
    write(*, *)
    write(*, "(A,I0)") "Check fmm_l2l_fast vs fmm_l2l_baseline for p=", p
    write(*, "(A)") "======================================================="
    v = 0
    call fmm_l2l_fast(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_l, &
        & v)
    write(*,*) "|| L(1)2L(2) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph2_l, 1) / dnrm2(nbasis, sph2_l, 1)
    v = 0
    call fmm_l2l_fast(sph1-sph3, sph1_r, sph3_r, p, vscales, sph1_l, &
        & v)
    write(*,*) "|| L(1)2L(3) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph3_l, 1) / dnrm2(nbasis, sph3_l, 1)
    v = 0
    call fmm_l2l_fast(sph1-sph4, sph1_r, sph4_r, p, vscales, sph1_l, &
        & v)
    write(*,*) "|| L(1)2L(4) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph4_l, 1) / dnrm2(nbasis, sph4_l, 1)
    ! Check fmm_l2l_test_mat
    write(*, *)
    write(*, "(A,I0)") "Check fmm_l2l_test_mat vs fmm_l2l_baseline for p=", p
    write(*, "(A)") "======================================================="
    v = 0
    call fmm_l2l_test_mat(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_l, &
        & v)
    write(*,*) "|| L(1)2L(2) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph2_l, 1) / dnrm2(nbasis, sph2_l, 1)
    v = 0
    call fmm_l2l_test_mat(sph1-sph3, sph1_r, sph3_r, p, vscales, sph1_l, &
        & v)
    write(*,*) "|| L(1)2L(3) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph3_l, 1) / dnrm2(nbasis, sph3_l, 1)
    v = 0
    call fmm_l2l_test_mat(sph1-sph4, sph1_r, sph4_r, p, vscales, sph1_l, &
        & v)
    write(*,*) "|| L(1)2L(4) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph4_l, 1) / dnrm2(nbasis, sph4_l, 1)
    ! Check fmm_l2l_get_mat + fmm_l2l_use_mat
    write(*, *)
    write(*, "(A,A,I0)") "Check fmm_l2l_get_mat+fmm_l2l_use_mat vs ", &
        & "fmm_l2l_baseline for p=", p
    write(*, "(A)") "======================================================="
    v = 0
    call fmm_l2l_get_mat(sph1-sph2, sph1_r, sph2_r, p, vscales, reflect_mat, &
        & ztrans_mat)
    call fmm_l2l_use_mat(sph1-sph2, sph1_r, sph2_r, p, reflect_mat, &
        & ztrans_mat, sph1_l, v)
    write(*,*) "|| L(1)2L(2) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph2_l, 1) / dnrm2(nbasis, sph2_l, 1)
    v = 0
    call fmm_l2l_get_mat(sph1-sph3, sph1_r, sph3_r, p, vscales, reflect_mat, &
        & ztrans_mat)
    call fmm_l2l_use_mat(sph1-sph3, sph1_r, sph3_r, p, reflect_mat, &
        & ztrans_mat, sph1_l, v)
    write(*,*) "|| L(1)2L(3) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph3_l, 1) / dnrm2(nbasis, sph3_l, 1)
    v = 0
    call fmm_l2l_get_mat(sph1-sph4, sph1_r, sph4_r, p, vscales, reflect_mat, &
        & ztrans_mat)
    call fmm_l2l_use_mat(sph1-sph4, sph1_r, sph4_r, p, reflect_mat, &
        & ztrans_mat, sph1_l, v)
    write(*,*) "|| L(1)2L(4) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasis, v-sph4_l, 1) / dnrm2(nbasis, sph4_l, 1)
    ! Now check performance
    write(*, *)
    write(*, "(A,I0)") "Check performance of L2L for p=", p
    write(*, "(A)") "======================================================="
    call cpu_time(start)
    do i = 1, 1000
        call fmm_l2l_baseline(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_l, &
            & sph2_l)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_l2l_baseline:", (finish-start)/1000, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 1000
        call fmm_l2l_fast(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_l, &
            & sph2_l)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_l2l_fast:", (finish-start)/1000, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 1000
        call fmm_l2l_test_mat(sph1-sph2, sph1_r, sph2_r, p, vscales, sph1_l, &
            & sph2_l)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_l2l_test_mat:", (finish-start)/1000, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 1000
        call fmm_l2l_get_mat(sph1-sph2, sph1_r, sph2_r, p, vscales, &
            & reflect_mat, ztrans_mat)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_l2l_get_mat:", (finish-start)/1000, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 1000
        call fmm_l2l_use_mat(sph1-sph2, sph1_r, sph2_r, p, reflect_mat, &
            & ztrans_mat, sph1_l, v)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_l2l_use_mat:", (finish-start)/1000, &
        & " seconds"
end subroutine check_l2l_improved

subroutine check_m2l_improved(pm, pl)
    use pcm_fmm
    implicit none
    integer, intent(in) :: pm, pl
    integer, parameter :: ngrid=1
    integer :: nbasism, nbasisl, i
    real(kind=8) :: vscales((pm+pl+1)*(pm+pl+1)), w(ngrid), grid(3, ngrid)
    real(kind=8) :: vgrid((pl+1)*(pl+1), ngrid)
    ! sph1 contains multipole expansion
    real(kind=8) :: sph1(3)=(/0.0, 0.0, 0.0/), sph1_r=2.0
    ! sph2 and sph3 contain local expansion
    real(kind=8) :: sph2(3)=(/-1.0, -1.0, -2.0/), sph2_r=3
    real(kind=8) :: sph3(3)=(/0.0, 0.0, -2.2/), sph3_r=1.6
    real(kind=8), external :: dnrm2
    real(kind=8) :: sph1_m((pm+1)*(pm+1))
    real(kind=8) :: sph2_l((pl+1)*(pl+1))
    real(kind=8) :: sph3_l((pl+1)*(pl+1))
    real(kind=8) :: v((pl+1)*(pl+1))
    real(kind=8) :: reflect_mat((max(pm,pl)+1) * (2*max(pm,pl)+1) &
        & * (2*max(pm,pl)+3) / 3)
    real(kind=8) :: ztrans_mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6)
    real(kind=8) :: start, finish
    ! Preliminaries
    nbasism = (pm+1) * (pm+1)
    nbasisl = (pl+1) * (pl+1)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics and number of Lebedev grid points
    call init_globals(pm, pl, vscales, ngrid, w, grid, vgrid)
    ! Init values of spherical harmonics
    sph1_m = 1
    ! Get coefficients by baseline M2L
    sph2_l = 0
    call fmm_m2l_baseline(sph1-sph2, sph1_r, sph2_r, pm, pl, vscales, sph1_m, &
        & sph2_l)
    sph3_l = 0
    call fmm_m2l_baseline(sph1-sph3, sph1_r, sph3_r, pm, pl, vscales, sph1_m, &
        & sph3_l)
    ! Check fmm_m2l_fast
    write(*, *)
    write(*, "(A,I0,A,I0)") "Check fmm_m2l_fast vs fmm_m2l_baseline for pm=", &
        & pm, " and pl=", pl
    write(*, "(A)") "======================================================="
    v = 0
    call fmm_m2l_fast(sph1-sph2, sph1_r, sph2_r, pm, pl, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2L(2) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasisl, v-sph2_l, 1) / dnrm2(nbasisl, sph2_l, 1)
    v = 0
    call fmm_m2l_fast(sph1-sph3, sph1_r, sph3_r, pm, pl, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2L(3) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasisl, v-sph3_l, 1) / dnrm2(nbasisl, sph3_l, 1)
    ! Check fmm_m2l_test_mat
    write(*, *)
    write(*, "(A,A,I0,A,I0)") "Check fmm_m2l_test_mat vs fmm_m2l_baseline ", &
        & "for pm=", pm, " and pl=", pl
    write(*, "(A)") "======================================================="
    v = 0
    call fmm_m2l_test_mat(sph1-sph2, sph1_r, sph2_r, pm, pl, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2L(2) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasisl, v-sph2_l, 1) / dnrm2(nbasisl, sph2_l, 1)
    v = 0
    call fmm_m2l_test_mat(sph1-sph3, sph1_r, sph3_r, pm, pl, vscales, sph1_m, &
        & v)
    write(*,*) "|| M(1)2L(3) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasisl, v-sph3_l, 1) / dnrm2(nbasisl, sph3_l, 1)
    ! Check fmm_m2l_get_mat + fmm_m2l_use_mat
    write(*, *)
    write(*, "(A,A,I0,A,I0)") "Check fmm_m2l_get_mat+fmm_m2l_use_mat vs ", &
        & "fmm_m2l_baseline for pm=", pm, " and pl=", pl
    write(*, "(A)") "======================================================="
    v = 0
    call fmm_m2l_get_mat(sph1-sph2, sph1_r, sph2_r, pm, pl, vscales, &
        & reflect_mat, ztrans_mat)
    call fmm_m2l_use_mat(sph1-sph2, sph1_r, sph2_r, pm, pl, reflect_mat, &
        & ztrans_mat, sph1_m, v)
    write(*,*) "|| M(1)2L(2) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasisl, v-sph2_l, 1) / dnrm2(nbasisl, sph2_l, 1)
    v = 0
    call fmm_m2l_get_mat(sph1-sph3, sph1_r, sph3_r, pm, pl, vscales, &
        & reflect_mat, ztrans_mat)
    call fmm_m2l_use_mat(sph1-sph3, sph1_r, sph3_r, pm, pl, reflect_mat, &
        & ztrans_mat, sph1_m, v)
    write(*,*) "|| M(1)2L(3) - baseline ||  /  || baseline || =", &
        & dnrm2(nbasisl, v-sph3_l, 1) / dnrm2(nbasisl, sph3_l, 1)
    ! Now check performance
    write(*, *)
    write(*, "(A,I0,A,I0)") "Check performance of M2L for pm=", pm, &
        & " and pl=", pl
    write(*, "(A)") "======================================================="
    call cpu_time(start)
    do i = 1, 1000
        call fmm_m2l_baseline(sph1-sph2, sph1_r, sph2_r, pm, pl, vscales, &
            & sph1_m, sph2_l)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_m2l_baseline:", (finish-start)/1000, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 1000
        call fmm_m2l_fast(sph1-sph2, sph1_r, sph2_r, pm, pl, vscales, &
            & sph1_m, sph2_l)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_m2l_fast:", (finish-start)/1000, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 1000
        call fmm_m2l_test_mat(sph1-sph2, sph1_r, sph2_r, pm, pl, vscales, &
            & sph1_m, sph2_l)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_m2l_test_mat:", (finish-start)/1000, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 1000
        call fmm_m2l_get_mat(sph1-sph2, sph1_r, sph2_r, pm, pl, vscales, &
            & reflect_mat, ztrans_mat)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_m2l_get_mat:", (finish-start)/1000, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 1000
        call fmm_m2l_use_mat(sph1-sph2, sph1_r, sph2_r, pm, pl, reflect_mat, &
            & ztrans_mat, sph1_m, v)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " fmm_m2l_use_mat:", (finish-start)/1000, &
        & " seconds"
end subroutine check_m2l_improved

