program test_pcm_fmm
!use ddcosmo
implicit none
integer :: p1=6, p2=10, ngrid1=6, ngrid2=110, nsph1=1000

call check_p2m_m2p_m2m_baseline(p1)
call check_p2m_m2p_m2m_baseline(p2)
call check_m2l_l2p_l2l_baseline(p1, p1)
call check_m2l_l2p_l2l_baseline(p1, p2)
call check_m2l_l2p_l2l_baseline(p2, p1)
call check_m2m_improved(p1)
call check_m2m_improved(p2)
call check_l2l_improved(p1)
call check_l2l_improved(p2)
call check_m2l_improved(p1, p1)
call check_m2l_improved(p1, p2)
call check_m2l_improved(p2, p1)
call check_tree_m2m_l2l(nsph1, p1)
call check_tree_m2m_l2l(nsph1, p2)
call check_tree_m2l(nsph1, p1, p1)
call check_tree_m2l(nsph1, p1, p2)
call check_tree_m2l(nsph1, p2, p1)

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

subroutine check_tree_m2m_l2l(nsph, p)
    use pcm_fmm
    implicit none
    integer, intent(in) :: nsph, p
    integer, parameter :: ngrid=1
    real(kind=8) :: csph(3, nsph), rsph(nsph), ui(ngrid, nsph)
    real(kind=8) :: coef_sph((p+1)*(p+1), nsph)
    real(kind=8) :: coef_sph_mat((p+1)*(p+1), nsph)
    real(kind=8) :: vscales((p+1)*(p+1)), w(ngrid), grid(3, ngrid)
    real(kind=8) :: vgrid((p+1)*(p+1), ngrid)
    integer :: iseed(4), ind(nsph), i, nclusters, height
    integer :: cluster(2, 2*nsph-1), children(2, 2*nsph-1)
    integer :: parent(2*nsph-1)
    real(kind=8) :: cnode(3, 2*nsph-1), rnode(2*nsph-1), start, finish
    integer :: snode(nsph), nnnear, nnfar
    integer, allocatable :: nfar(:), nnear(:)
    integer, allocatable :: far(:), near(:)
    integer, allocatable :: sfar(:), snear(:)
    integer :: lwork, iwork, jwork
    integer, allocatable :: work(:, :)
    real(kind=8), allocatable :: reflect_mat(:, :)
    real(kind=8), allocatable :: ztrans_mat(:, :)
    real(kind=8), allocatable :: coef_node(:, :), coef_node_mat(:, :)
    external :: dlarnv
    real(kind=8), external :: dnrm2
    ! Preliminaries
    iseed = (/0, 0, 0, 1/)
    call dlarnv(3, iseed, nsph*3, csph)
    rsph = 0.1
    coef_sph = 1
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics and number of Lebedev grid points
    call init_globals(0, p, vscales, ngrid, w, grid, vgrid)
    ! Build a tree
    do i = 1, nsph
        ind(i) = i
    end do
    call btree_init(nsph, csph, rsph, ind, cluster, children, parent, &
        cnode, rnode, snode)
    nclusters = 2*nsph-1
    ! Get its height
    call tree_get_height(nclusters, parent, height)
    ! Allocate arrays of coefficients of spherical harmonics
    allocate(coef_node((p+1)*(p+1), nclusters))
    allocate(coef_node_mat((p+1)*(p+1), nclusters))
    ! Allocate reflection and ztranslation matrices
    allocate(reflect_mat((p+1)*(2*p+1)*(2*p+3)/3, nclusters-1))
    allocate(ztrans_mat((p+1)*(p+2)*(2*p+3)/6, nclusters-1))
    ! Compute M2M for entire tree by fast code
    coef_sph = 1
    call tree_m2m_fast(nsph, nclusters, p, vscales, coef_sph, ind, cluster, &
        & children, cnode, rnode, coef_node)
    ! Compute M2M for entire tree by get_mat+use_mat code
    call tree_m2m_get_mat(nclusters, children, cnode, rnode, p, vscales, &
        & reflect_mat, ztrans_mat)
    call tree_m2m_use_mat(nsph, nclusters, p, coef_sph, ind, cluster, &
        & children, cnode, rnode, reflect_mat, ztrans_mat, coef_node_mat)
    write(*, *)
    write(*, "(A,A,I0)") "Check tree_m2m_get_mat+tree_m2m_use_mat vs ", &
        & "tree_m2m_fast for p=", p
    write(*, "(A)") "==================================="
    write(*, "(A,ES24.16E3)") "relative error is", dnrm2(nclusters*(p+1)*(p+1), &
        & coef_node_mat-coef_node, 1) / dnrm2(nclusters*(p+1)*(p+1), &
        & coef_node, 1)
    ! Check performance of M2M
    write(*, *)
    write(*, "(A,I0)") "Check performance of tree M2M for p=", p
    write(*, "(A)") "======================================="
    write(*, "(A,I0,A,I0,A)") "Tree has ", nclusters, " nodes and ", height, &
        & " levels"
    call cpu_time(start)
    do i = 1, 100
        call tree_m2m_fast(nsph, nclusters, p, vscales, coef_sph, ind, &
            & cluster, children, cnode, rnode, coef_node)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " tree_m2m_fast:", (finish-start)/100, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 100
        call tree_m2m_get_mat(nclusters, children, cnode, rnode, p, &
            & vscales, reflect_mat, ztrans_mat)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " tree_m2m_get_mat:", (finish-start)/100, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 100
        call tree_m2m_use_mat(nsph, nclusters, p, coef_sph, ind, cluster, &
            & children, cnode, rnode, reflect_mat, ztrans_mat, &
            & coef_node_mat)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " tree_m2m_use_mat:", (finish-start)/100, &
        & " seconds"
    ! Compute L2L for entire tree by fast code
    coef_node = 1
    call tree_l2l_fast(nsph, nclusters, p, vscales, coef_node, ind, cluster, &
        & children, cnode, rnode, coef_sph)
    ! Compute M2M for entire tree by get_mat+use_mat code
    coef_node_mat = 1
    call tree_l2l_get_mat(nclusters, children, cnode, rnode, p, vscales, &
        & reflect_mat, ztrans_mat)
    call tree_l2l_use_mat(nsph, nclusters, p, coef_node_mat, ind, cluster, &
        & children, cnode, rnode, reflect_mat, ztrans_mat, coef_sph_mat)
    write(*, *)
    write(*, "(A,A,I0)") "Check tree_l2l_get_mat+tree_l2l_use_mat vs ", &
        & "tree_l2l_fast for p=", p
    write(*, "(A)") "==================================="
    write(*, "(A,ES24.16E3)") "relative error is", &
        & dnrm2(nclusters*(p+1)*(p+1), coef_node_mat-coef_node, 1) &
        & / dnrm2(nclusters*(p+1)*(p+1), coef_node, 1)
    ! Check performance of L2L
    write(*, *)
    write(*, "(A,I0)") "Check performance of tree L2L for p=", p
    write(*, "(A)") "======================================="
    write(*, "(A,I0,A,I0,A)") "Tree has ", nclusters, " nodes and ", height, &
        & " levels"
    call cpu_time(start)
    do i = 1, 100
        call tree_l2l_fast(nsph, nclusters, p, vscales, coef_node, ind, &
            & cluster, children, cnode, rnode, coef_sph)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " tree_l2l_fast:", (finish-start)/100, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 100
        call tree_l2l_get_mat(nclusters, children, cnode, rnode, p, &
            & vscales, reflect_mat, ztrans_mat)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " tree_l2l_get_mat:", (finish-start)/100, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 100
        call tree_l2l_use_mat(nsph, nclusters, p, coef_node_mat, ind, &
            & cluster, children, cnode, rnode, reflect_mat, ztrans_mat, &
            & coef_sph_mat)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " tree_l2l_use_mat:", (finish-start)/100, &
        & " seconds"
    ! Deallocate all temporaries
    deallocate(ztrans_mat)
    deallocate(reflect_mat)
    deallocate(coef_node_mat)
    deallocate(coef_node)
end subroutine check_tree_m2m_l2l

subroutine check_tree_m2l(nsph, pm, pl)
    use pcm_fmm
    implicit none
    integer, intent(in) :: nsph, pm, pl
    integer, parameter :: ngrid=1
    real(kind=8) :: csph(3, nsph), rsph(nsph), ui(ngrid, nsph)
    real(kind=8) :: vscales((pm+pl+1)*(pm+pl+1)), w(ngrid), grid(3, ngrid)
    real(kind=8) :: vgrid((pl+1)*(pl+1), ngrid)
    integer :: iseed(4), ind(nsph), i, nclusters, height
    integer :: cluster(2, 2*nsph-1), children(2, 2*nsph-1)
    integer :: parent(2*nsph-1)
    real(kind=8) :: cnode(3, 2*nsph-1), rnode(2*nsph-1), start, finish
    integer :: snode(nsph), nnnear, nnfar
    integer, allocatable :: nfar(:), nnear(:)
    integer, allocatable :: far(:), near(:)
    integer, allocatable :: sfar(:), snear(:)
    integer :: lwork, iwork, jwork
    integer, allocatable :: work(:, :)
    real(kind=8), allocatable :: reflect_mat(:, :)
    real(kind=8), allocatable :: ztrans_mat(:, :)
    real(kind=8), allocatable :: coef_node_m(:, :), coef_node_l(:, :)
    real(kind=8), allocatable :: coef_node_l_mat(:, :)
    external :: dlarnv
    real(kind=8), external :: dnrm2
    ! Preliminaries
    iseed = (/0, 0, 0, 1/)
    call dlarnv(3, iseed, nsph*3, csph)
    rsph = 0.1
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics and number of Lebedev grid points
    call init_globals(pm, pl, vscales, ngrid, w, grid, vgrid)
    ! Build a tree
    do i = 1, nsph
        ind(i) = i
    end do
    call btree_init(nsph, csph, rsph, ind, cluster, children, parent, &
        cnode, rnode, snode)
    nclusters = 2*nsph-1
    ! Get its height
    call tree_get_height(nclusters, parent, height)
    ! Allocate arrays of coefficients of spherical harmonics
    allocate(coef_node_m((pm+1)*(pm+1), nclusters))
    allocate(coef_node_l((pl+1)*(pl+1), nclusters))
    allocate(coef_node_l_mat((pl+1)*(pl+1), nclusters))
    coef_node_m = 1
    ! Allocate space to find all admissible pairs
    allocate(nfar(nclusters), nnear(nclusters))
    ! Try to find all admissibly far and near pairs of tree nodes
    lwork = nclusters*200 ! Some magic constant which might be changed
    allocate(work(3, lwork))
    iwork = 0 ! init with zero for first call to tree_get_farnear_work
    call tree_get_farnear_work(nclusters, children, cnode, rnode, lwork, &
        iwork, jwork, work, nnfar, nfar, nnnear, nnear)
    ! Increase size of work array if needed and run again. Function
    ! tree_get_farnear_work uses previously computed work array, so it will not
    ! do the same work several times.
    if (iwork .ne. jwork+1) then
        write(*,*) 'Please increase lwork'
        stop
    end if
    ! Allocate arrays for admissible far and near pairs
    allocate(far(nnfar), near(nnnear), sfar(nclusters+1), snear(nclusters+1))
    ! Get list of admissible pairs from temporary work array
    ! This is needed only by M2L
    call tree_get_farnear(jwork, lwork, work, nclusters, nnfar, nfar, sfar, &
        far, nnnear, nnear, snear, near)
    ! Print some info about tree
    ! Allocate reflection and ztranslation matrices
    allocate(reflect_mat((max(pm,pl)+1)*(2*max(pm,pl)+1)*(2*max(pm,pl)+3)/3, &
        & nnfar))
    allocate(ztrans_mat((min(pm,pl)+1)*(min(pm,pl)+2) &
        & *(3*max(pm,pl)+3-min(pm,pl))/6, nnfar))
    ! Compute M2L for entire tree by fast code
    call tree_m2l_fast(nclusters, cnode, rnode, nnfar, sfar, far, pm, pl, &
        & vscales, coef_node_m, coef_node_l)
    ! Compute M2L for entire tree by get_mat+use_mat code
    call tree_m2l_get_mat(nclusters, cnode, rnode, nnfar, sfar, far, pm, &
        & pl, vscales, reflect_mat, ztrans_mat)
    call tree_m2l_use_mat(nclusters, cnode, rnode, nnfar, sfar, far, pm, &
        pl, reflect_mat, ztrans_mat, coef_node_m, coef_node_l_mat)
    write(*, *)
    write(*, "(A,A,I0,A,I0)") "Check tree_m2l_get_mat+tree_m2l_use_mat vs ", &
        & "tree_m2l_fast for pm=", pm, " and pl=", pl
    write(*, "(A)") "==================================="
    write(*, "(A,ES24.16E3)") "relative error is", &
        & dnrm2(nclusters*(pl+1)*(pl+1), coef_node_l_mat-coef_node_l, 1) &
        & / dnrm2(nclusters*(pl+1)*(pl+1), coef_node_l, 1)
    ! Check performance of M2L
    write(*, *)
    write(*, "(A,I0,A,I0)") "Check performance of tree M2L for pm=", pm, &
        & " and pl=", pl
    write(*, "(A)") "======================================="
    write(*, "(A,I0,A,I0,A,I0,A)") "Tree has ", nclusters, " nodes, ", &
        & height, " levels and ", nnfar, " far admissible pairs of nodes"
    call cpu_time(start)
    do i = 1, 10
        call tree_m2l_fast(nclusters, cnode, rnode, nnfar, sfar, far, pm, &
            & pl, vscales, coef_node_m, coef_node_l)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " tree_m2l_fast:", (finish-start)/10, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 10
        call tree_m2l_get_mat(nclusters, cnode, rnode, nnfar, sfar, far, &
            & pm, pl, vscales, reflect_mat, ztrans_mat)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " tree_m2l_get_mat:", (finish-start)/10, &
        & " seconds"
    call cpu_time(start)
    do i = 1, 10
        call tree_m2l_use_mat(nclusters, cnode, rnode, nnfar, sfar, far, &
            & pm, pl, reflect_mat, ztrans_mat, coef_node_m, coef_node_l_mat)
    end do
    call cpu_time(finish)
    write(*, "(A,ES9.2,A)") " tree_m2l_use_mat:", (finish-start)/10, &
        & " seconds"
    ! Deallocate all temporaries
    deallocate(ztrans_mat)
    deallocate(reflect_mat)
    deallocate(far, near, sfar, snear)
    deallocate(work)
    deallocate(nfar, nnear)
    deallocate(coef_node_l_mat)
    deallocate(coef_node_l)
    deallocate(coef_node_m)
end subroutine check_tree_m2l

