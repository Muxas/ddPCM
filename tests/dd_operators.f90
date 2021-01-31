!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/dd_operators.f90
!! Tests for dd_operators module
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-01-31

program test_dd_operators
use dd_operators
implicit none

integer :: i, iprint = 1
! Some implementation do not work for alpha=1d-307, so range of test values is
! reduced. This can be fixed by enforcing proper order of multiplications like
! a*b*c into a*(b*c).
real(dp) :: alpha(4)=(/1d0, -1d0, 1d-300, 1d+300/)

! Check tree M2M
do i = 1, size(alpha)
    call check_tree_m2m(0, alpha(i), iprint, 10d0*epsilon(one))
    call check_tree_m2m(1, alpha(i), iprint, 10d0*epsilon(one))
    call check_tree_m2m(10, alpha(i), iprint, 10d0*epsilon(one))
end do

! Check tree L2L
do i = 1, size(alpha)
    call check_tree_l2l(0, alpha(i), iprint, 10d0*epsilon(one))
    call check_tree_l2l(1, alpha(i), iprint, 10d0*epsilon(one))
    call check_tree_l2l(10, alpha(i), iprint, 10d0*epsilon(one))
end do

contains

subroutine check_tree_m2m(p, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: p, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
    integer, parameter :: nsph = 10
    integer :: i, j, k, order(nsph), istatus
    real(dp) :: csph(3, nsph), rsph(nsph), src(3, nsph), csph2(3, nsph), &
        & rsph2(nsph), sph_m((p+1)**2, nsph), node_m((p+1)**2, 2*nsph-1), &
        & node_m2((p+1)**2, 2*nsph-1), vscales((p+1)**2), vfact(2*p+1), &
        & rel_err, sph_m2((p+1)**2, nsph), full_norm, diff_norm, &
        & node_m3((p+1)**2, 2*nsph-1)
    type(dd_data_type) :: dd_data
    logical :: ok
    real(dp), external :: dnrm2
    real(dp) :: full_mat((p+1)**2, 2*nsph-1, (p+1)**2, nsph)
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
    rsph = abs(alpha) * (/2d-1, 3d-1, 4d-1, 5d-1, 6d-1, 7d-1, 8d-1, 9d-1, &
        & 1d0, 1.1d0/)
    src(:, 1) = alpha * (/1.1d0, 0.9d0, 1d0/)
    src(:, 2) = alpha * (/2.1d0, 1.9d0, 2.1d0/)
    src(:, 3) = alpha * (/1.1d0, 1.1d0, 3.1d0/)
    src(:, 4) = alpha * (/1.9d0, 1.9d0, 4.1d0/)
    src(:, 5) = alpha * (/1.1d0, 0.9d0, 5d0/)
    src(:, 6) = alpha * (/1.9d0, 2.1d0, 5.9d0/)
    src(:, 7) = alpha * (/1.1d0, 0.9d0, 7d0/)
    src(:, 8) = alpha * (/2.1d0, 2.1d0, 8d0/)
    src(:, 9) = alpha * (/1.1d0, 0.9d0, 9d0/)
    src(:, 10) = alpha * (/2.1d0, 1.9d0, 9.9d0/)
    ! Reorder inputs
    order = (/4, 2, 8, 5, 7, 1, 9, 6, 10, 3/)
    do i = 1, nsph
        csph2(:, i) = csph(:, order(i))
        rsph2(i) = rsph(order(i))
    end do
    dd_data % nclusters = 2*nsph-1
    dd_data % pm = p
    dd_data % nsph = nsph
    ! Allocate space for a tree
    allocate(dd_data % cluster(2, 2*nsph-1), dd_data % children(2, 2*nsph-1), &
        & dd_data % parent(2*nsph-1), dd_data % cnode(3, 2*nsph-1), &
        & dd_data % rnode(2*nsph-1), dd_data % snode(nsph), &
        & dd_data % order(nsph), dd_data % vscales((p+1)**2), &
        & dd_data % vfact(2*p+1))
    ! Build a recursive inertial binary tree
    call tree_rib_build(nsph, csph2, rsph2, dd_data % order, &
        & dd_data % cluster, dd_data % children, dd_data % parent, &
        & dd_data % cnode, dd_data % rnode, dd_data % snode)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(p, dd_data % vscales)
    dd_data % vfact(1) = one
    do i = 2, 2*p+1
        dd_data % vfact(i) = dd_data % vfact(i-1) * sqrt(dble(i-1))
    end do
    ! Init input harmonics
    do i = 1, nsph
        call fmm_p2m(src(:, i)-csph(:, i), one, rsph(i), p, &
            & dd_data % vscales, zero, sph_m(:, i))
    end do
    ! Get reference result of M2M operation
    do i = 1, 2*nsph-1
        node_m(:, i) = zero
        do j = dd_data % cluster(1, i), dd_data % cluster(2, i)
            k = dd_data % order(j)
            call fmm_m2m_rotation(csph2(:, k)-dd_data % cnode(:, i), &
                & rsph2(k), dd_data % rnode(i), p, dd_data % vscales, &
                & dd_data % vfact, one, sph_m(:, k), one, node_m(:, i))
        end do
    end do
    ! Check tree_m2m_rotation
    node_m2 = zero
    do i = 1, nsph
        node_m2(:, dd_data % snode(i)) = sph_m(:, i)
    end do
    call tree_m2m_rotation(dd_data, node_m2)
    rel_err = dnrm2(((p+1)**2)*(2*nsph-1), node_m-node_m2, 1) / &
        & dnrm2(((p+1)**2)*(2*nsph-1), node_m, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, *) ok, "tree_m2m_rotation rel_err=", rel_err
    end if
    if (.not. ok) stop 1
    ! Check tree_m2m_reflection
    node_m2 = zero
    do i = 1, nsph
        node_m2(:, dd_data % snode(i)) = sph_m(:, i)
    end do
    call tree_m2m_reflection(dd_data, node_m2)
    rel_err = dnrm2(((p+1)**2)*(2*nsph-1), node_m-node_m2, 1) / &
        & dnrm2(((p+1)**2)*(2*nsph-1), node_m, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, *) ok, "tree_m2m_reflection rel_err=", rel_err
    end if
    if (.not. ok) stop 1
    ! Allocate space for transfer matrices
    dd_data % m2m_reflect_mat_size = (p+1)*(2*p+1)*(2*p+3)/3
    dd_data % m2m_ztranslate_mat_size = (p+1)*(p+2)*(p+3)/6
    allocate(dd_data % m2m_reflect_mat(dd_data % m2m_reflect_mat_size, &
        & dd_data % nclusters-1), stat=istatus)
    if (istatus .ne. 0) stop "Allocation failed"
    allocate(dd_data % m2m_ztranslate_mat(dd_data % m2m_ztranslate_mat_size, &
        & dd_data % nclusters-1), stat=istatus)
    if (istatus .ne. 0) stop "Allocation failed"
    ! Compute transfer matrices
    call tree_m2m_reflection_get_mat(dd_data)
    ! Check tree_m2m_reflection_use_mat
    node_m2 = zero
    do i = 1, nsph
        node_m2(:, dd_data % snode(i)) = sph_m(:, i)
    end do
    call tree_m2m_reflection_use_mat(dd_data, node_m2)
    rel_err = dnrm2(((p+1)**2)*(2*nsph-1), node_m-node_m2, 1) / &
        & dnrm2(((p+1)**2)*(2*nsph-1), node_m, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, *) ok, "tree_m2m_reflection_use_mat rel_err=", rel_err
    end if
    if (.not. ok) stop 1
    ! Check tree_m2m_reflection_use_mat_adj
    do i = 1, nsph
        do j = 1, (p+1)**2
            node_m2 = zero
            node_m2(j, dd_data % snode(i)) = one
            call tree_m2m_reflection_use_mat(dd_data, node_m2)
            full_mat(:, :, j, i) = node_m2
        end do
    end do
    full_norm = dnrm2((2*nsph-1)*nsph*((p+1)**4), full_mat, 1)
    do i = 1, 2*nsph-1
        do j = 1, (p+1)**2
            node_m2 = zero
            node_m2(j, i) = one
            call tree_m2m_reflection_use_mat_adj(dd_data, node_m2)
            do k = 1, nsph
                full_mat(j, i, :, k) = full_mat(j, i, :, k) - &
                    & node_m2(:, dd_data % snode(k))
            end do
        end do
    end do
    diff_norm = dnrm2((2*nsph-1)*nsph*((p+1)**4), full_mat, 1)
    rel_err = diff_norm / full_norm
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, *) ok, "tree_m2m_reflection_use_mat_adj rel_err=", rel_err
    end if
    if (.not. ok) stop 1
    ! Check tree_m2m_rotation_adj
    node_m2 = node_m
    call tree_m2m_reflection_use_mat_adj(dd_data, node_m2)
    node_m3 = node_m
    call tree_m2m_rotation_adj(dd_data, node_m3)
    rel_err = dnrm2((2*nsph-1)*((p+1)**2), node_m2-node_m3, 1) / &
        dnrm2((2*nsph-1)*((p+1)**2), node_m2, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, *) ok, "tree_m2m_rotation_adj rel_err=", rel_err
    end if
    if (.not. ok) stop 1
    ! Check tree_m2m_rotation_adj
    node_m3 = node_m
    call tree_m2m_reflection_adj(dd_data, node_m3)
    rel_err = dnrm2((2*nsph-1)*((p+1)**2), node_m2-node_m3, 1) / &
        dnrm2((2*nsph-1)*((p+1)**2), node_m2, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, *) ok, "tree_m2m_reflection_adj rel_err=", rel_err
    end if
    if (.not. ok) stop 1
    ! Deallocate tree
    call ddfree(dd_data)
end subroutine check_tree_m2m

subroutine check_tree_l2l(p, alpha, iprint, threshold)
    ! Inputs
    integer, intent(in) :: p, iprint
    real(dp), intent(in) :: alpha, threshold
    ! Local variables
    integer, parameter :: nsph = 10
    integer :: i, j, k, order(nsph), istatus
    real(dp) :: csph(3, nsph), rsph(nsph), src(3, nsph), csph2(3, nsph), &
        & rsph2(nsph), sph_l((p+1)**2, nsph), node_l((p+1)**2, 2*nsph-1), &
        & node_l2((p+1)**2, 2*nsph-1), vscales((p+1)**2), vfact(2*p+1), &
        & rel_err, sph_l2((p+1)**2, nsph), full_norm, diff_norm, &
        & node_l3((p+1)**2, 2*nsph-1)
    type(dd_data_type) :: dd_data
    logical :: ok
    real(dp), external :: dnrm2
    real(dp) :: full_mat((p+1)**2, nsph, (p+1)**2, 2*nsph-1)
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
    rsph = abs(alpha) * (/2d-1, 3d-1, 4d-1, 5d-1, 6d-1, 7d-1, 8d-1, 9d-1, &
        & 1d0, 1.1d0/)
    src(:, 1) = alpha * (/1.1d0, 0.9d0, 1d0/)
    src(:, 2) = alpha * (/2.1d0, 1.9d0, 2.1d0/)
    src(:, 3) = alpha * (/1.1d0, 1.1d0, 3.1d0/)
    src(:, 4) = alpha * (/1.9d0, 1.9d0, 4.1d0/)
    src(:, 5) = alpha * (/1.1d0, 0.9d0, 5d0/)
    src(:, 6) = alpha * (/1.9d0, 2.1d0, 5.9d0/)
    src(:, 7) = alpha * (/1.1d0, 0.9d0, 7d0/)
    src(:, 8) = alpha * (/2.1d0, 2.1d0, 8d0/)
    src(:, 9) = alpha * (/1.1d0, 0.9d0, 9d0/)
    src(:, 10) = alpha * (/2.1d0, 1.9d0, 9.9d0/)
    ! Reorder inputs
    order = (/4, 2, 8, 5, 7, 1, 9, 6, 10, 3/)
    do i = 1, nsph
        csph2(:, i) = csph(:, order(i))
        rsph2(i) = rsph(order(i))
    end do
    dd_data % nclusters = 2*nsph-1
    dd_data % pl = p
    dd_data % pm = p
    dd_data % nsph = nsph
    ! Allocate space for a tree
    allocate(dd_data % cluster(2, 2*nsph-1), dd_data % children(2, 2*nsph-1), &
        & dd_data % parent(2*nsph-1), dd_data % cnode(3, 2*nsph-1), &
        & dd_data % rnode(2*nsph-1), dd_data % snode(nsph), &
        & dd_data % order(nsph), dd_data % vscales((p+1)**2), &
        & dd_data % vfact(2*p+1))
    ! Build a recursive inertial binary tree
    call tree_rib_build(nsph, csph2, rsph2, dd_data % order, &
        & dd_data % cluster, dd_data % children, dd_data % parent, &
        & dd_data % cnode, dd_data % rnode, dd_data % snode)
    ! Get constants, corresponding to given maximum degree of spherical
    ! harmonics
    call ylmscale(p, dd_data % vscales)
    dd_data % vfact(1) = one
    do i = 2, 2*p+1
        dd_data % vfact(i) = dd_data % vfact(i-1) * sqrt(dble(i-1))
    end do
    ! Init input harmonics
    node_l = zero
    do i = 1, nsph
        call fmm_p2m(src(:, i)-csph(:, i), one, rsph(i), p, &
            & dd_data % vscales, zero, node_l(:, dd_data % snode(i)))
    end do
    call tree_m2m_rotation(dd_data, node_l)
    ! Get reference result of L2L operation
    sph_l = zero
    do i = 1, 2*nsph-1
        do j = dd_data % cluster(1, i), dd_data % cluster(2, i)
            k = dd_data % order(j)
            call fmm_l2l_rotation(dd_data % cnode(:, i)-csph2(:, k), &
                & dd_data % rnode(i), rsph2(k), p, dd_data % vscales, &
                & dd_data % vfact, one, node_l(:, i), one, sph_l(:, k))
        end do
    end do
    ! Check tree_l2l_rotation
    node_l2 = node_l
    call tree_l2l_rotation(dd_data, node_l2)
    do i = 1, nsph
        sph_l2(:, i) = node_l2(:, dd_data % snode(i))
    end do
    rel_err = dnrm2(((p+1)**2)*nsph, sph_l-sph_l2, 1) / &
        & dnrm2(((p+1)**2)*nsph, sph_l, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, *) ok, "tree_l2l_rotation rel_err=", rel_err
    end if
    if (.not. ok) stop 1
    ! Check tree_l2l_reflection
    node_l2 = node_l
    call tree_l2l_reflection(dd_data, node_l2)
    do i = 1, nsph
        sph_l2(:, i) = node_l2(:, dd_data % snode(i))
    end do
    rel_err = dnrm2(((p+1)**2)*nsph, sph_l-sph_l2, 1) / &
        & dnrm2(((p+1)**2)*nsph, sph_l, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, *) ok, "tree_l2l_reflection rel_err=", rel_err
    end if
    if (.not. ok) stop 1
    ! Allocate space for transfer matrices
    dd_data % l2l_reflect_mat_size = (p+1)*(2*p+1)*(2*p+3)/3
    dd_data % l2l_ztranslate_mat_size = (p+1)*(p+2)*(p+3)/6
    allocate(dd_data % l2l_reflect_mat(dd_data % l2l_reflect_mat_size, &
        & dd_data % nclusters-1), stat=istatus)
    if (istatus .ne. 0) stop "Allocation failed"
    allocate(dd_data % l2l_ztranslate_mat(dd_data % l2l_ztranslate_mat_size, &
        & dd_data % nclusters-1), stat=istatus)
    if (istatus .ne. 0) stop "Allocation failed"
    ! Compute transfer matrices
    call tree_l2l_reflection_get_mat(dd_data)
    ! Check tree_m2m_reflection_use_mat
    node_l2 = node_l
    call tree_l2l_reflection_use_mat(dd_data, node_l2)
    do i = 1, nsph
        sph_l2(:, i) = node_l2(:, dd_data % snode(i))
    end do
    rel_err = dnrm2(((p+1)**2)*nsph, sph_l-sph_l2, 1) / &
        & dnrm2(((p+1)**2)*nsph, sph_l, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, *) ok, "tree_l2l_reflection_use_mat rel_err=", rel_err
    end if
    if (.not. ok) stop 1
    ! Check tree_l2l_reflection_use_mat_adj
    do i = 1, 2*nsph-1
        do j = 1, (p+1)**2
            node_l2 = zero
            node_l2(j, i) = one
            call tree_l2l_reflection_use_mat(dd_data, node_l2)
            do k = 1, nsph
                full_mat(:, k, j, i) = node_l2(:, dd_data % snode(k))
            end do
        end do
    end do
    full_norm = dnrm2((2*nsph-1)*nsph*((p+1)**4), full_mat, 1)
    do i = 1, nsph
        do j = 1, (p+1)**2
            node_l2 = zero
            node_l2(j, dd_data % snode(i)) = one
            call tree_l2l_reflection_use_mat_adj(dd_data, node_l2)
            full_mat(j, i, :, :) = full_mat(j, i, :, :) - node_l2
        end do
    end do
    diff_norm = dnrm2((2*nsph-1)*nsph*((p+1)**4), full_mat, 1)
    rel_err = diff_norm / full_norm
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, *) ok, "tree_l2l_reflection_use_mat_adj rel_err=", rel_err
    end if
    if (.not. ok) stop 1
    ! Check tree_l2l_rotation_adj
    node_l2 = zero
    do i = 1, nsph
        node_l2(:, dd_data % snode(i)) = node_l(:, dd_data % snode(i))
    end do
    call tree_l2l_reflection_use_mat_adj(dd_data, node_l2)
    node_l3 = zero
    do i = 1, nsph
        node_l3(:, dd_data % snode(i)) = node_l(:, dd_data % snode(i))
    end do
    call tree_l2l_rotation_adj(dd_data, node_l3)
    rel_err = dnrm2((2*nsph-1)*((p+1)**2), node_l2-node_l3, 1) / &
        dnrm2((2*nsph-1)*((p+1)**2), node_l2, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, *) ok, "tree_l2l_rotation_adj rel_err=", rel_err
    end if
    if (.not. ok) stop 1
    ! Check tree_l2l_reflection_adj
    node_l3 = zero
    do i = 1, nsph
        node_l3(:, dd_data % snode(i)) = node_l(:, dd_data % snode(i))
    end do
    call tree_l2l_reflection_adj(dd_data, node_l3)
    rel_err = dnrm2((2*nsph-1)*((p+1)**2), node_l2-node_l3, 1) / &
        dnrm2((2*nsph-1)*((p+1)**2), node_l2, 1)
    ok = rel_err .le. threshold
    if (iprint .gt. 0) then
        write(*, *) ok, "tree_l2l_reflection_adj rel_err=", rel_err
    end if
    if (.not. ok) stop 1
    ! Deallocate tree
    call ddfree(dd_data)
end subroutine check_tree_l2l

end program

