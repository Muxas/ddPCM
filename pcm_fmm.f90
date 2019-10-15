module pcm_fmm
implicit none
real(kind=8) :: sqrt_2=0, sqrt_four_pi=0
contains

! Init global constants
subroutine init_globals()
    sqrt_2 = sqrt(dble(2))
    sqrt_four_pi = 4*sqrt(atan(dble(1)))
end subroutine init_globals

! Compute cos(mx) and sin(mx) from a given cos(x) and sin(x)
subroutine trgev(c, s, p, vcos, vsin)
! Parameters:
!   c: cos(x)
!   s: sin(x)
!   p: maximum value of m, for which to compute cos(mx) and sin(mx)
!   vcos: values of cos(mx) for m from 0 to p
!   vsin: values of sin(mx) for m from 0 to p
    real(kind=8), intent(in) :: c, s
    integer, intent(in) :: p
    real(kind=8), intent(out) :: vcos(p+1), vsin(p+1)
    integer :: m
    vcos(1) = 1
    vsin(1) = 0
    do m = 2, p+1
        vcos(m) = vcos(m-1)*c - vsin(m-1)*s
        vsin(m) = vcos(m-1)*s + vsin(m-1)*c
    end do
end subroutine trgev

! Compute associated Legendre polynomials P(l,m)(x)
! Uses recurrence formula
!   (l-m)P(l,m) = x(2l-1)P(l-1,m) - (l+m-1)P(l-2,m)
! Input x must be in range [-1;1]
! This function is simply copied from ddCOSMO
subroutine polleg(x, y, p, vplm)
! Parameters:
!   x: coordinate in interval [-1;1]
!   y: sqrt(1-x*2)
!   p: maximum degree of polynomials to compute
!   vplm: values of associated legendre polynomials P(l,m)(x)
    real(kind=8), intent(in) :: x, y
    integer, intent(in) :: p
    real(kind=8), intent(out) :: vplm((p+1)*(p+1))
    integer :: m, ind, l, ind2
    real(kind=8) :: fact, pmm, pmm1, pmmo, pll, fm, fl
    fact = 1
    pmm = 1
    do m = 0, p
        ind = (m + 1)*(m + 1)
        vplm(ind) = pmm
        if (m .eq. p) then
            return
        end if
        fm = dble(m)
        pmm1 = x * (2*fm+1) * pmm
        ind2 = ind + 2*m + 2
        vplm(ind2) = pmm1
        pmmo = pmm
        do l = m+2, p
            fl = dble(l)
            pll = x*(2*fl-1)*pmm1 - (fl+fm-1)*pmm
            pll = pll / (fl-fm)
            ind = l*l + l + 1
            vplm(ind+m) = pll
            pmm = pmm1
            pmm1 = pll
        end do
        pmm = -pmmo * fact * y
        fact = fact + 2
    end do
end subroutine polleg

! Compute scaling factors for normalized real spherical harmonics
subroutine scales_real_normal(p, vscales)
! Parameters:
!   p: maximum degree of spherical harmonics
!   vscale: values of scaling factors
    integer, intent(in) :: p
    real(kind=8), intent(out) :: vscales((p+1)*(p+1))
    real(kind=8) :: tmp
    integer :: l, ind, m
    do l = 0, p
        ! m = 0
        ind = l*l + l + 1
        vscales(ind) = sqrt(dble(2*l+1)) / sqrt_four_pi
        tmp = vscales(ind) * sqrt_2
        ! m != 0
        do m = 1, l
            tmp = -tmp / sqrt(dble((l-m+1)*(l+m)))
            vscales(ind+m) = tmp
        end do
    end do
end subroutine scales_real_normal

! Compute multipole coefficients for particle of unit charge
! Based on normalized scaled real spherical harmonics of given radius
! This function is not needed for pcm, but it is useful for testing purposes
subroutine fmm_p2m(c, r, p, m)
! Parameters:
!   c: coordinates of charged particle (relative to center of harmonics)
!   r: radius of spherical harmonics
!   p: maximum degree of multipole basis functions
!   m: multipole coefficients
    real(kind=8), intent(in) :: c(3), r
    integer, intent(in) :: p
    real(kind=8), intent(out) :: m((p+1)*(p+1))
    real(kind=8) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=8) :: vplm((p+1)*(p+1)), t, vscales((p+1)*(p+1)), tmp, rcoef
    integer :: n, k, ind
    stheta = c(1)*c(1) + c(2)*c(2)
    rho = sqrt(c(3)*c(3) + stheta)
    if (rho .ne. 0) then
        ctheta = c(3) / rho
        if (stheta .ne. 0) then
            stheta = sqrt(stheta)
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / rho
            call trgev(cphi, sphi, p, vcos, vsin)
        else
            cphi = 1
            sphi = 0
            vcos = 1
            vsin = 0
        end if
        call polleg(ctheta, stheta, p, vplm)
        call scales_real_normal(p, vscales)
        ! Now build harmonics to fill multipole coefficients
        rcoef = rho / r
        t = 1
        do n = 0, p
            ind = n*n + n + 1
            m(ind) = t * vscales(ind) * vplm(ind)
            do k = 1, n
                tmp = t * vscales(ind+k) * vplm(ind+k)
                m(ind+k) = tmp * vcos(k+1)
                m(ind-k) = -(tmp * vsin(k+1))
            end do
            t = t * rcoef
        end do
    else
        m(1) = 1
        m(2:) = 0
    end if
end subroutine fmm_p2m

! Compute potential, induced by spherical harmonics
! Based on normalized scaled real spherical harmonics of given radius
subroutine fmm_m2p(c, r, p, m, v)
! Parameters:
!   c: relative distance from center of harmonics to point of potential
!   r: radius of spherical harmonics
!   p: maximum degree of multipole basis functions
!   m: multipole expansion at origin
!   v: value of induced potential
    real(kind=8), intent(in) :: c(3), r
    real(kind=8), intent(in) :: m((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: v
    real(kind=8) :: tmp, tmp1, tmp2
    real(kind=8) :: rho, t, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=8) :: vscales((p+1)*(p+1)), vplm((p+1)*(p+1)), rcoef
    integer :: n, k, ind
    t = 1 / r
    stheta = c(1)*c(1) + c(2)*c(2)
    rho = sqrt(c(3)*c(3) + stheta)
    v = 0
    if (rho .eq. 0) then
        return
    end if
    ! rho is always > 0
    ctheta = c(3) / rho
    if (stheta .ne. 0) then
        stheta = sqrt(stheta)
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        stheta = stheta / rho
        call trgev(cphi, sphi, p, vcos, vsin)
    else
        cphi = 1
        sphi = 0
        vcos = 1
        vsin = 0
    end if
    call polleg(ctheta, stheta, p, vplm)
    call scales_real_normal(p, vscales)
    rcoef = r / rho
    do n = 0, p
        t = t * rcoef
        ind = n*n + n + 1
        ! k = 0
        tmp = m(ind) * vplm(ind) / vscales(ind)
        ! k != 0
        do k = 1, n
            tmp1 = vplm(ind+k) * vscales(ind+k) / vscales(ind)**2
            tmp2 = m(ind+k)*vcos(k+1) - m(ind-k)*vsin(k+1)
            tmp = tmp + tmp1*tmp2
        end do
        v = v + t*tmp
    end do
end subroutine fmm_m2p

! M2M baseline translation (p^4 operations)
! Baseline in terms of operation count: p^4
subroutine fmm_m2m_baseline(c, src_r, dst_r, p, src_m, dst_m)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r, src_m((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(inout) :: dst_m((p+1)*(p+1))
    real(kind=8) :: r, r1, r2, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=8) :: vplm((p+1)*(p+1)), vscales((p+1)*(p+1))
    real(kind=8) :: fact(2*p+1), tmpk1, tmpk2, tmpk3, tmp1, tmp2
    real(kind=8) :: pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, m, indj, indm, indn, indjn
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    if (r .ne. 0) then
        ctheta = c(3) / r
        if (stheta .ne. 0) then
            stheta = sqrt(stheta)
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / r
            call trgev(cphi, sphi, p, vcos, vsin)
        else
            cphi = 1
            sphi = 0
            vcos = 1
            vsin = 0
        end if
        call polleg(ctheta, stheta, p, vplm)
        call scales_real_normal(p, vscales)
        r1 = src_r / dst_r
        r2 = r / dst_r
        pow_r1(1) = 1
        pow_r2(1) = 1
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        ! Fill square roots of factorials
        fact(1) = 1
        do j = 2, 2*p+1
            fact(j) = sqrt(dble(j-1)) * fact(j-1)
        end do
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = fact(j-k+1) * fact(j+k+1)
                if (k .ne. 0) then
                    tmp1 = tmp1 * sqrt_2
                end if
                do n = 0, j
                    indn = n*n + n + 1
                    indjn = (j-n)**2 + (j-n) + 1
                    tmp2 = tmp1 * pow_r1(j-n+1) * pow_r2(n+1) * &
                        & vscales(indj) / vscales(indjn) / vscales(indn)
                    do m = max(k+n-j, -n), min(k+j-n, n)
                        indm = indn + abs(m)
                        cphi = vcos(1+abs(m))
                        sphi = vsin(1+abs(m))
                        tmpk1 = tmp2 / fact(n-m+1) / fact(n+m+1) / &
                            & fact(j-n-k+m+1) / fact(j-n+k-m+1) * &
                            & vplm(indm) * vscales(indm)
                        if (mod(abs(k-abs(m)-abs(k-m)), 4) .eq. 2) then
                            tmpk1 = -tmpk1
                        end if
                        tmpk2 = src_m(indjn+abs(k-m)) * cphi
                        if ((m .lt. 0) .or. (m .gt. k)) then
                            sphi = -sphi
                        end if
                        tmpk3 = -src_m(indjn+abs(k-m)) * sphi
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt_2
                            tmpk2 = tmpk2 + src_m(indjn-abs(k-m))*sphi
                            tmpk3 = tmpk3 + src_m(indjn-abs(k-m))*cphi
                        end if
                        if (m .gt. k) then
                            tmpk3 = -tmpk3
                        end if
                        dst_m(indj+k) = dst_m(indj+k) + tmpk1*tmpk2
                        if (k .ne. 0) then
                            dst_m(indj-k) = dst_m(indj-k) + tmpk1*tmpk3
                        end if
                    end do
                end do
            end do
        end do
    else
        r1 = src_r / dst_r
        tmpk1 = 1
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_m(k) = dst_m(k) + src_m(k)*tmpk1
            end do
            tmpk1 = tmpk1 * r1
        end do
    end if
end subroutine fmm_m2m_baseline

! Divide given cluster of spheres into two subclusters by inertial bisection
subroutine cluster_divide(nsph, csph, n, ind, div)
! Parameters:
!   nsph: Number of all spheres
!   csph: Centers of all spheres
!   n: Number of spheres in given cluster
!   ind: Indexes of spheres in given cluster (sorted on exit)
!   div: Break point between two clusters. ind(1:div) belong to first cluster
!       and ind(div+1:n) belongs to second cluster
    integer, intent(in) :: nsph, n
    real(kind=8), intent(in) :: csph(3, nsph)
    integer, intent(inout) :: ind(n)
    integer, intent(out) :: div
    real(kind=8) :: c(3), tmp_csph(3, n), a(3, 3), w(3), work(9), scal(n)
    real(kind=8) :: alpha=1, beta=0
    integer :: i, l, r, lwork=9, info, tmp_ind(n)
    c = 0
    do i = 1, n
        c = c + csph(:, ind(i))
    end do
    c = c / n
    do i = 1, n
        tmp_csph(:, i) = csph(:, ind(i)) - c
    end do
    call dgemm('N', 'T', 3, 3, n, alpha, tmp_csph, 3, tmp_csph, 3, beta, a, 3)
    call dsyev('V', 'L', 3, a, 3, w, work, lwork, info)
    call dgemv('T', 3, n, alpha, tmp_csph, 3, a(:, 3), 1, beta, scal, 1)
    l = 1
    r = n
    do i = 1, n
        if (scal(i) .ge. 0) then
            tmp_ind(l) = ind(i)
            l = l + 1
        else
            tmp_ind(r) = ind(i)
            r = r - 1
        end if
    end do
    div = r
    ind = tmp_ind
end subroutine cluster_divide

! Prepare tree (divide and compute bounding spheres)
! Number of clusters is always 2*nsph-1
subroutine tree_init(nsph, csph, rsph, ind, cluster, children, parent, cnode, &
        & rnode)
! Parameters:
!   nsph: Number of all spheres
!   csph: Centers of all spheres
!   rsph: Radiuses of all spheres
!   ind: Permutation of spheres (to localize them)
!   cluster: first and last spheres (from ind array), belonging to each cluster
!   children: children of each cluster. 0 means no children
!   parent: parent of each cluster. 0 means no parent
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
    integer, intent(in) :: nsph
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph)
    integer, intent(inout) :: ind(nsph)
    integer, intent(out) :: cluster(2, 2*nsph-1), children(2, 2*nsph-1)
    integer, intent(out) :: parent(2*nsph-1)
    real(kind=8), intent(out) :: cnode(3, 2*nsph-1), rnode(2*nsph-1)
    integer :: i, j, n, s, e, div
    real(kind=8) :: r, r1, r2, c(3), c1(3), c2(3), d
    cluster(1, 1) = 1
    cluster(2, 1) = nsph
    parent(1) = 0
    j = 2
    ! Divide tree until leaves with single spheres inside
    do i = 1, 2*nsph-1
        s = cluster(1, i)
        e = cluster(2, i)
        n = e - s + 1
        if (n .gt. 1) then
            call cluster_divide(nsph, csph, n, ind(s:e), div)
            cluster(1, j) = s
            cluster(2, j) = s + div - 1
            cluster(1, j+1) = s + div
            cluster(2, j+1) = e
            children(1, i) = j
            children(2, i) = j + 1
            parent(j) = i
            parent(j+1) = i
            j = j + 2
        else
            children(:, i) = 0
        end if
    end do
    ! Compute bounding spheres
    do i = 2*nsph-1, 1, -1
        if (children(1, i) .eq. 0) then
            j = cluster(1, i)
            cnode(:, i) = csph(:, ind(j))
            rnode(i) = rsph(ind(j))
        else
            j = children(1, i)
            c1 = cnode(:, j)
            r1 = rnode(j)
            j = children(2, i)
            c2 = cnode(:, j)
            r2 = rnode(j)
            c = c1 - c2
            d = sqrt(c(1)*c(1) + c(2)*c(2) + c(3)*c(3))
            if ((r1-r2) .ge. d) then
                c = c1
                r = r1
            else if((r2-r1) .ge. d) then
                c = c2
                r = r2
            else
                r = (r1+r2+d) / 2
                c = c2 + c/d*(r-r2)
            end if
            cnode(:, i) = c
            rnode(i) = r
        end if
    end do
end subroutine tree_init

! Compute multipole coefficients for each node of tree
subroutine tree_m2m_baseline(nsph, p, coef_sph, ind, cluster, children, &
        & cnode, rnode, coef_node)
! Parameters:
!   nsph: Number of all spheres
!   p: maximum degree of multipole basis functions
!   coef_sph: multipole coefficients of input spheres
!   ind: permutation of all spheres (to localize sequential spheres)
!   cluster: first and last spheres (from ind array), belonging to each cluster
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   coef_node: multipole coefficients of bounding spheres of nodes
    integer, intent(in) :: nsph, p, ind(nsph), cluster(2, 2*nsph-1)
    real(kind=8), intent(in) :: coef_sph((p+1)*(p+1), nsph)
    integer, intent(in) :: children(2, 2*nsph-1)
    real(kind=8), intent(in) :: cnode(3, 2*nsph-1), rnode(2*nsph-1)
    real(kind=8), intent(out) :: coef_node((p+1)*(p+1), 2*nsph-1)
    integer :: i, j(2)
    real(kind=8) :: c1(3), c2(3), c(3), r1, r2, r
    do i = 2*nsph-1, 1, -1
        j = children(:, i)
        if (j(1) .eq. 0) then
            coef_node(:, i) = coef_sph(:, ind(cluster(1, i)))
        else
            c = cnode(:, i)
            r = rnode(i)
            c1 = cnode(:, j(1))
            r1 = rnode(j(1))
            c2 = cnode(:, j(2))
            r2 = rnode(j(2))
            coef_node(:, i) = 0
            call fmm_m2m_baseline(c1-c, r1, r, p, coef_node(:, j(1)), &
                & coef_node(:, i))
            call fmm_m2m_baseline(c2-c, r2, r, p, coef_node(:, j(2)), &
                & coef_node(:, i))
        end if
    end do
end subroutine tree_m2m_baseline

end module

