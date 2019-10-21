module pcm_fmm
implicit none
real(kind=8) :: sqrt_2=0, sqrt_four_pi=0
contains

! Init global constants
subroutine init_globals(p, vscales, ngrid, w, grid, vgrid)
! Parameters:
!   p: maximum degree of polynomials to compute
!   vscales: values of scaling factors for spherical harmonics
!   ngrid: number of Lebedev grid points
!   w: weights of Lebedev grid points
!   grid: coordinates of Lebedev grid points on unit sphere
!   vgrid: values of spherical harmonics at grid points
    integer, intent(in) :: p, ngrid
    real(kind=8), intent(out) :: vscales((p+1)*(p+1)), w(ngrid), grid(3, ngrid)
    real(kind=8), intent(out) :: vgrid((p+1)*(p+1), ngrid)
    integer :: i, n, m, indn, indm
    real(kind=8) :: c(3), ctheta, stheta, cphi, sphi, vplm((p+1)*(p+1))
    real(kind=8) :: vcos(p+1), vsin(p+1), tmp
    sqrt_2 = sqrt(dble(2))
    sqrt_four_pi = 4*sqrt(atan(dble(1)))
    call scales_real_normal(p, vscales)
    call llgrid(ngrid, w, grid)
    do i = 1, ngrid
        c = grid(:, i)
        ctheta = c(3)
        stheta = sqrt(c(1)*c(1) + c(2)*c(2))
        if (stheta .ne. 0) then
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            call trgev(cphi, sphi, p, vcos, vsin)
        else
            cphi = 1
            sphi = 0
            vcos = 1
            vsin = 0
        end if
        call polleg(ctheta, stheta, p, vplm)
        do n = 0, p
            indn = n*n + n + 1
            vgrid(indn, i) = vscales(indn) * vplm(indn)
            do m = 1, n
                indm = indn + m
                tmp = vscales(indm) * vplm(indm)
                vgrid(indn+m, i) = tmp * vcos(1+m)
                vgrid(indn-m, i) = tmp * vsin(1+m)
            end do
        end do
    end do
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
!   vscales: values of scaling factors
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
subroutine fmm_p2m(c, r, p, vscales, m)
! Parameters:
!   c: coordinates of charged particle (relative to center of harmonics)
!   r: radius of spherical harmonics
!   p: maximum degree of multipole basis functions
!   vscales: normalization constants for Y_lm
!   m: multipole coefficients
    real(kind=8), intent(in) :: c(3), r, vscales((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: m((p+1)*(p+1))
    real(kind=8) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=8) :: vplm((p+1)*(p+1)), t, tmp, rcoef
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
        ! Now build harmonics to fill multipole coefficients
        rcoef = rho / r
        t = 1
        do n = 0, p
            ind = n*n + n + 1
            m(ind) = t * vscales(ind) * vplm(ind)
            do k = 1, n
                tmp = t * vscales(ind+k) * vplm(ind+k)
                m(ind+k) = tmp * vcos(k+1)
                m(ind-k) = tmp * vsin(k+1)
            end do
            t = t * rcoef
        end do
    else
        m(1) = 1
        m(2:) = 0
    end if
end subroutine fmm_p2m

! Compute potential, induced by multipole spherical harmonics
! Based on normalized scaled real spherical harmonics of given radius
subroutine fmm_m2p(c, r, p, vscales, m, v)
! Parameters:
!   c: relative distance from center of harmonics to point of potential
!   r: radius of spherical harmonics
!   p: maximum degree of multipole basis functions
!   vscales: normalization constants for Y_lm
!   m: multipole expansion at origin
!   v: value of induced potential
    real(kind=8), intent(in) :: c(3), r, vscales((p+1)*(p+1)), m((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: v
    real(kind=8) :: tmp, tmp1, tmp2
    real(kind=8) :: rho, t, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=8) :: vplm((p+1)*(p+1)), rcoef
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
    rcoef = r / rho
    do n = 0, p
        t = t * rcoef
        ind = n*n + n + 1
        ! k = 0
        tmp = m(ind) * vplm(ind) * vscales(ind)
        ! k != 0
        do k = 1, n
            tmp1 = vplm(ind+k) * vscales(ind+k)
            tmp2 = m(ind+k)*vcos(k+1) + m(ind-k)*vsin(k+1)
            tmp = tmp + tmp1*tmp2
        end do
        v = v + t*tmp/vscales(ind)**2
    end do
end subroutine fmm_m2p

! M2M baseline translation (p^4 operations)
! Baseline in terms of operation count: p^4
subroutine fmm_m2m_baseline(c, src_r, dst_r, p, vscales, src_m, dst_m)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   vscales: normalization constants for Y_lm
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1))
    real(kind=8), intent(in) :: src_m((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(inout) :: dst_m((p+1)*(p+1))
    real(kind=8) :: r, r1, r2, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=8) :: vplm((p+1)*(p+1)), fact(2*p+1), tmpk1, tmpk2, tmpk3, tmp1
    real(kind=8) :: tmp2, pow_r1(p+1), pow_r2(p+1)
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
                tmp1 = vscales(indj) * fact(j-k+1) * fact(j+k+1)
                if (k .ne. 0) then
                    tmp1 = tmp1 * sqrt_2
                end if
                do n = 0, j
                    indn = n*n + n + 1
                    indjn = (j-n)**2 + (j-n) + 1
                    tmp2 = tmp1 * pow_r1(j-n+1) * pow_r2(n+1) / &
                        & vscales(indjn) / vscales(indn)
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
                        if ((m .ge. 0) .and. (m .le. k)) then
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

! Compute potential, induced by local spherical harmonics
! Based on normalized scaled real spherical harmonics of given radius
subroutine fmm_l2p(c, r, p, vscales, l, v)
! Parameters:
!   c: relative distance from center of harmonics to point of potential
!   r: radius of spherical harmonics
!   p: maximum degree of local basis functions
!   vscales: normalization constants for Y_lm
!   l: local expansion at origin
!   v: value of induced potential
    real(kind=8), intent(in) :: c(3), r, vscales((p+1)*(p+1)), l((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: v
    real(kind=8) :: tmp, tmp1, tmp2
    real(kind=8) :: rho, t, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=8) :: vplm((p+1)*(p+1)), rcoef
    integer :: n, k, ind
    t = 1 / r
    stheta = c(1)*c(1) + c(2)*c(2)
    rho = sqrt(c(3)*c(3) + stheta)
    if (rho .eq. 0) then
        ! Only first term
        v = t * l(1) / vscales(1)
        return
    end if
    v = 0
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
    rcoef = rho / r
    do n = 0, p
        ind = n*n + n + 1
        ! k = 0
        tmp = l(ind) * vplm(ind) * vscales(ind)
        ! k != 0
        do k = 1, n
            tmp1 = vplm(ind+k) * vscales(ind+k)
            tmp2 = l(ind+k)*vcos(k+1) + l(ind-k)*vsin(k+1)
            tmp = tmp + tmp1*tmp2
        end do
        v = v + t*tmp/vscales(ind)**2
        t = t * rcoef
    end do
end subroutine fmm_l2p

! Compute local expansion by given multipole expansion
subroutine fmm_m2l_baseline(c, src_r, dst_r, pm, pl, vscales, src_m, dst_l)
! Parameters:
!   c: radius-vector from new (local) to old (multipole) centers of harmonics
!   src_r: radius of old (multipole) harmonics
!   dst_r: radius of new (local) harmonics
!   pm: maximum degree of multipole spherical harmonics
!   pl: maximum degree of local spherical harmonics
!   vscales: normalization constants for Y_lm (of degree up to pl+pm)
!   src_m: expansion in old (multipole) harmonics
!   dst_l: expansion in new (local) harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    real(kind=8), intent(in) :: vscales((pm+pl+1)*(pm+pl+1))
    real(kind=8), intent(in) :: src_m((pm+1)*(pm+1))
    integer, intent(in) :: pm, pl
    real(kind=8), intent(inout) :: dst_l((pl+1)*(pl+1))
    real(kind=8) :: r, r1, r2, ctheta, stheta, cphi, sphi
    real(kind=8) :: vcos(pm+pl+1), vsin(pm+pl+1)
    real(kind=8) :: vplm((pm+pl+1)*(pm+pl+1)), fact(2*(pm+pl)+1), tmpk1, tmpk2
    real(kind=8) :: tmpk3, tmp1, tmp2, pow_r1(pm+1), pow_r2(pl+1)
    integer :: j, k, n, m, indj, indmk, indn, indjn
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    ! r cannot be zero, as input sphere (multipole) must not intersect with
    ! output sphere (local)
    if (r .eq. 0) then
        return
    end if
    ctheta = c(3) / r
    if (stheta .ne. 0) then
        stheta = sqrt(stheta)
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        stheta = stheta / r
        call trgev(cphi, sphi, pm+pl, vcos, vsin)
    else
        cphi = 1
        sphi = 0
        vcos = 1
        vsin = 0
    end if
    call polleg(ctheta, stheta, pm+pl, vplm)
    r1 = src_r / r
    r2 = dst_r / r
    pow_r1(1) = 1
    pow_r2(1) = r2
    do j = 2, pm+1
        pow_r1(j) = pow_r1(j-1) * r1
    end do
    do j = 2, pl+1
        pow_r2(j) = pow_r2(j-1) * r2
    end do
    ! Fill square roots of factorials
    fact(1) = 1
    do j = 2, 2*(pm+pl)+1
        fact(j) = sqrt(dble(j-1)) * fact(j-1)
    end do
    do j = 0, pl
        indj = j*j + j + 1
        do k = 0, j
            tmp1 = vscales(indj) * pow_r2(j+1) / fact(j-k+1) / fact(j+k+1)
            if (k .ne. 0) then
                tmp1 = tmp1 * sqrt_2
            end if
            do n = 0, pm
                indn = n*n + n + 1
                indjn = (j+n)**2 + (j+n) + 1
                tmp2 = tmp1 * pow_r1(n+1) / vscales(indjn) / vscales(indn)
                if (mod(n, 2) .eq. 1) then
                    tmp2 = -tmp2
                end if
                do m = -n, n
                    indmk = indjn + abs(m-k)
                    cphi = vcos(1+abs(m-k))
                    sphi = vsin(1+abs(m-k))
                    tmpk1 = tmp2 / fact(n-m+1) / fact(n+m+1) * &
                        & fact(j+n-m+k+1) * fact(j+n+m-k+1) * vplm(indmk) * &
                        & vscales(indmk)
                    if (mod(abs(k+abs(m)-abs(k-m)), 4) .eq. 2) then
                        tmpk1 = -tmpk1
                    end if
                    tmpk2 = src_m(indn+abs(m)) * cphi
                    if ((m .ge. 0) .and. (m .le. k)) then
                        sphi = -sphi
                    end if
                    tmpk3 = -src_m(indn+abs(m)) * sphi
                    if (m .ne. k) then
                        tmpk1 = tmpk1 / sqrt_2
                    end if
                    if (m .ne. 0) then
                        tmpk1 = tmpk1 / sqrt_2
                        tmpk2 = tmpk2 + src_m(indn-abs(m))*sphi
                        tmpk3 = tmpk3 + src_m(indn-abs(m))*cphi
                    end if
                    if (m .lt. 0) then
                        tmpk3 = -tmpk3
                    end if
                    dst_l(indj+k) = dst_l(indj+k) + tmpk1*tmpk2
                    if (k .ne. 0) then
                        dst_l(indj-k) = dst_l(indj-k) + tmpk1*tmpk3
                    end if
                end do
            end do
        end do
    end do
end subroutine fmm_m2l_baseline

! Translate local expansion to another sphere
subroutine fmm_l2l_baseline(c, src_r, dst_r, p, vscales, src_l, dst_l)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of local spherical harmonics
!   vscales: normalization constants for Y_lm (of degree up to p)
!   src_l: expansion in old harmonics
!   dst_l: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1))
    real(kind=8), intent(in) :: src_l((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(inout) :: dst_l((p+1)*(p+1))
    real(kind=8) :: r, r1, r2, ctheta, stheta, cphi, sphi
    real(kind=8) :: vcos(p+1), vsin(p+1)
    real(kind=8) :: vplm((p+1)*(p+1)), fact(2*p+1), tmpk1, tmpk2, tmpk3
    real(kind=8) :: tmp1, tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, m, indj, indmk, indn, indjn
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
        r1 = r / src_r
        r2 = dst_r / r
        pow_r1(1) = r1
        pow_r2(1) = r2
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
                tmp1 = pow_r2(j+1) / fact(j-k+1) / fact(j+k+1)
                if (k .ne. 0) then
                    tmp1 = tmp1 * sqrt_2
                end if
                do n = j, p
                    indn = n*n + n + 1
                    indjn = (n-j)**2 + (n-j) + 1
                    tmp2 = tmp1 * pow_r1(n+1) * vscales(indj) / &
                        & vscales(indjn) / vscales(indn)
                    if (mod(n+j, 2) .eq. 1) then
                        tmp2 = -tmp2
                    end if
                    do m = k+j-n, k+n-j
                        indmk = indjn + abs(m-k)
                        cphi = vcos(1+abs(m-k))
                        sphi = vsin(1+abs(m-k))
                        tmpk1 = tmp2 * fact(n-m+1) * fact(n+m+1) / &
                            & fact(n-j-m+k+1) / fact(n-j+m-k+1) * &
                            & vplm(indmk) * vscales(indmk)
                        if (mod(abs(k+abs(m-k)-abs(m)), 4) .eq. 2) then
                            tmpk1 = -tmpk1
                        end if
                        tmpk2 = src_l(indn+abs(m)) * cphi
                        if ((m .ge. 0) .and. (m .le. k)) then
                            sphi = -sphi
                        end if
                        tmpk3 = -src_l(indn+abs(m)) * sphi
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt_2
                            tmpk2 = tmpk2 + src_l(indn-abs(m))*sphi
                            tmpk3 = tmpk3 + src_l(indn-abs(m))*cphi
                        end if
                        if (m .lt. 0) then
                            tmpk3 = -tmpk3
                        end if
                        dst_l(indj+k) = dst_l(indj+k) + tmpk1*tmpk2
                        if (k .ne. 0) then
                            dst_l(indj-k) = dst_l(indj-k) + tmpk1*tmpk3
                        end if
                    end do
                end do
            end do
        end do
    else
        r1 = dst_r / src_r
        tmpk1 = r1
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_l(k) = dst_l(k) + src_l(k)*tmpk1
            end do
            tmpk1 = tmpk1 * r1
        end do
    end if
end subroutine fmm_l2l_baseline

! Integrate spherical harmonics (grid -> coefficients)
subroutine int_grid(p, ngrid, w, vgrid, x, xlm)
! Parameters:
!   p: maximum degree of spherical harmonics
!   ngrid: number of Lebedev grid points
!   w: weights of Lebedev grid
!   vgrid: values of spherical harmonics at Lebedev grid points
!   x: values at grid points
!   xlm: resulting weights of spherical harmonics
    integer, intent(in) :: p, ngrid
    real(kind=8), intent(in) :: w(ngrid), vgrid((p+1)*(p+1), ngrid), x(ngrid)
    real(kind=8), intent(out) :: xlm((p+1)*(p+1))
    integer :: i
    xlm = 0
    do i = 1, ngrid
        xlm = xlm + x(i)*w(i)*vgrid(:,i)
    end do
end subroutine int_grid

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
        & rnode, snode)
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
!   snode: which node is leaf and contains only given sphere
    integer, intent(in) :: nsph
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph)
    integer, intent(inout) :: ind(nsph)
    integer, intent(out) :: cluster(2, 2*nsph-1), children(2, 2*nsph-1)
    integer, intent(out) :: parent(2*nsph-1)
    real(kind=8), intent(out) :: cnode(3, 2*nsph-1), rnode(2*nsph-1)
    integer, intent(out) :: snode(nsph)
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
            snode(ind(s)) = i
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

! Find near and far admissible pairs of tree nodes
subroutine tree_get_m2l(nsph, children, cnode, rnode, lwork, iwork, &
        & jwork, work, nnfar, nfar, nnnear, nnear)
! Parameters:
!   nsph: number of leaf nodes in a tree. total number of nodes is 2*nsph-1
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   lwork: size of work array in dimension 2
!   iwork: index of current pair of nodes that needs to be checked for
!       admissibility. must be 0 for first call of this subroutine. if on exit
!       iwork is less or equal to jwork, that means lwork was too small, please
!       reallocate work array and copy all the values into new array and then
!       run procedure again.
!   jwork: amount of stored possible admissible pairs of nodes.
!       Please read iwork comments.
!   work: all the far and near pairs will be stored here
!   nnfar: total amount of far admissible pairs. valid only if iwork is
!       greater than jwork on exit.
!   nfar: amount of far admissible pairs for each node. valid only if iwork is
!       greater than jwork on exit.
!   nnnear: total amount of near admissible pairs. valid only if iwork is
!       greater than jwork on exit
!   nnear: amount of near admissible pairs for each node. valid only if iwork
!       is greater than jwork on exit
    integer, intent(in) :: nsph, children(2, 2*nsph-1), lwork
    real(kind=8), intent(in) :: cnode(3, 2*nsph-1), rnode(2*nsph-1)
    integer, intent(inout) :: iwork, jwork, work(3, lwork)
    integer, intent(out) :: nnfar, nfar(2*nsph-1), nnnear, nnear(2*nsph-1)
    integer :: j(2)
    real(kind=8) :: c(3), r, d
    if (iwork .eq. 0) then
        work(1, 1) = 1
        work(2, 1) = 1
        iwork = 1
        jwork = 1
    end if
    do while (iwork .le. jwork)
        j = work(1:2, iwork)
        c = cnode(:, j(1)) - cnode(:, j(2))
        r = rnode(j(1)) + rnode(j(2))
        d = sqrt(c(1)*c(1) + c(2)*c(2) + c(3)*c(3))
        if (d .ge. 2*r) then
            ! Mark as far admissible pair
            !write(*,*) "FAR:", j
            work(3, iwork) = 1
        else if ((children(1, j(1)) .eq. 0) .or. &
            & (children(1, j(2)) .eq. 0)) then
            ! Mark as near admissible pair if one of nodes is a leaf node
            !write(*,*) "NEAR:", j
            work(3, iwork) = 2
        else if (jwork+4 .gt. lwork) then
            ! Exit procedure, since work array was too small
            !write(*,*) "SMALL LWORK"
            return
        else
            ! Mark as non-admissible pair and check all pairs of children nodes
            work(3, iwork) = 0
            work(1, jwork+1) = children(1, j(1))
            work(2, jwork+1) = children(1, j(2))
            work(1, jwork+2) = children(1, j(1))
            work(2, jwork+2) = children(2, j(2))
            work(1, jwork+3) = children(2, j(1))
            work(2, jwork+3) = children(1, j(2))
            work(1, jwork+4) = children(2, j(1))
            work(2, jwork+4) = children(2, j(2))
            jwork = jwork + 4
            !write(*,*) "NON:", j
        end if
        iwork = iwork + 1
    end do
    nfar = 0
    nnear = 0
    do iwork = 1, jwork
        if (work(3, iwork) .eq. 1) then
            nfar(work(1, iwork)) = nfar(work(1, iwork)) + 1
        else if (work(3, iwork) .eq. 2) then
            nnear(work(1, iwork)) = nnear(work(1, iwork)) + 1
        end if
    end do
    iwork = jwork + 1
    nnfar = sum(nfar)
    nnnear = sum(nnear)
end subroutine tree_get_m2l

! Get near and far admissible pairs from work array of tree_get_m2l
subroutine tree_get_far_near(jwork, lwork, work, nsph, nfar, sfar, far, &
        nnear, snear, near)
! Parameters:
    integer, intent(in) :: jwork, lwork, work(3, lwork), nsph, nfar, nnear
    integer, intent(out) :: sfar(nsph+1), far(nfar), snear(nsph+1), near(nnear)
    integer :: i
    do i = 1, jwork
        if (work(3, i) .eq. 1) then
            ! Far
        else if (work(3, i) .eq. 2) then
            ! Near
        end if
    end do
end subroutine tree_get_far_near

! Compute multipole coefficients for each node of tree
subroutine tree_m2m_baseline(nsph, p, vscales, coef_sph, ind, cluster, &
        & children, cnode, rnode, coef_node)
! Parameters:
!   nsph: Number of all spheres
!   p: maximum degree of multipole basis functions
!   vscales: normalization constants for Y_lm
!   coef_sph: multipole coefficients of input spheres
!   ind: permutation of all spheres (to localize sequential spheres)
!   cluster: first and last spheres (from ind array), belonging to each cluster
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   coef_node: multipole coefficients of bounding spheres of nodes
    integer, intent(in) :: nsph, p, ind(nsph), cluster(2, 2*nsph-1)
    real(kind=8), intent(in) :: vscales((p+1)*(p+1))
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
            call fmm_m2m_baseline(c1-c, r1, r, p, vscales, &
                & coef_node(:, j(1)), coef_node(:, i))
            call fmm_m2m_baseline(c2-c, r2, r, p, vscales, &
                & coef_node(:, j(2)), coef_node(:, i))
        end if
    end do
end subroutine tree_m2m_baseline

! Apply M2P from entire tree to a given point
subroutine tree_m2p_baseline(c, leaf, p, vscales, nclusters, children, cnode, &
        & rnode, coef_node, v)
! Parameters:
!   c: coordinate where to compute potential
!   leaf: which leaf node contains given sphere with grid point c
!   p: maximum degree of multipole basis functions
!   vscales: normalization constants for Y_lm
!   nclusters: number of nodes in a tree
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   coef_node: multipole coefficients of bounding spheres of nodes
!   v: value of potential
    real(kind=8), intent(in) :: c(3), vscales((p+1)*(p+1)), cnode(3, nclusters)
    integer, intent(in) :: p, leaf, nclusters, children(2, nclusters)
    real(kind=8), intent(in) :: rnode(nclusters)
    real(kind=8), intent(in) :: coef_node((p+1)*(p+1), nclusters)
    real(kind=8), intent(out) :: v
    integer :: i, far(nclusters), j(2)
    real(kind=8) :: d(3), r, tmp_v
    far = 0
    ! far(i): 0 if M2P must be ignored, 1 if M2P must be applied and 2 if need
    ! to decide if M2P is applicable
    far(1) = 2
    v = 0
    do i = 1, nclusters
        ! If M2P is not applicable, ignore node
        if (far(i) .eq. 0) then
            cycle
        end if
        j = children(:, i)
        d = c - cnode(:, i)
        ! If c belongs to the origin leaf, ignore M2P
        if (leaf .eq. i) then
            far(i) = 0
            cycle
        ! Apply M2P for other leaf nodes (c is always outside for them)
        else if (j(1) .eq. 0) then
            far(i) = 1
        ! Check if point is outside sphere then apply M2P otherwise check
        ! hierarchically
        else
            r = sqrt(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
            ! constant 2 guarantees good convergence and accuracy
            if (r .gt. 2*rnode(i)) then
                far(i) = 1
            else
                far(i) = 0
                far(j(1)) = 2
                far(j(2)) = 2
            end if
        end if
        ! If M2P is needed
        if (far(i) .eq. 1) then
            call fmm_m2p(d, rnode(i), p, vscales, coef_node(:, i), tmp_v)
            v = v + tmp_v
        end if
    end do
end subroutine tree_m2p_baseline

! Apply matvec for ddPCM spherical harmonics
subroutine pcm_matvec_grid(nsph, csph, rsph, ngrid, grid, w, vgrid, ui, p, &
        & vscales, ind, cluster, children, cnode, rnode, snode, coef_sph, &
        & coef_out)
! Parameters:
!   nsph: number of all spheres
!   csph: centers of all spheres
!   rsph: radiuses of all spheres
!   ngrid: number of Lebedev grid points on each sphere
!   grid: coordinates of Lebedev grid points on a unit sphere
!   ui: "outside" factor of grid points (0 is inside, 1 is "fully" outside)
!   p: maximum degree of multipole basis functions
!   vscales: normalization constants for Y_lm
!   ind: permutation of all spheres (to localize sequential spheres)
!   cluster: first and last spheres (from ind array), belonging to each cluster
!   children: children of each cluster. 0 means no children
!   cnode: center of bounding sphere of each cluster (node) of tree
!   rnode: radius of bounding sphere of each cluster (node) of tree
!   snode: which node is leaf and contains only given sphere
!   coef_sph: multipole coefficients of bounding spheres of nodes
!   coef_out: output multipole coefficients of spherical harmonics
    integer, intent(in) :: nsph, ngrid, p, ind(nsph), cluster(2, 2*nsph-1)
    integer, intent(in) :: children(2, 2*nsph-1), snode(nsph)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph), grid(3, ngrid)
    real(kind=8), intent(in) :: w(ngrid), vgrid((p+1)*(p+1), ngrid)
    real(kind=8), intent(in) :: ui(ngrid, nsph), vscales((p+1)*(p+1))
    real(kind=8), intent(in) :: cnode(3, 2*nsph-1), rnode(2*nsph-1)
    real(kind=8), intent(in) :: coef_sph((p+1)*(p+1), nsph)
    real(kind=8), intent(out) :: coef_out((p+1)*(p+1), nsph)
    real(kind=8) :: coef_sph_scaled((p+1)*(p+1), nsph)
    real(kind=8) :: coef_node((p+1)*(p+1), 2*nsph-1), c(3), x(ngrid)
    integer :: i, j, leaf, indi, indj
    do i = 0, p
        indi = i*i + i + 1
        do j = -i, i
            indj = indi + j
            coef_sph_scaled(indj, :) = i * coef_sph(indj, :) * rsph
        end do
    end do
    call tree_m2m_baseline(nsph, p, vscales, coef_sph_scaled, ind, cluster, &
        & children, cnode, rnode, coef_node)
    do i = 1, nsph
        leaf = snode(i)
        do j = 1, ngrid
            if (ui(j, i) .eq. 0) then
                x(j) = 0
            else
                c = csph(:, i) + rsph(i)*grid(:,j)
                call tree_m2p_baseline(c, leaf, p, vscales, 2*nsph-1, &
                    & children, cnode, rnode, coef_node, x(j))
                x(j) = ui(j, i) * x(j)
            end if
        end do
        call int_grid(p, ngrid, w, vgrid, x, coef_out(:, i))
    end do
end subroutine

end module

