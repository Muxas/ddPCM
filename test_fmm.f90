module fmm_pcm
implicit none
contains

! Compute cos(nx) and sin(nx) from a given cos(x) and sin(x)
subroutine trgev(c, s, p, ncos, nsin)
! Parameters:
!   c: cosinus of phi angle
!   s: sinus of phi angle
!   p: maximum degree of polynomials
!   ncos: values of cos(n*phi)
!   nsin: values of sin(n*phi)
    real(kind=8), intent(in) :: c, s
    integer, intent(in) :: p
    real(kind=8), intent(out) :: ncos(p+1), nsin(p+1)
    integer :: n

    ncos(1) = 1
    nsin(1) = 0
    do n = 2, p+1
        ncos(n) = ncos(n-1)*c - nsin(n-1)*s
        nsin(n) = ncos(n-1)*s + nsin(n-1)*c
    end do
    ! Later code corresponds to how it is done in ddCOSMO/ddPCM
    !ncos(2) = c
    !nsin(2) = s
    !do n = 3, p+1
    !    ncos(n) = 2*c*ncos(n-1) - ncos(n-2)
    !    nsin(n) = 2*c*nsin(n-1) - nsin(n-2)
    !end do
end subroutine trgev

! Compute associated Legendre polynomials P(l,m)(c)
! Uses recurrence formula
!   (l-m)p(l,m) = x(2l-1)p(l-1,m) - (l+m-1)p(l-2,m)
subroutine polleg(c, s, p, plm)
! Parameters:
!   c: cosinus of theta angle (where to compute polynomials)
!   s: sinus of theta angle
!   p: maximum degree of polynomials
!   plm: values of associated legendre polynomials P(l,m)(c)
    real(kind=8), intent(in) :: c, s
    integer, intent(in) :: p
    real(kind=8), intent(out) :: plm((p+1)*(p+1))
    integer :: m, ind, l, ind2
    real(kind=8) :: fact, pmm, pmm1, pmmo, pll, fm, fl

    fact  = 1
    pmm   = 1
    do m = 0, p 
        ind      = (m + 1)*(m + 1)
        plm(ind) = pmm
        if(m.eq.p) return
        fm = dble(m)
        pmm1 = c*(2*fm + 1)*pmm
        ind2 = ind + 2*m + 2
        plm(ind2) = pmm1
        pmmo = pmm
        do l = m+2, p
            fl = dble(l)
            pll   = (c*(2*fl - 1)*pmm1 - (fl + fm - 1)*pmm)/(fl - fm)
            ind = l*l + l + 1
            plm(ind+m) = pll
            pmm  = pmm1
            pmm1 = pll
        end do
        pmm  = -pmmo*fact*s
        fact = fact + 2
    end do
end subroutine polleg

! Compute scaling factors for non-normalized complex spherical harmonics
subroutine scales_complex(p, val)
! Parameters:
!   p: maximum degree of polynomial
!   val: corresponding scaling factors
    integer, intent(in) :: p
    real(kind=8), intent(out) :: val((p+1)*(p+1))
    real(kind=8) :: tmp
    integer :: l, ind, m
    do l = 0, p
        ! m = 0
        ind = l*l + l + 1
        val(ind) = 1
        tmp = 1
        ! m != 0
        do m = 1, l
            tmp = tmp / sqrt(dble((l-m+1)*(l+m)))
            val(ind+m) = tmp
            ! Negative , is ignored
            !val(ind-m) = tmp
        end do
    end do
end subroutine scales_complex

! Compute scaling factors for normalized complex spherical harmonics
subroutine scales_complex_normal(p, val)
! Parameters:
!   p: maximum degree of polynomial
!   val: corresponding scaling factors
    integer, intent(in) :: p
    real(kind=8), intent(out) :: val((p+1)*(p+1))
    real(kind=8) :: sqrt_four_pi, tmp
    integer :: l, ind, m
    sqrt_four_pi = 4*sqrt(atan(dble(1)))
    do l = 0, p
        ! m = 0
        ind = l*l + l + 1
        val(ind) = sqrt(dble(2*l+1)) / sqrt_four_pi
        tmp = val(ind)
        ! m != 0
        do m = 1, l
            tmp = -tmp / sqrt(dble((l-m+1)*(l+m)))
            val(ind+m) = tmp
            ! Negative , is ignored
            !val(ind-m) = tmp
        end do
    end do
end subroutine scales_complex_normal

! Compute scaling factors for non-normalized real spherical harmonics
subroutine scales_real(p, val)
! Parameters:
!   p: maximum degree of polynomial
!   val: corresponding scaling factors
    integer, intent(in) :: p
    real(kind=8), intent(out) :: val((p+1)*(p+1))
    real(kind=8) :: sqrt_2, tmp
    integer :: l, ind, m
    sqrt_2 = sqrt(dble(2))
    do l = 0, p
        ! m = 0
        ind = l*l + l + 1
        val(ind) = 1
        tmp = sqrt_2
        ! m != 0
        do m = 1, l
            tmp = tmp / sqrt(dble((l-m+1)*(l+m)))
            val(ind+m) = tmp
            ! Negative , is ignored
            !val(ind-m) = tmp
        end do
    end do
end subroutine scales_real

! Compute scaling factors for normalized real spherical harmonics
subroutine scales_real_normal(p, val)
! Parameters:
!   p: maximum degree of polynomial
!   val: corresponding scaling factors
    integer, intent(in) :: p
    real(kind=8), intent(out) :: val((p+1)*(p+1))
    real(kind=8) :: sqrt_2, sqrt_four_pi, tmp
    integer :: l, ind, m
    sqrt_2 = sqrt(dble(2))
    sqrt_four_pi = 4*sqrt(atan(dble(1)))
    do l = 0, p
        ! m = 0
        ind = l*l + l + 1
        val(ind) = sqrt(dble(2*l+1)) / sqrt_four_pi
        tmp = val(ind) * sqrt_2
        ! m != 0
        do m = 1, l
            tmp = -tmp / sqrt(dble((l-m+1)*(l+m)))
            val(ind+m) = tmp
            ! Negative , is ignored
            !val(ind-m) = tmp
        end do
    end do
end subroutine scales_real_normal

! Compute multipole coefficients for particle of unit charge
! Based on non-normalized complex spherical harmonics
subroutine p2m_complex(c, p, m)
! parameters:
!   c: coordinates of charged particle (relative to center of harmonics)
!   p: maximum power multipole basis function
!   m: multipole coefficients
    real(kind=8), intent(in) :: c(3)
    integer, intent(in) :: p
    complex(kind=8), intent(out) :: m((p+1)*(p+1))
    real(kind=8) :: r, ctheta, stheta, cphi, sphi, ncos(p+1), nsin(p+1)
    real(kind=8) :: plm((p+1)*(p+1)), t, vscales((p+1)*(p+1)), tmp
    integer :: n, k, ind
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    if (r .ne. 0) then
        ctheta = c(3) / r
        if (stheta .ne. 0) then
            stheta = sqrt(stheta)
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / r
            call trgev(cphi, sphi, p, ncos, nsin)
        else
            cphi = 1
            sphi = 0
            ncos = 1
            nsin = 0
        end if
        call polleg(ctheta, stheta, p, plm)
        call scales_complex(p, vscales)
        ! Now build harmonics to fill multipole coefficients
        t = 1
        do n = 0, p
            ind = n*n + n + 1
            m(ind) = t * vscales(ind) * plm(ind)
            do k = 1, n
                tmp = t * vscales(ind+k) * plm(ind+k)
                m(ind-k) = tmp * cmplx(ncos(k+1), nsin(k+1), 8)
                m(ind+k) = conjg(m(ind-k))
                !m(ind+k) = tmp * cmplx(ncos(k+1), -nsin(k+1), 8)
            end do
            t = t * r
        end do
    else
        m(1) = 1
        m(2:) = 0
    end if
end subroutine p2m_complex

! Compute potential, induced by non-normalized complex spherical harmonics
subroutine m2p_complex(c, p, m, v)
! Parameters:
!   c: relative distance from center of harmonics to point of potential
!   p: maximum degree of polynomials
!   m: multipole expansion at origin
!   v: value of induced potential
    real(kind=8), intent(in) :: c(3)
    complex(kind=8), intent(in) :: m((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: v
    complex(kind=8) :: cmplx_v, tmp, tmp1
    real(kind=8) :: r, t, ctheta, stheta, cphi, sphi, ncos(p+1), nsin(p+1)
    real(kind=8) :: vscales((p+1)*(p+1)), plm((p+1)*(p+1))
    integer :: n, k, ind
    t = 1
    cmplx_v = 0
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    ! r is always > 0
    ctheta = c(3) / r
    if (stheta .ne. 0) then
        stheta = sqrt(stheta)
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        stheta = stheta / r
        call trgev(cphi, sphi, p, ncos, nsin)
    else
        cphi = 1
        sphi = 0
        ncos = 1
        nsin = 0
    end if
    call polleg(ctheta, stheta, p, plm)
    call scales_complex(p, vscales)
    do n = 0, p
        t = t / r
        ind = n*n + n + 1
        ! k = 0
        tmp = m(ind) * plm(ind) * vscales(ind)
        ! k != 0
        do k = 1, n
            tmp1 = plm(ind+k) * vscales(ind+k) * &
                & cmplx(ncos(k+1), nsin(k+1), 8)
            tmp = tmp + m(ind+k)*tmp1 + m(ind-k)*conjg(tmp1)
        end do
        cmplx_v = cmplx_v + t*tmp
    end do
    v = real(cmplx_v)
end subroutine m2p_complex

! Compute multipole coefficients for particle of unit charge
! Based on normalized complex spherical harmonics
subroutine p2m_complex_normal(c, p, m)
! parameters:
!   c: coordinates of charged particle (relative to center of harmonics)
!   p: maximum power multipole basis function
!   m: multipole coefficients
    real(kind=8), intent(in) :: c(3)
    integer, intent(in) :: p
    complex(kind=8), intent(out) :: m((p+1)*(p+1))
    real(kind=8) :: r, ctheta, stheta, cphi, sphi, ncos(p+1), nsin(p+1)
    real(kind=8) :: plm((p+1)*(p+1)), t, vscales((p+1)*(p+1)), tmp
    integer :: n, k, ind
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    if (r .ne. 0) then
        ctheta = c(3) / r
        if (stheta .ne. 0) then
            stheta = sqrt(stheta)
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / r
            call trgev(cphi, sphi, p, ncos, nsin)
        else
            cphi = 1
            sphi = 0
            ncos = 1
            nsin = 0
        end if
        call polleg(ctheta, stheta, p, plm)
        call scales_complex_normal(p, vscales)
        ! Now build harmonics to fill multipole coefficients
        t = 1
        do n = 0, p
            ind = n*n + n + 1
            m(ind) = t * vscales(ind) * plm(ind)
            do k = 1, n
                tmp = t * vscales(ind+k) * plm(ind+k)
                m(ind-k) = tmp * cmplx(ncos(k+1), nsin(k+1), 8)
                m(ind+k) = conjg(m(ind-k))
                !m(ind+k) = tmp * cmplx(ncos(k+1), -nsin(k+1), 8)
            end do
            t = t * r
        end do
    else
        m(1) = vscales(1)
        m(2:) = 0
    end if
end subroutine p2m_complex_normal

! Compute potential, induced by normalized complex spherical harmonics
subroutine m2p_complex_normal(c, p, m, v)
! Parameters:
!   c: relative distance from center of harmonics to point of potential
!   p: maximum degree of polynomials
!   m: multipole expansion at origin
!   v: value of induced potential
    real(kind=8), intent(in) :: c(3)
    complex(kind=8), intent(in) :: m((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: v
    complex(kind=8) :: cmplx_v, tmp, tmp1
    real(kind=8) :: r, t, ctheta, stheta, cphi, sphi, ncos(p+1), nsin(p+1)
    real(kind=8) :: vscales((p+1)*(p+1)), plm((p+1)*(p+1))
    integer :: n, k, ind
    t = 1
    cmplx_v = 0
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    ! r is always > 0
    ctheta = c(3) / r
    if (stheta .ne. 0) then
        stheta = sqrt(stheta)
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        stheta = stheta / r
        call trgev(cphi, sphi, p, ncos, nsin)
    else
        cphi = 1
        sphi = 0
        ncos = 1
        nsin = 0
    end if
    call polleg(ctheta, stheta, p, plm)
    call scales_complex_normal(p, vscales)
    do n = 0, p
        t = t / r
        ind = n*n + n + 1
        ! k = 0
        tmp = m(ind) * plm(ind) / vscales(ind)
        ! k != 0
        do k = 1, n
            tmp1 = plm(ind+k) * vscales(ind+k) / vscales(ind)**2 * &
                & cmplx(ncos(k+1), nsin(k+1), 8)
            tmp = tmp + m(ind+k)*tmp1 + m(ind-k)*conjg(tmp1)
        end do
        cmplx_v = cmplx_v + t*tmp
    end do
    v = real(cmplx_v)
end subroutine m2p_complex_normal

! Compute multipole coefficients for particle of unit charge
! Based on non-normalized real spherical harmonics
subroutine p2m_real(c, p, m)
! parameters:
!   c: coordinates of charged particle (relative to center of harmonics)
!   p: maximum power multipole basis function
!   m: multipole coefficients
    real(kind=8), intent(in) :: c(3)
    integer, intent(in) :: p
    real(kind=8), intent(out) :: m((p+1)*(p+1))
    real(kind=8) :: r, ctheta, stheta, cphi, sphi, ncos(p+1), nsin(p+1)
    real(kind=8) :: plm((p+1)*(p+1)), t, vscales((p+1)*(p+1)), tmp
    integer :: n, k, ind
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    if (r .ne. 0) then
        ctheta = c(3) / r
        if (stheta .ne. 0) then
            stheta = sqrt(stheta)
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / r
            call trgev(cphi, sphi, p, ncos, nsin)
        else
            cphi = 1
            sphi = 0
            ncos = 1
            nsin = 0
        end if
        call polleg(ctheta, stheta, p, plm)
        call scales_real(p, vscales)
        ! Now build harmonics to fill multipole coefficients
        t = 1
        do n = 0, p
            ind = n*n + n + 1
            m(ind) = t * vscales(ind) * plm(ind)
            do k = 1, n
                tmp = t * vscales(ind+k) * plm(ind+k)
                m(ind+k) = tmp * ncos(k+1)
                m(ind-k) = -(tmp * nsin(k+1))
            end do
            t = t * r
        end do
    else
        m(1) = 1
        m(2:) = 0
    end if
end subroutine p2m_real

! Compute potential, induced by non-normalized real spherical harmonics
subroutine m2p_real(c, p, m, v)
! Parameters:
!   c: relative distance from center of harmonics to point of potential
!   p: maximum degree of polynomials
!   m: multipole expansion at origin
!   v: value of induced potential
    real(kind=8), intent(in) :: c(3)
    real(kind=8), intent(in) :: m((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: v
    real(kind=8) :: cmplx_v, tmp, tmp1
    real(kind=8) :: r, t, ctheta, stheta, cphi, sphi, ncos(p+1), nsin(p+1)
    real(kind=8) :: vscales((p+1)*(p+1)), plm((p+1)*(p+1))
    integer :: n, k, ind
    t = 1
    cmplx_v = 0
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    ! r is always > 0
    ctheta = c(3) / r
    if (stheta .ne. 0) then
        stheta = sqrt(stheta)
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        stheta = stheta / r
        call trgev(cphi, sphi, p, ncos, nsin)
    else
        cphi = 1
        sphi = 0
        ncos = 1
        nsin = 0
    end if
    call polleg(ctheta, stheta, p, plm)
    call scales_real(p, vscales)
    do n = 0, p
        t = t / r
        ind = n*n + n + 1
        ! k = 0
        tmp = m(ind) * plm(ind) * vscales(ind)
        ! k != 0
        do k = 1, n
            tmp1 = plm(ind+k) * vscales(ind+k)
            tmp = tmp + m(ind+k)*tmp1*ncos(k+1) - m(ind-k)*tmp1*nsin(k+1)
        end do
        cmplx_v = cmplx_v + t*tmp
    end do
    v = cmplx_v
end subroutine m2p_real

! Compute multipole coefficients for particle of unit charge
! Based on normalized real spherical harmonics
subroutine p2m_real_normal(c, p, m)
! parameters:
!   c: coordinates of charged particle (relative to center of harmonics)
!   p: maximum power multipole basis function
!   m: multipole coefficients
    real(kind=8), intent(in) :: c(3)
    integer, intent(in) :: p
    real(kind=8), intent(out) :: m((p+1)*(p+1))
    real(kind=8) :: r, ctheta, stheta, cphi, sphi, ncos(p+1), nsin(p+1)
    real(kind=8) :: plm((p+1)*(p+1)), t, vscales((p+1)*(p+1)), tmp
    integer :: n, k, ind
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    if (r .ne. 0) then
        ctheta = c(3) / r
        if (stheta .ne. 0) then
            stheta = sqrt(stheta)
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / r
            call trgev(cphi, sphi, p, ncos, nsin)
        else
            cphi = 1
            sphi = 0
            ncos = 1
            nsin = 0
        end if
        call polleg(ctheta, stheta, p, plm)
        call scales_real_normal(p, vscales)
        ! Now build harmonics to fill multipole coefficients
        t = 1
        do n = 0, p
            ind = n*n + n + 1
            m(ind) = t * vscales(ind) * plm(ind)
            do k = 1, n
                tmp = t * vscales(ind+k) * plm(ind+k)
                m(ind+k) = tmp * ncos(k+1)
                m(ind-k) = -(tmp * nsin(k+1))
            end do
            t = t * r
        end do
    else
        m(1) = 1
        m(2:) = 0
    end if
end subroutine p2m_real_normal

! Compute potential, induced by normalized real spherical harmonics
subroutine m2p_real_normal(c, p, m, v)
! Parameters:
!   c: relative distance from center of harmonics to point of potential
!   p: maximum degree of polynomials
!   m: multipole expansion at origin
!   v: value of induced potential
    real(kind=8), intent(in) :: c(3)
    real(kind=8), intent(in) :: m((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: v
    real(kind=8) :: cmplx_v, tmp, tmp1
    real(kind=8) :: r, t, ctheta, stheta, cphi, sphi, ncos(p+1), nsin(p+1)
    real(kind=8) :: vscales((p+1)*(p+1)), plm((p+1)*(p+1))
    integer :: n, k, ind
    t = 1
    cmplx_v = 0
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    ! r is always > 0
    ctheta = c(3) / r
    if (stheta .ne. 0) then
        stheta = sqrt(stheta)
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        stheta = stheta / r
        call trgev(cphi, sphi, p, ncos, nsin)
    else
        cphi = 1
        sphi = 0
        ncos = 1
        nsin = 0
    end if
    call polleg(ctheta, stheta, p, plm)
    call scales_real_normal(p, vscales)
    do n = 0, p
        t = t / r
        ind = n*n + n + 1
        ! k = 0
        tmp = m(ind) * plm(ind) / vscales(ind)
        ! k != 0
        do k = 1, n
            tmp1 = plm(ind+k) * vscales(ind+k) / vscales(ind)**2
            tmp = tmp + m(ind+k)*tmp1*ncos(k+1) - m(ind-k)*tmp1*nsin(k+1)
        end do
        cmplx_v = cmplx_v + t*tmp
    end do
    v = cmplx_v
end subroutine m2p_real_normal

! Compute multipole coefficients for particle of unit charge
! Based on normalized real spherical harmonics or given radius
subroutine p2o(c, rho, p, m)
! parameters:
!   c: coordinates of charged particle (relative to center of harmonics)
!   rho: radius of spherical harmonics
!   p: maximum power multipole basis function
!   m: multipole coefficients
    real(kind=8), intent(in) :: c(3), rho
    integer, intent(in) :: p
    real(kind=8), intent(out) :: m((p+1)*(p+1))
    real(kind=8) :: r, ctheta, stheta, cphi, sphi, ncos(p+1), nsin(p+1)
    real(kind=8) :: plm((p+1)*(p+1)), t, vscales((p+1)*(p+1)), tmp, rcoef
    integer :: n, k, ind
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    if (r .ne. 0) then
        ctheta = c(3) / r
        if (stheta .ne. 0) then
            stheta = sqrt(stheta)
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / r
            call trgev(cphi, sphi, p, ncos, nsin)
        else
            cphi = 1
            sphi = 0
            ncos = 1
            nsin = 0
        end if
        call polleg(ctheta, stheta, p, plm)
        call scales_real_normal(p, vscales)
        ! Now build harmonics to fill multipole coefficients
        rcoef = r / rho
        t = 1
        do n = 0, p
            ind = n*n + n + 1
            m(ind) = t * vscales(ind) * plm(ind)
            do k = 1, n
                tmp = t * vscales(ind+k) * plm(ind+k)
                m(ind+k) = tmp * ncos(k+1)
                m(ind-k) = -(tmp * nsin(k+1))
            end do
            t = t * rcoef
        end do
    else
        m(1) = 1
        m(2:) = 0
    end if
end subroutine p2o

! Compute potential, induced by normalized real spherical harmonics with radius
subroutine o2p(c, rho, p, m, v)
! Parameters:
!   c: relative distance from center of harmonics to point of potential
!   rho: radius of spherical harmonics
!   p: maximum degree of polynomials
!   m: multipole expansion at origin
!   v: value of induced potential
    real(kind=8), intent(in) :: c(3), rho
    real(kind=8), intent(in) :: m((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: v
    real(kind=8) :: cmplx_v, tmp, tmp1, tmp2
    real(kind=8) :: r, t, ctheta, stheta, cphi, sphi, ncos(p+1), nsin(p+1)
    real(kind=8) :: vscales((p+1)*(p+1)), plm((p+1)*(p+1)), rcoef
    integer :: n, k, ind
    t = 1 / rho
    cmplx_v = 0
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    ! r is always > 0
    ctheta = c(3) / r
    if (stheta .ne. 0) then
        stheta = sqrt(stheta)
        cphi = c(1) / stheta
        sphi = c(2) / stheta
        stheta = stheta / r
        call trgev(cphi, sphi, p, ncos, nsin)
    else
        cphi = 1
        sphi = 0
        ncos = 1
        nsin = 0
    end if
    call polleg(ctheta, stheta, p, plm)
    call scales_real_normal(p, vscales)
    rcoef = rho / r
    do n = 0, p
        t = t * rcoef
        ind = n*n + n + 1
        ! k = 0
        tmp = m(ind) * plm(ind) / vscales(ind)
        ! k != 0
        do k = 1, n
            tmp1 = plm(ind+k) * vscales(ind+k) / vscales(ind)**2
            tmp2 = m(ind+k)*ncos(k+1) - m(ind-k)*nsin(k+1)
            tmp = tmp + tmp1*tmp2
        end do
        cmplx_v = cmplx_v + t*tmp
    end do
    v = cmplx_v
end subroutine o2p


! M2M baseline translation
subroutine m2m_baseline(c, p, src_m, dst_m)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   p: maximum degree of spherical harmonics
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    real(kind=8), intent(in) :: c(3)
    integer, intent(in) :: p
    complex(kind=8), intent(in) :: src_m((p+1)*(p+1))
    complex(kind=8), intent(out) :: dst_m((p+1)*(p+1))
    real(kind=8) :: r, ctheta, stheta, cphi, sphi, ncos(p+1), nsin(p+1)
    real(kind=8) :: plm((p+1)*(p+1)), vscales((p+1)*(p+1))
    real(kind=8) :: fact(2*p+1), pow_r(p+1)
    complex(kind=8) :: tmpk, tmpm
    integer :: j, k, n, m, indj, indk, indm, indn, indjn
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    if (r .ne. 0) then
        ctheta = c(3) / r
        if (stheta .ne. 0) then
            stheta = sqrt(stheta)
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / r
            call trgev(cphi, sphi, p, ncos, nsin)
        else
            cphi = 1
            sphi = 0
            ncos = 1
            nsin = 0
        end if
        call polleg(ctheta, stheta, p, plm)
        call scales_complex(p, vscales)
        pow_r(1) = 1
        do j = 2, p+1
            pow_r(j) = pow_r(j-1) * r
        end do
        ! Fill square roots of factorials
        fact(1) = 1
        do j = 2, 2*p+1
            fact(j) = sqrt(dble(j-1)) * fact(j-1)
        end do
        do j = 0, p
            indj = j*j + j + 1
            do k = -j, j
                indk = indj + k
                tmpk = 0
                do n = 0, j
                    indn = n*n + n + 1
                    indjn = (j-n)**2 + (j-n) + 1
                    do m = -n, n
                        indm = indn + abs(m)
                        if (j-n .lt. abs(k-m)) then
                            cycle
                        end if
                        if (m .le. 0) then
                            tmpm = cmplx(ncos(1-m), nsin(1-m), 8)
                        else
                            tmpm = cmplx(ncos(m+1), -nsin(m+1), 8)
                        end if
                        if (mod(abs(abs(k)-abs(m)-abs(k-m)), 4) .eq. 2) then
                            tmpm = -tmpm
                        end if
                        tmpk = tmpk + src_m(indjn+k-m)*fact(j-k+1)* &
                            & fact(j+k+1)/fact(n+m+1)/fact(n-m+1)/ &
                            & fact(j-n-k+m+1)/fact(j-n+k-m+1)*plm(indm)*tmpm* &
                            & pow_r(n+1)*vscales(indm)
                    end do
                end do
                dst_m(indk) = tmpk
            end do
        end do
    else
        dst_m = src_m
    end if
end subroutine m2m_baseline

! O2O baseline translation
! Basline in terms of operation count: p^4, while it can be done in p^3
subroutine o2o_baseline(c, src_r, dst_r, p, src_m, dst_m)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   src_m: expansion in old harmonics
!   dst_m: expansion in new harmonics
    real(kind=8), intent(in) :: c(3), src_r, dst_r, src_m((p+1)*(p+1))
    integer, intent(in) :: p
    real(kind=8), intent(out) :: dst_m((p+1)*(p+1))
    real(kind=8) :: r, r1, r2, ctheta, stheta, cphi, sphi, ncos(p+1), nsin(p+1)
    real(kind=8) :: plm((p+1)*(p+1)), vscales((p+1)*(p+1))
    real(kind=8) :: fact(2*p+1), tmpk1, tmpk2, tmpk3, sqrt_2
    real(kind=8) :: pow_r1(p+1), pow_r2(p+1), sqrt_four_pi
    integer :: j, k, n, m, indj, indm, indn, indjn
    sqrt_2 = sqrt(dble(2))
    sqrt_four_pi = 4*sqrt(atan(dble(1)))
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    if (r .ne. 0) then
        ctheta = c(3) / r
        if (stheta .ne. 0) then
            stheta = sqrt(stheta)
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / r
            call trgev(cphi, sphi, p, ncos, nsin)
        else
            cphi = 1
            sphi = 0
            ncos = 1
            nsin = 0
        end if
        call polleg(ctheta, stheta, p, plm)
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
                dst_m(indj+k) = 0
                dst_m(indj-k) = 0
                do n = 0, j
                    indn = n*n + n + 1
                    indjn = (j-n)**2 + (j-n) + 1
                    do m = max(k+n-j,-n), -1
                        indm = indn + abs(m)
                        cphi = ncos(1+abs(m))
                        sphi = nsin(1+abs(m))
                        tmpk1 = fact(j-k+1) * fact(j+k+1) / fact(n-m+1) / &
                            & fact(n+m+1) / fact(j-n-k+m+1) / &
                            & fact(j-n+k-m+1) * plm(indm) * pow_r1(j-n+1) * &
                            & pow_r2(n+1) * vscales(indm) * vscales(indj) / &
                            & vscales(indjn) / vscales(indn)
                        if (mod(abs(abs(k)-abs(m)-abs(k-m)), 4) .eq. 2) then
                            tmpk1 = -tmpk1
                        end if
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        if (k .ne. 0) then
                            tmpk1 = tmpk1 * sqrt_2
                        end if
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        tmpk2 = src_m(indjn+abs(k-m))*cphi - &
                            & src_m(indjn-abs(k-m))*sphi
                        dst_m(indj+k) = dst_m(indj+k) + tmpk1*tmpk2
                        if (k .ne. 0) then
                            tmpk3 = src_m(indjn-abs(k-m))*cphi + &
                                & src_m(indjn+abs(k-m))*sphi
                            dst_m(indj-k) = dst_m(indj-k) + tmpk1*tmpk3
                        end if
                    end do
                    do m = max(0,k+n-j), min(k-1,n)
                        indm = indn + abs(m)
                        cphi = ncos(1+abs(m))
                        sphi = nsin(1+abs(m))
                        tmpk1 = fact(j-k+1) * fact(j+k+1) / fact(n-m+1) / &
                            & fact(n+m+1) / fact(j-n-k+m+1) / &
                            & fact(j-n+k-m+1) * plm(indm) * pow_r1(j-n+1) * &
                            & pow_r2(n+1) * vscales(indm) * vscales(indj) / &
                            & vscales(indjn) / vscales(indn)
                        if (mod(abs(abs(k)-abs(m)-abs(k-m)), 4) .eq. 2) then
                            tmpk1 = -tmpk1
                        end if
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        if (k .ne. 0) then
                            tmpk1 = tmpk1 * sqrt_2
                        end if
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        tmpk2 = src_m(indjn+abs(k-m))*cphi + &
                            & src_m(indjn-abs(k-m))*sphi
                        dst_m(indj+k) = dst_m(indj+k) + tmpk1*tmpk2
                        if (k .ne. 0) then
                            tmpk3 = src_m(indjn-abs(k-m))*cphi - &
                                & src_m(indjn+abs(k-m))*sphi
                            dst_m(indj-k) = dst_m(indj-k) + tmpk1*tmpk3
                        end if
                    end do
                    if (k .le. n) then
                        m = k
                        indm = indn + abs(m)
                        cphi = ncos(1+abs(m))
                        sphi = nsin(1+abs(m))
                        tmpk1 = fact(j-k+1) * fact(j+k+1) / fact(n-m+1) / &
                            & fact(n+m+1) / fact(j-n-k+m+1) / &
                            & fact(j-n+k-m+1) * plm(indm) * pow_r1(j-n+1) * &
                            & pow_r2(n+1) * vscales(indm) * vscales(indj) / &
                            & vscales(indjn) / vscales(indn)
                        if (mod(abs(abs(k)-abs(m)-abs(k-m)), 4) .eq. 2) then
                            tmpk1 = -tmpk1
                        end if
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        if (k .ne. 0) then
                            tmpk1 = tmpk1 * sqrt_2
                        end if
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        tmpk2 = src_m(indjn)*cphi
                        dst_m(indj+k) = dst_m(indj+k) + tmpk1*tmpk2
                        if (k .ne. 0) then
                            tmpk3 = src_m(indjn)*sphi
                            dst_m(indj-k) = dst_m(indj-k) - tmpk1*tmpk3
                        end if
                    end if
                    do m = k+1, min(j-n+k, n)
                        indm = indn + abs(m)
                        cphi = ncos(1+abs(m))
                        sphi = nsin(1+abs(m))
                        tmpk1 = fact(j-k+1) * fact(j+k+1) / fact(n-m+1) / &
                            & fact(n+m+1) / fact(j-n-k+m+1) / &
                            & fact(j-n+k-m+1) * plm(indm) * pow_r1(j-n+1) * &
                            & pow_r2(n+1) * vscales(indm) * vscales(indj) / &
                            & vscales(indjn) / vscales(indn)
                        if (mod(abs(abs(k)-abs(m)-abs(k-m)), 4) .eq. 2) then
                            tmpk1 = -tmpk1
                        end if
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        if (k .ne. 0) then
                            tmpk1 = tmpk1 * sqrt_2
                        end if
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        tmpk2 = src_m(indjn+abs(k-m))*cphi - &
                            & src_m(indjn-abs(k-m))*sphi
                        dst_m(indj+k) = dst_m(indj+k) + tmpk1*tmpk2
                        if (k .ne. 0) then
                            tmpk3 = src_m(indjn-abs(k-m))*cphi + &
                                & src_m(indjn+abs(k-m))*sphi
                            dst_m(indj-k) = dst_m(indj-k) - tmpk1*tmpk3
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
                dst_m(k) = src_m(k) * tmpk1
            end do
            tmpk1 = tmpk1 * r1
        end do
    end if
end subroutine o2o_baseline

! O2O baseline translation matrix
! Basline in terms of operation count: p^4, while it can be done in p^3
subroutine o2o_matrix(c, src_r, dst_r, p, mat)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   mat: m2m matrix
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    integer, intent(in) :: p
    real(kind=8), intent(out) :: mat((p+1)*(p+1), (p+1)*(p+1))
    real(kind=8) :: r, r1, r2, ctheta, stheta, cphi, sphi, ncos(p+1), nsin(p+1)
    real(kind=8) :: plm((p+1)*(p+1)), vscales((p+1)*(p+1))
    real(kind=8) :: fact(2*p+1), tmpk1, sqrt_2
    real(kind=8) :: pow_r1(p+1), pow_r2(p+1), sqrt_four_pi
    integer :: j, k, n, m, indj, indn, indm, indjn
    sqrt_2 = sqrt(dble(2))
    sqrt_four_pi = 4*sqrt(atan(dble(1)))
    stheta = c(1)*c(1) + c(2)*c(2)
    r = sqrt(c(3)*c(3) + stheta)
    mat = 0
    if (r .ne. 0) then
        ctheta = c(3) / r
        if (stheta .ne. 0) then
            stheta = sqrt(stheta)
            cphi = c(1) / stheta
            sphi = c(2) / stheta
            stheta = stheta / r
            call trgev(cphi, sphi, p, ncos, nsin)
        else
            cphi = 1
            sphi = 0
            ncos = 1
            nsin = 0
        end if
        call polleg(ctheta, stheta, p, plm)
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
                do n = 0, j
                    indn = n*n + n + 1
                    indjn = (j-n)**2 + (j-n) + 1
                    do m = max(k+n-j,-n), -1
                        indm = indn + abs(m)
                        cphi = ncos(1+abs(m))
                        sphi = nsin(1+abs(m))
                        tmpk1 = fact(j-k+1) * fact(j+k+1) / fact(n-m+1) / &
                            & fact(n+m+1) / fact(j-n-k+m+1) / &
                            & fact(j-n+k-m+1) * plm(indm) * pow_r1(j-n+1) * &
                            & pow_r2(n+1) * vscales(indm) * vscales(indj) / &
                            & vscales(indjn) / vscales(indn)
                        if (mod(abs(abs(k)-abs(m)-abs(k-m)), 4) .eq. 2) then
                            tmpk1 = -tmpk1
                        end if
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        if (k .ne. 0) then
                            tmpk1 = tmpk1 * sqrt_2
                        end if
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        !tmpk2 = src_m(indjn+abs(k-m))*cphi - &
                        !    & src_m(indjn-abs(k-m))*sphi
                        !dst_m(indj+k) = dst_m(indj+k) + tmpk1*tmpk2
                        mat(indj+k, indjn+abs(k-m)) = cphi*tmpk1
                        mat(indj+k, indjn-abs(k-m)) = -sphi*tmpk1
                        if (k .ne. 0) then
                            !tmpk3 = src_m(indjn-abs(k-m))*cphi + &
                            !    & src_m(indjn+abs(k-m))*sphi
                            !dst_m(indj-k) = dst_m(indj-k) + tmpk1*tmpk3
                            mat(indj-k, indjn+abs(k-m)) = sphi*tmpk1
                            mat(indj-k, indjn-abs(k-m)) = cphi*tmpk1
                        end if
                    end do
                    do m = max(0,k+n-j), min(k-1,n)
                        indm = indn + abs(m)
                        cphi = ncos(1+abs(m))
                        sphi = nsin(1+abs(m))
                        tmpk1 = fact(j-k+1) * fact(j+k+1) / fact(n-m+1) / &
                            & fact(n+m+1) / fact(j-n-k+m+1) / &
                            & fact(j-n+k-m+1) * plm(indm) * pow_r1(j-n+1) * &
                            & pow_r2(n+1) * vscales(indm) * vscales(indj) / &
                            & vscales(indjn) / vscales(indn)
                        if (mod(abs(abs(k)-abs(m)-abs(k-m)), 4) .eq. 2) then
                            tmpk1 = -tmpk1
                        end if
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        if (k .ne. 0) then
                            tmpk1 = tmpk1 * sqrt_2
                        end if
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        !tmpk2 = src_m(indjn+abs(k-m))*cphi + &
                        !    & src_m(indjn-abs(k-m))*sphi
                        !dst_m(indj+k) = dst_m(indj+k) + tmpk1*tmpk2
                        mat(indj+k, indjn+abs(k-m)) = cphi*tmpk1
                        mat(indj+k, indjn-abs(k-m)) = sphi*tmpk1
                        if (k .ne. 0) then
                            !tmpk3 = src_m(indjn-abs(k-m))*cphi - &
                            !    & src_m(indjn+abs(k-m))*sphi
                            !dst_m(indj-k) = dst_m(indj-k) + tmpk1*tmpk3
                            !mat(indj-k, indjn+abs(k-m)) = &
                            !    & mat(indj-k, indjn+abs(k-m)) - sphi*tmpk1
                            !mat(indj-k, indjn-abs(k-m)) = &
                            !    & mat(indj-k, indjn-abs(k-m)) + cphi*tmpk1
                            mat(indj-k, indjn+abs(k-m)) = -sphi*tmpk1
                            mat(indj-k, indjn-abs(k-m)) = cphi*tmpk1
                        end if
                    end do
                    if (k .le. n) then
                        m = k
                        indm = indn + abs(m)
                        cphi = ncos(1+abs(m))
                        sphi = nsin(1+abs(m))
                        tmpk1 = fact(j-k+1) * fact(j+k+1) / fact(n-m+1) / &
                            & fact(n+m+1) / fact(j-n-k+m+1) / &
                            & fact(j-n+k-m+1) * plm(indm) * pow_r1(j-n+1) * &
                            & pow_r2(n+1) * vscales(indm) * vscales(indj) / &
                            & vscales(indjn) / vscales(indn)
                        if (mod(abs(abs(k)-abs(m)-abs(k-m)), 4) .eq. 2) then
                            tmpk1 = -tmpk1
                        end if
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        if (k .ne. 0) then
                            tmpk1 = tmpk1 * sqrt_2
                        end if
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        !tmpk2 = src_m(indjn)*cphi
                        !dst_m(indj+k) = dst_m(indj+k) + tmpk1*tmpk2
                        mat(indj+k, indjn) = cphi*tmpk1
                        if (k .ne. 0) then
                            !tmpk3 = src_m(indjn)*sphi
                            !dst_m(indj-k) = dst_m(indj-k) - tmpk1*tmpk3
                            mat(indj-k, indjn) = -sphi*tmpk1
                        end if
                    end if
                    do m = k+1, min(j-n+k, n)
                        indm = indn + abs(m)
                        cphi = ncos(1+abs(m))
                        sphi = nsin(1+abs(m))
                        tmpk1 = fact(j-k+1) * fact(j+k+1) / fact(n-m+1) / &
                            & fact(n+m+1) / fact(j-n-k+m+1) / &
                            & fact(j-n+k-m+1) * plm(indm) * pow_r1(j-n+1) * &
                            & pow_r2(n+1) * vscales(indm) * vscales(indj) / &
                            & vscales(indjn) / vscales(indn)
                        if (mod(abs(abs(k)-abs(m)-abs(k-m)), 4) .eq. 2) then
                            tmpk1 = -tmpk1
                        end if
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        if (k .ne. 0) then
                            tmpk1 = tmpk1 * sqrt_2
                        end if
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        !tmpk2 = src_m(indjn+abs(k-m))*cphi - &
                        !    & src_m(indjn-abs(k-m))*sphi
                        !dst_m(indj+k) = dst_m(indj+k) + tmpk1*tmpk2
                        mat(indj+k, indjn+abs(k-m)) = &
                            & mat(indj+k, indjn+abs(k-m)) + cphi*tmpk1
                        mat(indj+k, indjn-abs(k-m)) = &
                            & mat(indj+k, indjn-abs(k-m)) - sphi*tmpk1
                        if (k .ne. 0) then
                            !tmpk3 = src_m(indjn-abs(k-m))*cphi + &
                            !    & src_m(indjn+abs(k-m))*sphi
                            !dst_m(indj-k) = dst_m(indj-k) - tmpk1*tmpk3
                            mat(indj-k, indjn+abs(k-m)) = &
                                & mat(indj-k, indjn+abs(k-m)) - sphi*tmpk1
                            mat(indj-k, indjn-abs(k-m)) = &
                                & mat(indj-k, indjn-abs(k-m)) - cphi*tmpk1
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
                mat(k, k) = tmpk1
            end do
            tmpk1 = tmpk1 * r1
        end do
    end if
end subroutine

! Divide given cluster of spheres into two clusters by inertial bisection
subroutine geometry_divide_cluster(nsph, csph, n, ind, div)
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
    real(kind=8) :: c(3), tmpcsph(3, n), a(3, 3), w(3), work(9), scal(n)
    real(kind=8) :: alpha=1, beta=0
    integer :: i, l, r, lwork=9, info, tmp_ind(n)
    c = 0
    do i = 1, n
        c = c + csph(:, ind(i))
    end do
    c = c / n
    do i = 1, n
        tmpcsph(:, i) = csph(:, ind(i)) - c
    end do
    call dgemm('N', 'T', 3, 3, n, alpha, tmpcsph, 3, tmpcsph, 3, beta, a, 3)
    call dsyev('V', 'L', 3, a, 3, w, work, lwork, info)
    call dgemv('T', 3, n, alpha, tmpcsph, 3, a(:, 3), 1, beta, scal, 1)
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
end subroutine geometry_divide_cluster

! Divide hierarchically until cluster consists of a single sphere
! Number of clusters is always 2*nsph-1
subroutine geometry_divide_hierarchically(nsph, csph, ind, cluster, children, &
        & parent)
    integer, intent(in) :: nsph
    real(kind=8), intent(in) :: csph(3, nsph)
    integer, intent(inout) :: ind(nsph)
    integer, intent(out) :: cluster(2, 2*nsph-1), children(2, 2*nsph-1)
    integer, intent(out) :: parent(2*nsph-1)
    integer :: i, j, n, s, e, div
    cluster(1, 1) = 1
    cluster(2, 1) = nsph
    parent(1) = 0
    j = 2
    do i = 1, 2*nsph-1
        s = cluster(1, i)
        e = cluster(2, i)
        n = e - s + 1
        if (n .gt. 1) then
            call geometry_divide_cluster(nsph, csph, n, ind(s:e), div)
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
end subroutine geometry_divide_hierarchically

! Compute 
subroutine tree_o2o(nsph, p, csph, rsph, coef_sph, ind, cluster, children, &
        & ccluster, rcluster, coef_cluster)
    integer, intent(in) :: nsph, p, ind(nsph), cluster(2, 2*nsph-1)
    integer, intent(in) :: children(2, 2*nsph-1)
    real(kind=8), intent(in) :: csph(3, nsph), rsph(nsph)
    real(kind=8), intent(in) :: coef_sph((p+1)*(p+1), nsph)
    real(kind=8), intent(out) :: ccluster(3, 2*nsph-1), rcluster(2*nsph-1)
    real(kind=8), intent(out) :: coef_cluster((p+1)*(p+1), 2*nsph-1)
    integer :: i, j
    real(kind=8) :: c1(3), c2(3), c(3), r1, r2, r, d, tmp_coef((p+1)*(p+1))
    do i = 2*nsph-1, 1, -1
        if (children(1, i) .eq. 0) then
            j = cluster(1, i)
            ccluster(:, i) = csph(:, ind(j))
            rcluster(i) = rsph(ind(j))
            coef_cluster(:, i) = coef_sph(:, ind(j))
        else
            c1 = ccluster(:, children(1, i))
            r1 = rcluster(children(1, i))
            c2 = ccluster(:, children(2, i))
            r2 = rcluster(children(2, i))
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
            ccluster(:, i) = c
            rcluster(i) = r
            call o2o_baseline(c1-c, r1, r, p, &
                & coef_cluster(:, children(1, i)), tmp_coef)
            coef_cluster(:, i) = tmp_coef
            call o2o_baseline(c2-c, r2, r, p, &
                & coef_cluster(:, children(2, i)), tmp_coef)
            coef_cluster(:, i) = coef_cluster(:, i) + tmp_coef
        end if
    end do
end subroutine tree_o2o

end module
