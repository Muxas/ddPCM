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
            tmp = -tmp / sqrt(dble((l-m+1)*(l+m)))
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
            tmp = -tmp / sqrt(dble((l-m+1)*(l+m)))
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
    real(kind=8) :: plm((p+1)*(p+1)), t, vscales((p+1)*(p+1))
    real(kind=8) :: fact(2*p+1), pow_r(p+1)
    complex(kind=8) :: tmpk, tmpm
    integer :: j, k, n, m, indj, indk, indm, indn, indjn, indjna
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
                        if ((j .eq. 1) .and. (k .eq. -1)) then
                            !write(*,*) j, k, n, m, src_m(indjn+k-m), &
                            !    & fact(j-k+1), fact(j+k+1), 1/fact(n+m+1), &
                            !    & 1/fact(n-m+1), &
                            !    & 1/fact(j-n-k+m+1), 1/fact(j-n+k-m+1), &
                            !    & plm(indm), tmpm, pow_r(n+1), vscales(indm)
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
subroutine o2o_baseline_old(c, src_r, dst_r, p, src_m, dst_m)
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
    real(kind=8) :: plm((p+1)*(p+1)), t, vscales((p+1)*(p+1)), tmp, rcoef
    real(kind=8) :: fact(2*p+1), tmpk, tmpk1, tmpk2, tmpm1, tmpm2, sqrt_2
    real(kind=8) :: pow_r1(p+1), pow_r2(p+1), sqrt_four_pi
    integer :: j, k, n, m, indj, indk, indm, indn, indjn, indjna
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
                        if (k .ge. 0) then
                            tmpm1 = ncos(1+abs(m))
                            tmpm2 = nsin(1+abs(m))
                        else
                            tmpm1 = -nsin(1+abs(m))
                            tmpm2 = ncos(1+abs(m))
                        end if
                        if (mod(abs(abs(k)-abs(m)-abs(k-m)), 4) .eq. 2) then
                            tmpm1 = -tmpm1
                            tmpm2 = -tmpm2
                        end if
                        tmpk1 = fact(j-k+1) * fact(j+k+1) / fact(n-m+1) / &
                            & fact(n+m+1) / fact(j-n-k+m+1) / &
                            & fact(j-n+k-m+1) * plm(indm) * pow_r1(j-n+1) * &
                            & pow_r2(n+1) * vscales(indm) * vscales(indj) / &
                            & vscales(indjn) / vscales(indn)
                        if (m .ne. 0) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        if (k .ne. 0) then
                            tmpk1 = tmpk1 * sqrt_2
                        end if
                        if (m .ne. k) then
                            tmpk1 = tmpk1 / sqrt_2
                        end if
                        if ((m .ge. 0) .and. (k .ge. 0) .and. (k .lt. m)) then
                            tmpm2 = -tmpm2
                        end if
                        if (k .eq. m) then
                            tmpk2 = src_m(indjn) * tmpm1
                        else
                            tmpk2 = src_m(indjn+abs(k-m))*tmpm1 + &
                                & src_m(indjn-abs(k-m))*tmpm2
                        end if
                        tmpk = tmpk + tmpk1*tmpk2
                        if ((j .eq. 2) .and. ((k .eq. 2))) then
                            write(*,*) j, k, n, m, tmpk1, tmpk2, src_m(indjn+abs(k-m)), &
                                & tmpm1, src_m(indjn-abs(k-m)), tmpm2
                        end if
                    end do
                end do
                dst_m(indk) = tmpk
            end do
        end do
    else
        dst_m = src_m
    end if
end subroutine o2o_baseline_old

! O2O baseline translation
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
    real(kind=8) :: plm((p+1)*(p+1)), t, vscales((p+1)*(p+1)), tmp, rcoef
    real(kind=8) :: fact(2*p+1), tmpk, tmpk1, tmpk2, tmpk3, sqrt_2
    real(kind=8) :: pow_r1(p+1), pow_r2(p+1), sqrt_four_pi
    integer :: j, k, n, m, indj, indk, indm, indn, indjn, indjna
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
                            !sphi = -sphi
                            !write(*,*) j, k, n, m, tmpk1
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
                            if ((j .eq. 3) .and. (k .eq. 1)) then
                                write(*,*) j, k, n, m, 1, tmpk1*tmpk3
                            end if
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
                            !sphi = -sphi
                            !write(*,*) j, k, n, m, tmpk1
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
                            if ((j .eq. 3) .and. (k .eq. 1)) then
                                write(*,*) j, k, n, m, 2, tmpk1*tmpk3
                            end if
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
                            !sphi = -sphi
                            !write(*,*) j, k, n, m, tmpk1
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
                            if ((j .eq. 3) .and. (k .eq. 1)) then
                                write(*,*) j, k, n, m, 3, -tmpk1*tmpk3
                            end if
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
                            !sphi = -sphi
                            !write(*,*) j, k, n, m, tmpk1
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
                            if ((j .eq. 3) .and. (k .eq. 1)) then
                                write(*,*) j, k, n, m, 4, -tmpk1*tmpk3
                            end if
                        end if
                    end do
                end do
            end do
        end do
    else
        dst_m = src_m
    end if
end subroutine o2o_baseline

! Generate M2M matrix for spherical harmonics
subroutine o2o_matrix(c, src_r, dst_r, p, mat)
! Parameters:
!   c: radius-vector from new to old centers of harmonics
!   src_r: radius of old harmonics
!   dst_r: radius of new harmonics
!   p: maximum degree of spherical harmonics
!   mat: output matrix
    real(kind=8), intent(in) :: c(3), src_r, dst_r
    integer, intent(in) :: p
    real(kind=8), intent(out) :: mat((p+1)*(p+1), (p+1)*(p+1))
    real(kind=8) :: r, r1, r2, ctheta, stheta, cphi, sphi, ncos(p+1), nsin(p+1)
    real(kind=8) :: plm((p+1)*(p+1)), t, vscales((p+1)*(p+1)), tmp, rcoef
    real(kind=8) :: fact(2*p+1), tmpk, tmpm, sqrt_2, sqrt_four_pi
    real(kind=8) :: pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, m, indj, indk, indm, indn, indjn, indjna
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
                            tmpm = ncos(1-m)
                        else
                            tmpm = -nsin(m+1)
                        end if
                        if (mod(abs(abs(k)-abs(m)-abs(k-m)), 4) .eq. 2) then
                            tmpm = -tmpm
                        end if
                        mat(indk, indjn+k-m) = fact(j-k+1)*fact(j+k+1)/ &
                            & fact(n+m+1)**2/fact(j-n-k+m+1)/fact(j-n+k-m+1)* &
                            & plm(indm)*tmpm*pow_r2(n+1)*pow_r1(j-n+1)/ &
                            & vscales(indjn)*vscales(indj)
                    end do
                end do
            end do
        end do
    else
        do j = 1, (p+1)*(p+1)
            mat(j, j) = 1
        end do
    end if
end subroutine

end module
