module fmm_pcm
implicit none
contains

!-----------------------------------------------------------------------------
! Purpose : service routine for computation of spherical harmonics
!-----------------------------------------------------------------------------
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
    !ncos(2) = c
    !nsin(2) = s
    !do n = 3, p+1
    !    ncos(n) = 2*c*ncos(n-1) - ncos(n-2)
    !    nsin(n) = 2*c*nsin(n-1) - nsin(n-2)
    !end do
end subroutine trgev

!------------------------------------------------------------------------------
! Purpose : compute the l,m associated legendre polynomial for -1 <= x <= 1
!           using the recurrence formula
!               (l-m)p(l,m) = x(2l-1)p(l-1,m) - (l+m-1)p(l-2,m)
! COPIED FROM DDCOSMO
!------------------------------------------------------------------------------
subroutine polleg(c, s, p, plm)
! Parameters:
!   c: cosinus of theta angle
!   s: sinus of theta angle
!   p: maximum degree of polynomials
!   plm: values of associated legendre polynomials p(l,m)
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

! Precompute normalization constants for harmonics
subroutine scales_complex(p, val)
! Parameters:
!   p: maximum degree of polynomial
!   val: values of normalization constants
    integer, intent(in) :: p
    real(kind=8), intent(out) :: val((p+1)*(p+1))
    real(kind=8) :: sqrt_2, sqrt_four_pi, tmp
    integer :: l, ind, m
    sqrt_2 = sqrt(dble(2))
    sqrt_four_pi = 4*sqrt(atan(dble(1)))
    do l = 0, p
        ! m = 0
        ind = l*l + l + 1
        val(ind) = 1 !sqrt(dble(2*l+1)) / sqrt_four_pi
        !val(ind) = 1
        tmp = val(ind)! * sqrt_2
        ! m != 0
        do m = 1, l
            tmp = tmp / sqrt(dble((l-m+1)*(l+m)))
            val(ind+m) = tmp
            ! Negative , is ignored
            !val(ind-m) = tmp
        end do
    end do
end subroutine scales_complex

! Precompute normalization constants for harmonics
subroutine scales_real(p, val)
! Parameters:
!   p: maximum degree of polynomial
!   val: values of normalization constants
    integer, intent(in) :: p
    real(kind=8), intent(out) :: val((p+1)*(p+1))
    real(kind=8) :: sqrt_2, sqrt_four_pi, tmp
    integer :: l, ind, m
    sqrt_2 = sqrt(dble(2))
    sqrt_four_pi = 4*sqrt(atan(dble(1)))
    do l = 0, p
        ! m = 0
        ind = l*l + l + 1
        val(ind) = 1 !sqrt(dble(2*l+1)) / sqrt_four_pi
        !val(ind) = 1
        tmp = val(ind) * sqrt_2
        ! m != 0
        do m = 1, l
            tmp = tmp / sqrt(dble((l-m+1)*(l+m)))
            val(ind+m) = tmp
            ! Negative , is ignored
            !val(ind-m) = tmp
        end do
    end do
end subroutine scales_real

! Compute multipole coefficients for given charged particle (unit charge)
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

! Compute potential, induced by multipole harmonics
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

! Compute multipole coefficients for given charged particle (unit charge)
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
                m(ind-k) = tmp * nsin(k+1)
            end do
            t = t * r
        end do
    else
        m(1) = 1
        m(2:) = 0
    end if
end subroutine p2m_real

! Compute potential, induced by multipole harmonics
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
            tmp = tmp + m(ind+k)*tmp1*ncos(k+1) + m(ind-k)*tmp1*nsin(k+1)
        end do
        cmplx_v = cmplx_v + t*tmp
    end do
    v = cmplx_v
end subroutine m2p_real

end module
