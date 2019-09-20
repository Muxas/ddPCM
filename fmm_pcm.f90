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

! M2M translation of non-normalized complex spherical harmonics in p^4 flops
subroutine m2m_complex_baseline(c, p, src_m, dst_m)
! Parameters:
!   c: coordinate of old center relative to new center of harmonics
!   p: maximum degree of spherical harmonics
!   src_m: multipole coefficients at source
!   dst_m: multipole coefficients at destination
    real(kind=8), intent(in) :: c(3)
    integer, intent(in) :: p
    complex(kind=8), intent(in) :: src_m((p+1)*(p+1))
    complex(kind=8), intent(out) :: dst_m((p+1)*(p+1))
    real(kind=8) :: r, ctheta, stheta, cphi, sphi, ncos(p+1), nsin(p+1)
    real(kind=8) :: vscales(p+1), plm((p+1)*(p+1))
    integer :: j, k, n, m, ind, ind_neg, ind_pos
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
        do j = 0, p
            ! k = 0
            ! k > 0
            do k = 1, j
                do n = 0, j
                    ! m = 0
                    ! m > 0
                    do m = 1, n
                    end do
                end do
            end do
        end do
    else
        dst_m = src_m
    end if
end subroutine m2m_complex_baseline

! M2M translation of non-normalized complex spherical harmonics along OZ axis
! NOT YET WORKING WELL
subroutine m2m_complex_translate_z(z, p, src_m, dst_m)
! Parameters:
!   z: change in Z coordinate
!   p: maximum degree of spherical harmonics
!   src_m: multipole coefficients at source
!   dst_m: multipole coefficients at destination
    real(kind=8), intent(in) :: z
    integer, intent(in) :: p
    complex(kind=8), intent(in) :: src_m((p+1)*(p+1))
    complex(kind=8), intent(out) :: dst_m((p+1)*(p+1))
    integer :: l, j, i, m, ind, ind_neg, ind_pos
    real(kind=8) :: coef, sq2=sqrt(dble(2)), plm((p+1)*(p+1))
    do l = 0, p
        ! m = 0
        ind = l*l + l + 1
        dst_m(ind) = 0
        do j = 0, l
            coef = 1
            do i = 1, l-j
                coef = coef * z / i
            end do
            dst_m(ind) = dst_m(ind) + coef*src_m(j*j+j+1)
        end do
        ! m > 0
        do m = 1, l
            ind_neg = ind - m
            ind_pos = ind + m
            dst_m(ind_neg) = 0
            dst_m(ind_pos) = 0
            do j = m, l
                coef = 1
                do i = 1, l-j
                    coef = coef * z / i
                end do
                dst_m(ind_neg) = dst_m(ind_neg) + coef*src_m(j*j+j+1-m)
                dst_m(ind_pos) = dst_m(ind_pos) + coef*src_m(j*j+j+1+m)
            end do
        end do
    end do
end subroutine

end module
