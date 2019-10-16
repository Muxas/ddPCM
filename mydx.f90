module mydx
    use ddcosmo
    use ddpcm_lib, only : dx, dodiag
    implicit none
    real*8  :: tobohr
    real*8, parameter :: toang = 0.52917721092d0, tokcal = 627.509469d0
    real*8, allocatable :: x(:), y(:), z(:), rvdw(:), charge(:)

contains
    subroutine init
        implicit none
        integer :: i, n
        dodiag = .false.
        open (unit=100,file='Input.txt',form='formatted',access='sequential')
        read(100,*) iprint      ! printing flag
        read(100,*) nproc       ! number of openmp threads
        read(100,*) lmax        ! max angular momentum of spherical harmonics basis
        read(100,*) ngrid       ! number of lebedev points
        read(100,*) iconv       ! 10^(-iconv) is the convergence threshold for the iterative solver
        read(100,*) igrad       ! whether to compute (1) or not (0) forces
        read(100,*) eps         ! dielectric constant of the solvent
        read(100,*) eta         ! regularization parameter
        read(100,*) n           ! number of atoms
        allocate (x(n), y(n), z(n), rvdw(n), charge(n))
        do i = 1, n
            read(100,*) charge(i), x(i), y(i), z(i), rvdw(i)
        end do
        tobohr = 1.0d0 / toang
        x    = x * tobohr
        y    = y * tobohr
        z    = z * tobohr
        rvdw = rvdw *tobohr
        close (100)
        call ddinit(n, x, y, z, rvdw)
        deallocate (x, y, z, rvdw, charge)
    end subroutine init

    subroutine fini
        implicit none
        call memfree
    end subroutine fini

    subroutine get_sizes(nbasis_, ngrid_, nsph_)
        implicit none
        integer, intent(out) :: nbasis_, ngrid_, nsph_
        nbasis_ = nbasis
        ngrid_ = ngrid
        nsph_ = nsph
    end subroutine get_sizes

    subroutine get_spheres(nsph, csph_, rsph_)
        implicit none
        integer, intent(in) :: nsph
        real*8, intent(out) :: csph_(3, nsph), rsph_(nsph)
        csph_ = csph
        rsph_ = rsph
    end subroutine get_spheres

    subroutine get_ngrid_ext(ngrid_ext)
        implicit none
        integer, intent(out) :: ngrid_ext
        integer :: isph, igrid
        ngrid_ext = 0
        do isph = 1, nsph
            do igrid = 1, ngrid
                if(ui(igrid, isph) .ne. zero) then
                    ngrid_ext = ngrid_ext + 1
                end if
            end do
        end do
    end subroutine get_ngrid_ext

    subroutine get_grid(ngrid_ext, cgrid, grid_sph)
        implicit none
        integer, intent(in) :: ngrid_ext
        real*8, intent(out) :: cgrid(3, ngrid_ext)
        integer, intent(out) :: grid_sph(ngrid_ext)
        integer :: isph, igrid, igrid_ext
        igrid_ext = 0
        do isph = 1, nsph
            do igrid = 1, ngrid
                if(ui(igrid, isph) .ne. zero) then
                    igrid_ext = igrid_ext + 1
                    cgrid(:, igrid_ext) = csph(:, isph) + &
                        & rsph(isph)*grid(:, igrid)
                    grid_sph(igrid_ext) = isph
                end if
            end do
        end do
    end subroutine get_grid

    subroutine ddpcm_dx(nbasis, nsph, x, y)
        implicit none
        integer, intent(in) :: nbasis, nsph
        real*8, intent(in) :: x(nbasis, nsph)
        real*8, intent(out) :: y(nbasis, nsph)
        call dx(nbasis*nsph, x, y)
    end subroutine ddpcm_dx

    subroutine mypolleg(x, y, l, m, plm)
        implicit none
        real*8, intent(in) :: x, y
        integer, intent(in) :: l, m
        real*8, intent(out) :: plm
        integer :: ll, mm
        real*8 :: pmm, pll1m, pll2m, fact
        pmm = one
        fact = -one
        do mm = 1, m
            pmm = fact * y * pmm
            fact = fact - two
        end do
        if (m .eq. l) then
            plm = pmm
            return
        end if
        fact = -fact
        pll1m = x * fact * pmm
        if (m+1 .eq. l) then
            plm = pll1m
            return
        end if
        pll2m = pmm
        do ll = m+2, l
            fact = fact + two
            plm = x*fact*pll1m - dble(ll+m-1)*pll2m
            plm = plm / dble(ll-m)
            pll2m = pll1m
            pll1m = plm
        end do
    end subroutine mypolleg

    subroutine check_polleg()
        implicit none
        real*8 :: s(3), tmp, cthe, sthe
        real*8 :: basloc(nbasis), vplm(nbasis), vcos(lmax+1), vsin(lmax+1)
        integer :: l, m
        s = one / sqrt(3.0)
        cthe = s(3)
        sthe = sqrt(one - cthe*cthe)
        call ylmbas(s, basloc, vplm, vcos, vsin)
        do l = 0, lmax
            do m = 0, l
                call mypolleg(cthe, sthe, l, m, tmp)
                write (*, *) l, m, vplm(l*l+l+1+m), tmp
            end do
        end do
    end subroutine check_polleg

    subroutine kernel(src, src_r, dst, ibasis, mat)
        implicit none
        integer, intent(in) :: ibasis
        real*8, intent(in) :: src(3), src_r, dst(3)
        real*8, intent(out) :: mat
        integer :: i, l, ind, m
        real*8 :: v(3), vv, t, s(3), tt, f, cthe, sthe
        real*8 :: basloc(nbasis), vplm(nbasis), vcos(lmax+1), vsin(lmax+1)
        real*8 :: fourpi
        real*8 :: tmp
        complex*16 :: phi, mphi
        fourpi = four * pi
        v = dst - src
        vv = sqrt(dot_product(v, v))
        t = vv / src_r
        s = v / vv
        cthe = s(3)
        sthe = sqrt(s(1)*s(1) + s(2)*s(2))
        if (sthe .ne. zero) then
            phi = cmplx(s(1)/sthe, s(2)/sthe, 8)
        else
            phi = one
        end if
        ! build the local basis
        !call ylmbas(s, basloc, vplm, vcos, vsin)
        ! with all the required stuff, finally compute
        ! the "potential" at the point
        tt = one / t
        do l = 0, lmax
            ind = l*l + l + 1
            m = ibasis - ind
            if((m .ge. -l) .and. (m .le. l)) then
                f = fourpi * dble(l) / (two*dble(l) + one) * tt
                !mat = f * basloc(ibasis)
                !exit
                call mypolleg(cthe, sthe, l, abs(m), tmp)
                mphi = phi**abs(m)
                if(m .ge. zero) then
                    mat = f * tmp * facs(ibasis) * real(mphi)
                else
                    mat = f * tmp * facs(ibasis) * imag(mphi)
                end if
                exit
            end if
            tt = tt / t
        end do
    end subroutine kernel

    subroutine kernel_block(nharms, charms, rharms, iharms, harm_sph, nsrc, &
            & src, ngrid_ext, grid_ext, grid_sph, ndst, dst, mat)
        implicit none
        integer, intent(in) :: nharms, iharms(nharms), harm_sph(nharms), nsrc
        integer, intent(in) :: src(nsrc), ngrid_ext, grid_sph(ngrid_ext), ndst
        integer, intent(in) :: dst(ndst)
        real*8, intent(in) :: charms(3, nharms), rharms(nharms)
        real*8, intent(in) :: grid_ext(3, ngrid_ext)
        real*8, intent(out) :: mat(ndst, nsrc)
        integer :: idst, i, isrc, j
        do isrc = 1, nsrc
            i = src(isrc) + 1
            do idst = 1, ndst
                j = dst(idst) + 1
                !write(*, *) "OUT: ", i-1, j-1, grid_sph(j)
                if(harm_sph(i) .eq. grid_sph(j)) then
                    mat(idst, isrc) = zero
                else
                    call kernel(charms(:, i), rharms(i), grid_ext(:, j), &
                        & iharms(i), mat(idst, isrc))
                end if
            end do
        end do
    end subroutine kernel_block

    subroutine kernel_point(src, src_r, dst, lmax, nbasis, mat)
        implicit none
        integer, intent(in) :: lmax, nbasis
        real*8, intent(in) :: src(3), dst(3)
        real*8, intent(in) :: src_r
        real*8, intent(out) :: mat(nbasis)
        integer :: i, l, ind, m
        real*8 :: v(3), vv, t, s(3), tt, f
        real*8 :: basloc(nbasis), vplm(nbasis), vcos(lmax+1), vsin(lmax+1)
        real*8 :: fourpi
        fourpi = four * pi
        v = dst - src
        vv = sqrt(dot_product(v, v))
        t = vv / src_r
        s = v / vv
        ! build the local basis
        call ylmbas(s, basloc, vplm, vcos, vsin)
        ! with all the required stuff, finally compute
        ! the "potential" at the point
        tt = one / t
        do l = 0, lmax
            ind = l*l + l + 1
            f = fourpi * dble(l) / (two*dble(l) + one) * tt
            do m = -l, l
                mat(ind+m) = f * basloc(ind+m)
            end do
            tt = tt / t
        end do
    end subroutine kernel_point

! My ylmbas
subroutine ylmbas2( x, basloc, lmax, nbasis )
!        
      implicit none
      integer, intent(in) :: lmax, nbasis
      real*8, dimension(3), intent(in) :: x
      real*8, dimension(nbasis), intent(out) :: basloc
      real*8 :: vplm(nbasis)
      real*8, dimension(lmax+1) :: vcos, vsin
!
      integer :: l, m, ind
      real*8  :: cthe, sthe, cphi, sphi, plm
!      
!------------------------------------------------------------------------------------------------
!
!     get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
!     coordinates of x.
!
!     evaluate cos( theta ) ; sin( theta )  
      cthe = x(3)
      sthe = sqrt(x(1)*x(1) + x(2)*x(2))
!
!     evalutate cos( phi ) ; sin( phi )
      if (sthe.ne.zero) then
        cphi = x(1)/sthe
        sphi = x(2)/sthe
      else
        cphi = zero
        sphi = zero
      endif
!
!     evaluate cos(m*phi) and sin(m*phi) arrays. notice that this is 
!     pointless if z = 1, as the only non vanishing terms will be the 
!     ones with m=0.
      if(sthe.ne.zero) then
        call trgev(cphi,sphi,vcos,vsin)
      else
        vcos = one
        vsin = zero
      endif
!
!     evaluate the generalized legendre polynomials
      call polleg(cthe,sthe,vplm)
!
!     now build the spherical harmonics. we will distinguish m=0,
!     m>0 and m<0:
      do l = 0, lmax
        ind = l**2 + l + 1
!
!       m = 0
        basloc(ind) = facs(ind)*vplm(ind)
        !basloc(ind) = vplm(ind)
!
        do m = 1, l
!        
          plm = vplm(ind+m)
!
!         m > 0
          basloc(ind+m) = facs(ind+m)*plm*vcos(m+1)
          !basloc(ind+m) = plm*vcos(m+1)
!
!         m < 0
          basloc(ind-m) = facs(ind-m)*plm*vsin(m+1)
          !basloc(ind-m) = plm*vsin(m+1)
!
        enddo
      enddo
!      
!      
    end subroutine ylmbas2

    subroutine kernel2_point(src, src_r, dst, lmax, nbasis, mat)
        implicit none
        integer, intent(in) :: lmax, nbasis
        real*8, intent(in) :: src(3), dst(3)
        real*8, intent(in) :: src_r
        real*8, intent(out) :: mat(nbasis)
        integer :: i, l, ind, m
        real*8 :: v(3), vv, t, s(3), tt, f
        real*8 :: basloc(nbasis), vplm(nbasis), vcos(lmax+1), vsin(lmax+1)
        real*8 :: fourpi
        fourpi = four * pi
        v = dst - src
        vv = sqrt(dot_product(v, v))
        t = vv / src_r
        s = v / vv
        ! build the local basis
        call ylmbas(s, basloc, vplm, vcos, vsin)
        ! with all the required stuff, finally compute
        ! the "potential" at the point
        tt = one / t
        do l = 0, lmax
            ind = l*l + l + 1
            f = tt
            do m = -l, l
                mat(ind+m) = f * basloc(ind+m)
            end do
            tt = tt / t
        end do
    end subroutine kernel2_point

    subroutine kernel_point_block(nsph, csph, rsph, nsrc, src, ngrid_ext, &
            & grid_ext, grid_sph, ndst, dst, lmax, nbasis, mat)
        implicit none
        integer, intent(in) :: nsph, nsrc, src(nsrc), ngrid_ext, lmax
        integer, intent(in) :: grid_sph(ngrid_ext), ndst, dst(ndst), nbasis
        real*8, intent(in) :: csph(3, nsph), rsph(nsph), grid_ext(3, ngrid_ext)
        real*8, intent(out) :: mat(ndst, nbasis, nsrc)
        integer :: idst, i, isrc, j
        do isrc = 1, nsrc
            i = src(isrc) + 1
            do idst = 1, ndst
                j = dst(idst) + 1
                !write(*, *) "OUT: ", i-1, j-1, grid_sph(j)
                if(i-1 .eq. grid_sph(j)) then
                    mat(idst, :, isrc) = zero
                else
                    call kernel_point(csph(:, i), rsph(i), grid_ext(:, j), &
                        & lmax, nbasis, mat(idst, :, isrc))
                end if
            end do
        end do
    end subroutine kernel_point_block

    subroutine kernel_point_ngrid(src_c, src_r, dst_c, dst_r, lmax, nbasis, &
            & ngrid, mat)
        implicit none
        integer, intent(in) :: lmax, nbasis, ngrid
        real*8, intent(in) :: src_c(3), dst_c(3)
        real*8, intent(in) :: src_r, dst_r
        real*8, intent(out) :: mat(ngrid, nbasis)
        integer :: i, l, ind, m
        real*8 :: v(3), vi(3), vvi, ti, si(3), tt, f
        real*8 :: basloc(nbasis), vplm(nbasis), vcos(lmax+1), vsin(lmax+1)
        real*8 :: fourpi
        fourpi = four * pi
        do i = 1, ngrid
            vi = dst_c + dst_r*grid(:, i)
            call kernel_point(src_c, src_r, vi, lmax, nbasis, mat(i, :))
        end do
    end subroutine kernel_point_ngrid

    subroutine kernel_ngrid(src_c, src_r, dst_c, dst_r, lmax, nbasis, ngrid, &
            & mat)
        implicit none
        integer, intent(in) :: lmax, nbasis, ngrid
        real*8, intent(in) :: src_c(3), dst_c(3)
        real*8, intent(in) :: src_r, dst_r
        real*8, intent(out) :: mat(ngrid, nbasis)
        integer :: i, l, ind, m
        real*8 :: v(3), vi(3), vvi, ti, si(3), tt, f
        real*8 :: basloc(nbasis), vplm(nbasis), vcos(lmax+1), vsin(lmax+1)
        real*8 :: fourpi
        fourpi = four * pi
        v = dst_c - src_c
        do i = 1, ngrid
            vi = v + dst_r*grid(:, i)
            vvi = sqrt(dot_product(vi, vi))
            ti = vvi / src_r
            si = vi / vvi
            ! build the local basis
            call ylmbas(si, basloc, vplm, vcos, vsin)
            ! with all the required stuff, finally compute
            ! the "potential" at the point
            tt = one / ti
            do l = 0, lmax
                ind = l*l + l + 1
                f = fourpi * dble(l) / (two*dble(l) + one) * tt
                do m = -l, l
                    mat(i, ind+m) = f * basloc(ind+m)
                end do
                tt = tt / ti
            end do
        end do
    end subroutine kernel_ngrid

    subroutine gen_mat_kernel(nbasis, ngrid, nsph, mat)
        implicit none
        integer, intent(in) :: nbasis, ngrid, nsph
        real*8, intent(out) :: mat(ngrid, nsph, nbasis, nsph)
        real*8, allocatable :: vts(:), vplm(:), basloc(:), vcos(:), vsin(:)
        real*8 :: c(3), vij(3), sij(3)
        real*8 :: vvij, tij, fourpi, tt, f, f1
        integer :: its, isph, jsph, l, m, ind, lm, istatus

        allocate(vts(ngrid), vplm(nbasis), basloc(nbasis), vcos(lmax+1), &
            & vsin(lmax+1), stat=istatus)
        if(istatus .ne. 0) then
            write(6,*) 'gen_mat: allocation failed !'
            stop
        end if
        fourpi = four * pi

        do isph = 1, nsph
            do jsph = 1, nsph
                if(isph .eq. jsph) then
                    mat(:, isph, :, jsph) = zero
                    cycle
                end if
                call kernel_ngrid(csph(:, jsph), rsph(jsph), csph(:, isph), &
                    & rsph(isph), lmax, nbasis, ngrid, mat(:, isph, :, jsph))
                !do its = 1, ngrid
                !    mat(its, isph, :, jsph) = ui(its, isph) * &
                !        & mat(its, isph, :, jsph)
                !end do
            end do
        end do
        deallocate(vts, vplm, basloc, vcos, vsin, stat=istatus)
        if(istatus .ne. 0) then
            write(6,*) 'gen_mat: deallocation failed !'
            stop
        end if
    end subroutine gen_mat_kernel

    subroutine gen_mat_ngrid(nbasis, ngrid, nsph, mat)
        implicit none
        integer, intent(in) :: nbasis, ngrid, nsph
        real*8, intent(out) :: mat(ngrid, nsph, nbasis, nsph)
        real*8, allocatable :: vts(:), vplm(:), basloc(:), vcos(:), vsin(:)
        real*8 :: c(3), vij(3), sij(3)
        real*8 :: vvij, tij, fourpi, tt, f, f1
        integer :: its, isph, jsph, l, m, ind, lm, istatus

        allocate(vts(ngrid), vplm(nbasis), basloc(nbasis), vcos(lmax+1), &
            & vsin(lmax+1), stat=istatus)
        if(istatus .ne. 0) then
            write(6,*) 'gen_mat: allocation failed !'
            stop
        end if
        fourpi = four * pi

        do isph = 1, nsph
            do jsph = 1, nsph
                if(isph .eq. jsph) then
                    mat(:, isph, :, jsph) = zero
                    cycle
                end if
                do its = 1, ngrid
                    !if(ui(its, isph) .eq. zero) then
                    !    mat(its, isph, :, jsph) = zero
                    !    cycle
                    !end if
                    ! build the geometrical variables
                    c = csph(:, isph) + rsph(isph)*grid(:, its)
                    vij = c - csph(:, jsph)
                    vvij = sqrt(dot_product(vij, vij))
                    tij = vvij / rsph(jsph)
                    sij = vij / vvij 
                    ! build the local basis
                    call ylmbas(sij, basloc, vplm, vcos, vsin)
                    ! with all the required stuff, finally compute
                    ! the "potential" at the point
                    tt = one / tij
                    do l = 0, lmax
                        ind = l*l + l + 1
                        f = fourpi * dble(l) / (two*dble(l) + one) * tt! * &
                            !& ui(its, isph)
                        do m = -l, l
                            mat(its, isph, ind+m, jsph) = f * basloc(ind+m)
                            ! write(*, *) its, isph, ind+m, jsph, mat(its, isph, ind+m, jsph)
                        end do
                        tt = tt / tij
                    end do
                end do
            end do
        end do
        deallocate(vts, vplm, basloc, vcos, vsin, stat=istatus)
        if(istatus .ne. 0) then
            write(6,*) 'gen_mat: deallocation failed !'
            stop
        end if
    end subroutine gen_mat_ngrid

    subroutine gen_mat(nbasis, nsph, mat)
        implicit none
        integer, intent(in) :: nbasis, nsph
        real*8, intent(out) :: mat(nbasis, nsph, nbasis, nsph)
        real*8, allocatable :: vts(:,:), vplm(:), basloc(:), vcos(:), vsin(:)
        real*8 :: c(3), vij(3), sij(3)
        real*8 :: vvij, tij, fourpi, tt, f, f1
        integer :: its, isph, jsph, l, m, ind, lm, istatus

        allocate(vts(ngrid, nbasis), vplm(nbasis), basloc(nbasis), &
            & vcos(lmax+1), vsin(lmax+1), stat=istatus)
        if(istatus .ne. 0) then
            write(6,*) 'gen_mat: allocation failed !'
            stop
        end if
        fourpi = four * pi

        do isph = 1, nsph
            do jsph = 1, nsph
                if(isph .eq. jsph) then
                    mat(:, isph, :, jsph) = zero
                    cycle
                end if
                do its = 1, ngrid
                    if(ui(its, isph) .eq. zero) then
                        vts(its, :) = zero
                        !mat(its, isph, :, jsph) = zero
                        cycle
                    end if
                    ! build the geometrical variables
                    c = csph(:, isph) + rsph(isph)*grid(:, its)
                    vij = c - csph(:, jsph)
                    vvij = sqrt(dot_product(vij, vij))
                    tij = vvij / rsph(jsph)
                    sij = vij / vvij 
                    ! build the local basis
                    call ylmbas(sij, basloc, vplm, vcos, vsin)
                    ! with all the required stuff, finally compute
                    ! the "potential" at the point
                    tt = one / tij
                    do l = 0, lmax
                        ind = l*l + l + 1
                        f = fourpi * dble(l) / (two*dble(l) + one) * tt * &
                            & ui(its, isph)
                        do m = -l, l
                            vts(its, ind+m) = f * basloc(ind+m)
                            !mat(its, isph, ind+m, jsph) = f * basloc(ind+m)
                            ! write(*, *) its, isph, ind+m, jsph, mat(its, isph, ind+m, jsph)
                        end do
                        tt = tt / tij
                    end do
                end do
                do ind = 1, nbasis
                    call intrhs(isph, vts(:, ind), mat(:, isph, ind, jsph))
                end do
            end do
        end do
        deallocate(vts, vplm, basloc, vcos, vsin, stat=istatus)
        if(istatus .ne. 0) then
            write(6,*) 'gen_mat: deallocation failed !'
            stop
        end if
    end subroutine gen_mat

    subroutine result_integrate(nbasis, ngrid, nsph, y, z)
        implicit none
        integer, intent(in) :: nbasis, ngrid, nsph
        real*8, intent(in) :: y(ngrid, nsph)
        real*8, intent(out) :: z(nbasis, nsph)
        integer :: isph
        z = zero
        do isph = 1, nsph
            call intrhs(isph, y(:, isph), z(:, isph))
        end do
    end subroutine result_integrate

    subroutine result_integrate_ui(nbasis, ngrid, nsph, y, z)
        implicit none
        integer, intent(in) :: nbasis, ngrid, nsph
        real*8, intent(in) :: y(ngrid, nsph)
        real*8, intent(out) :: z(nbasis, nsph)
        real*8 :: tmp(ngrid)
        integer :: isph
        z = zero
        do isph = 1, nsph
            tmp = y(:, isph) * ui(:, isph)
            call intrhs(isph, tmp, z(:, isph))
        end do
    end subroutine result_integrate_ui

    subroutine result_integrate_ui_grid_ext(ngrid_ext, x, nbasis, nsph, y)
        implicit none
        integer, intent(in) :: ngrid_ext, nbasis, nsph
        real*8, intent(in) :: x(ngrid_ext)
        real*8, intent(out) :: y(nbasis, nsph)
        integer :: isph, igrid, igrid_ext
        real*8 :: tmp(ngrid)
        igrid_ext = 0
        do isph = 1, nsph
            do igrid = 1, ngrid
                if(ui(igrid, isph) .eq. zero) then
                    tmp(igrid) = zero
                else
                    igrid_ext = igrid_ext + 1
                    tmp(igrid) = ui(igrid, isph) * x(igrid_ext)
                end if
            end do
            call intrhs(isph, tmp, y(:, isph))
        end do
    end subroutine result_integrate_ui_grid_ext

    subroutine get_ui(nsph, ngrid, ui_out)
        integer, intent(in) :: nsph, ngrid
        real(kind=8), intent(out) :: ui_out(ngrid, nsph)
        ui_out = ui
    end subroutine

end module mydx
