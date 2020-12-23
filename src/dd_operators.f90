!> @copyright (c) 2020-2020 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/dd_operators.f90
!! Operators of ddX software
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2020-12-17

!> Operators shared among ddX methods
module dd_operators
! Use underlying core routines
use dd_core

contains

!> Apply double layer potential to spherical harmonics
!!
!! @param[in] dd_data
subroutine dx(dd_data, do_diag, x, y)
    ! Inputs
    type(dd_data_type), intent(in) :: dd_data
    integer, intent(in) :: do_diag
    real(kind=rp), intent(in) :: x((dd_data % lmax+1)**2, dd_data % nsph)
    ! Output
    real(kind=rp), intent(out) :: y((dd_data % lmax+1)**2, dd_data % nsph)
    ! Local variables
    real(kind=rp), allocatable :: vts(:), vplm(:), basloc(:), vcos(:), vsin(:)
    real(kind=rp) :: c(3), vij(3), sij(3)
    real(kind=rp) :: vvij, tij, fourpi, tt, f, f1, rho, ctheta, stheta, cphi, sphi
    integer :: its, isph, jsph, l, m, ind, lm, istatus
    ! Allocate temporaries
    allocate(vts(dd_data % ngrid), vplm((dd_data % lmax+1)**2), &
        & basloc((dd_data % lmax+1)**2),vcos(dd_data % lmax + 1), &
        & vsin(dd_data % lmax + 1), stat=istatus)
    if (istatus.ne.0) then
        write(6,*) 'dx: allocation failed !'
        stop
    end if
    y = zero
    ! this loop is easily parallelizable
    ! !!$omp parallel do default(none) schedule(dynamic) &
    ! !!$omp private(isph,its,jsph,basloc,vplm,vcos,vsin,vij, &
    ! !!$omp vvij,tij,sij,tt,l,ind,f,m,vts,c) &
    ! !!$omp shared(nsph,ngrid,ui,csph,rsph,grid, &
    ! !!$omp lmax,fourpi,dodiag,x,y,basis)
    do isph = 1, dd_data % nsph
        ! compute the "potential" from the other spheres
        ! at the exposed lebedv points of the i-th sphere 
        vts = zero
        do its = 1, dd_data % ngrid
            if (dd_data % ui(its,isph).gt.zero) then
                c = dd_data % csph(:,isph) + dd_data % rsph(isph)* &
                    & dd_data % cgrid(:,its)
                do jsph = 1, dd_data % nsph
                    if (jsph.ne.isph) then
                        ! build the geometrical variables
                        vij = c - dd_data % csph(:,jsph)
                        vvij = sqrt(dot_product(vij,vij))
                        tij = vvij / dd_data % rsph(jsph)
                        sij = vij/vvij 
                        ! build the local basis
                        call ylmbas(sij, rho, ctheta, stheta, cphi, sphi, &
                            & dd_data % lmax, dd_data % vscales, basloc, &
                            & vplm, vcos, vsin)
                        ! with all the required stuff, finally compute
                        ! the "potential" at the point 
                        tt = one/tij 
                        do l = 0, dd_data % lmax
                            ind = l*l + l + 1
                            f = fourpi*dble(l)/(two*dble(l) + one)*tt
                            do m = -l, l
                                vts(its) = vts(its) + f*x(ind + m,jsph) * &
                                    & basloc(ind + m)
                            end do
                            tt = tt/tij
                        end do
                    else if (do_diag .eq. 1) then
                        ! add the diagonal contribution
                        do l = 0, dd_data % lmax
                            ind = l*l + l + 1
                            f = (two*dble(l) + one)/fourpi
                            do m = -l, l
                                vts(its) = vts(its) - pt5*x(ind + m,isph) * &
                                    & dd_data % vgrid(ind + m,its)/f
                            end do
                        end do
                    end if 
                end do
                vts(its) = dd_data % ui(its,isph)*vts(its) 
            end if
        end do
        ! now integrate the potential to get its modal representation
        call intrhs(dd_data % iprint, dd_data % ngrid, dd_data % lmax, &
            & dd_data % vwgrid, isph, vts, y(:,isph))
    end do
    ! Clean up temporary data
    deallocate(vts,vplm,basloc,vcos,vsin,stat=istatus)
    if (istatus.ne.0) then
        write(6,*) 'dx: deallocation failed !'
        stop
    end if
end subroutine dx

!> Apply adjoint double layer operator
!!
!!
subroutine dstarx(dd_data, do_diag, x, y)
    ! Inputs
    type(dd_data_type), intent(in) :: dd_data
    integer, intent(in) :: do_diag
    real(kind=rp), intent(in) :: x((dd_data % lmax+1)**2, dd_data % nsph)
    ! Output
    real(kind=rp), intent(out) :: y((dd_data % lmax+1)**2, dd_data % nsph)
    ! Local variables
    real(kind=rp), allocatable :: vts(:), vplm(:), basloc(:), vcos(:), vsin(:)
    real(kind=rp) :: c(3), vji(3), sji(3)
    real(kind=rp) :: vvji, tji, fourpi, tt, f, f1, rho, ctheta, stheta, cphi, sphi
    integer :: its, isph, jsph, l, m, ind, lm, istatus
    ! Allocate temporaries
    allocate(vts(dd_data % ngrid), vplm((dd_data % lmax+1)**2), &
        & basloc((dd_data % lmax+1)**2),vcos(dd_data % lmax + 1), &
        & vsin(dd_data % lmax + 1), stat=istatus)
    if (istatus.ne.0) then
        write(6,*) 'dx: allocation failed !'
        stop
    end if
    y = zero
    ! this loop is easily parallelizable
    ! !!$omp parallel do default(none) schedule(dynamic) &
    ! !!$omp private(isph,its,jsph,basloc,vplm,vcos,vsin,vji, &
    ! !!$omp vvji,tji,sji,tt,l,ind,f,m,vts,c) &
    ! !!$omp shared(nsph,ngrid,ui,csph,rsph,grid,facl, &
    ! !!$omp lmax,fourpi,dodiag,x,y,basis,w,nbasis)
    do isph = 1, dd_data % nsph
        do jsph = 1, dd_data % nsph
            if (jsph.ne.isph) then
                do its = 1, dd_data % ngrid
                    if (dd_data % ui(its,jsph).gt.zero) then
                        ! build the geometrical variables
                        vji = dd_data % csph(:,jsph) + dd_data % rsph(jsph) * &
                            & dd_data % cgrid(:,its) - dd_data % csph(:,isph)
                        vvji = sqrt(dot_product(vji,vji))
                        tji = vvji/dd_data % rsph(isph)
                        sji = vji/vvji
                        ! build the local basis
                        call ylmbas(sji, rho, ctheta, stheta, cphi, sphi, &
                            & dd_data % lmax, dd_data % vscales, basloc, &
                            & vplm, vcos, vsin)
                        tt = dd_data % ui(its,jsph)*dot_product(dd_data % vwgrid(:,its),x(:,jsph))/tji
                        do l = 0, dd_data % lmax
                            ind = l*l + l + 1
                            f = dble(l)*tt/ dd_data % vscales(ind)**2
                            do m = -l, l
                                y(ind+m,isph) = y(ind+m,isph) + f*basloc(ind+m)
                            end do
                            tt = tt/tji
                        end do
                    end if
                end do
            else if (do_diag .eq. 1) then
                do its = 1, dd_data % ngrid
                    f = pt5*dd_data % ui(its,jsph)*dot_product(dd_data % vwgrid(:,its),x(:,jsph))
                    do ind = 1, (dd_data % lmax+1)**2
                        y(ind,isph) = y(ind,isph) - f*dd_data % vgrid(ind,its)/dd_data % vscales(ind)**2
                    end do
                end do
            end if
        end do
    end do
    deallocate(vts,vplm,basloc,vcos,vsin,stat=istatus)
    if (istatus.ne.0) then
        write(6,*) 'dx: deallocation failed !'
        stop
    end if
end subroutine dstarx

!> Apply \f$ R_\varepsilon \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R_\varepsilon x = (2\pi(\varepsilon + 1) / (\varepsilon
!! - 1) - D) x \f$.
!!
!! @param[in] dd_data:
!! @param[in] do_diag:
!! @param[in] x:
!! @param[out] y:
subroutine repsx(dd_data, do_diag, x, y)
    ! Inputs
    type(dd_data_type), intent(in) :: dd_data
    integer, intent(in) :: do_diag
    real(kind=rp), intent(in) :: x((dd_data % lmax+1)**2, dd_data % nsph)
    ! Output
    real(kind=rp), intent(out) :: y((dd_data % lmax+1)**2, dd_data % nsph)
    ! Local variables
    real(kind=rp) :: fac
    ! Output `y` is cleaned here
    call dx(dd_data, do_diag, x, y)
    y = -y
    ! Apply diagonal
    if (do_diag .eq. 1) then
        fac = twopi * (dd_data % eps + one) / (dd_data % eps - one)
        y = y + fac*x
    end if
end subroutine repsx

!> Apply \f$ R_\infty \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R_\infty x = (2\pi - D) x \f$.
!!
!! @param[in] dd_data:
!! @param[in] do_diag:
!! @param[in] x:
!! @param[out] y:
subroutine rinfx(dd_data, do_diag, x, y)
    ! Inputs
    type(dd_data_type), intent(in) :: dd_data
    integer, intent(in) :: do_diag
    real(kind=rp), intent(in) :: x((dd_data % lmax+1)**2, dd_data % nsph)
    ! Output
    real(kind=rp), intent(out) :: y((dd_data % lmax+1)**2, dd_data % nsph)
    ! Local variables
    real(kind=rp) :: fac
    ! Output `y` is cleaned here
    call dx(dd_data, do_diag, x, y)
    y = -y
    ! Apply diagonal
    if (do_diag .eq. 1) then
        y = y + twopi*x
    end if
end subroutine rinfx

!> Apply adjoint f$ R_\varepsilon \f$ operator to spherical harmonics
!!
!! Compute \f$ y = R_\varepsilon^* x = (2\pi(\varepsilon + 1) / (\varepsilon
!! - 1) - D^*) x \f$.
!!
!! @param[in] dd_data:
!! @param[in] do_diag:
!! @param[in] x:
!! @param[out] y:
subroutine rstarepsx(dd_data, do_diag, x, y)
    ! Inputs
    type(dd_data_type), intent(in) :: dd_data
    integer, intent(in) :: do_diag
    real(kind=rp), intent(in) :: x((dd_data % lmax+1)**2, dd_data % nsph)
    ! Output
    real(kind=rp), intent(out) :: y((dd_data % lmax+1)**2, dd_data % nsph)
    ! Local variables
    real(kind=rp) :: fac
    ! Output `y` is cleaned here
    call dstarx(dd_data, do_diag, x, y)
    y = -y
    ! Apply diagonal
    if (do_diag .eq. 1) then
        fac = twopi * (dd_data % eps + one) / (dd_data % eps - one)
        y = y + fac*x
    end if
end subroutine rstarepsx

!> Apply preconditioner for 
subroutine apply_rstarx_prec(dd_data, x, y)
    ! Inputs
    type(dd_data_type), intent(in) :: dd_data
    real(kind=rp), intent(in) :: x((dd_data % lmax+1)**2, dd_data % nsph)
    real(kind=rp), intent(inout) :: y((dd_data % lmax+1)**2, dd_data % nsph)
    integer :: isph
    ! simply do a matrix-vector product with the transposed preconditioner 
    ! !!$omp parallel do default(shared) schedule(dynamic) &
    ! !!$omp private(isph)
    do isph = 1, dd_data % nsph
        call dgemm('t','n',(dd_data % lmax+1)**2,1,(dd_data % lmax+1)**2,one,&
            & dd_data % rx_prc(:,:,isph),(dd_data % lmax+1)**2, &
            & x(:,isph),(dd_data % lmax+1)**2,zero,y(:,isph),(dd_data % lmax+1)**2)
    end do
end subroutine apply_rstarx_prec

!> Compute preconditioner
!!
!! assemble the diagonal blocks of the reps matrix
!! then invert them to build the preconditioner
subroutine mkprec(dd_data)
    ! Inouts
    type(dd_data_type), intent(inout) :: dd_data
    integer :: isph, lm, ind, l1, m1, ind1, its, istatus
    real*8  :: f, f1
    integer, allocatable :: ipiv(:)
    real*8,  allocatable :: work(:)
    ! Allocation of temporaries
    allocate(ipiv(nbasis),work(nbasis*nbasis),stat=istatus)
    if (istatus.ne.0) then
        write(*,*) 'mkprec : allocation failed !'
        stop
    endif
    ! Init
    dd_data % rx_prc = zero
    ! Off-diagonal part
    do isph = 1, dd_data % nsph
        do its = 1, dd_data % ngrid
            f = twopi* dd_data % ui(its,isph) * dd_data % wgrid(its)
            do l1 = 0, dd_data % lmax
                ind1 = l1*l1 + l1 + 1
                do m1 = -l1, l1
                    f1 = f*dd_data % vgrid(ind1 + m1,its)/(two*dble(l1) + one)
                    do lm = 1, dd_data % nbasis
                        dd_data % rx_prc(lm,ind1 + m1,isph) = dd_data % rx_prc(lm,ind1 + m1,isph) + &
                            & f1*dd_data % vgrid(lm,its)
                    end do
                end do
            end do
        end do
    end do
    ! add diagonal
    f = twopi*(dd_data % eps + one)/(dd_data % eps - one)
    do isph = 1, dd_data % nsph
        do lm = 1, dd_data % nbasis
            dd_data % rx_prc(lm,lm,isph) = dd_data % rx_prc(lm,lm,isph) + f
        end do
    end do
    ! invert the blocks
    do isph = 1, dd_data % nsph
        call dgetrf(dd_data % nbasis, dd_data % nbasis, &
            & dd_data % rx_prc(:,:,isph), dd_data % nbasis, ipiv, istatus)
        if (istatus.ne.0) then 
            write(6,*) 'lu failed in mkprc'
            stop
        end if
        call dgetri(dd_data % nbasis, dd_data % rx_prc(:,:,isph), &
            & dd_data % nbasis, ipiv, work, dd_data % nbasis**2, istatus)
        if (istatus.ne.0) then 
            write(6,*) 'inversion failed in mkprc'
            stop
        end if
    end do
    ! Cleanup temporaries
    deallocate(work, ipiv, stat=istatus)
    if (istatus.ne.0) then
        write(*,*) 'mkprec : deallocation failed !'
        stop
    end if
endsubroutine mkprec

end module dd_operators
