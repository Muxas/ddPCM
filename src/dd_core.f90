!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/dd_core.f90
!! Core routines and parameters of entire ddX software
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2020-12-17

!> Core routines and parameters of ddX software
module dd_core
implicit none

!> Precision of real numbers
!!
!! This parameter chooses single or double precision or maybe even quad
integer, parameter :: rp=selected_real_kind(15)

!! Compile-time constants
real(kind=rp), parameter :: zero = 0d0, one = 1d0, two = 2d0, three = 3d0, four = 4d0
real(kind=rp), parameter :: pt5 = 5d-1
real(kind=rp), parameter :: sqrt2 = sqrt(two)
real(kind=rp), parameter :: pi4 = atan(one)
real(kind=rp), parameter :: pi = four * pi4
real(kind=rp), parameter :: fourpi = four * pi
real(kind=rp), parameter :: twopi = two * pi
real(kind=rp), parameter :: sqrt4pi = four * sqrt(atan(one))
!> Number of supported Lebedev grids
integer, parameter :: nllg = 32
!> Number of grid points of each Lebedev grid
integer, parameter :: ng0(nllg) = (/ 6, 14, 26, 38, 50, 74, 86, 110, 146, &
    & 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, &
    & 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 /)
!> Shift of characteristic function
! AJ: Difference from ddLPB Changed here
real(kind=rp), parameter :: se = -1.0d0

!> Main ddX type that stores all required information
type dd_data_type
    !!!!!!!!!!!! Parameters
    !> Model to use 1 for cosmo, 2 for pcm, 3 for lpb
    integer :: model
    !> Printing flag
    integer :: iprint
    !> Number of atoms with their centers and Van der Waals radii
    integer :: nsph
    !> Centers of atoms of dimension (3, nsph)
    real(kind=rp), allocatable :: csph(:, :)
    !> Array of radii of atoms of dimension (nsph)
    real(kind=rp), allocatable :: rsph(:)
    !> Dielectric permittivity inside
    real(kind=rp) :: epsin
    !> Dielectric permittivity outside
    real(kind=rp) :: epsout
    !> Relative dielectric permittivity
    real(kind=rp) :: eps
    !> Debye Hückel parameter
    real(kind=rp) :: kappa
    !> Regularization parameter
    real(kind=rp) :: eta
    !> Relative error threshold for iterative solvers
    real(kind=rp) :: tol
    !> Maximal degree of modeling spherical harmonics
    integer :: lmax
    !> Number modeling spherical harmonics per sphere
    integer :: nbasis
    !> Total number of modeling spherical harmonics
    integer :: n
    !> Number of Lebedev grid points on each sphere
    integer :: ngrid
    !> Whether to compute (1) or not (0) forces
    integer :: force
    !> Enable (1) or disable (0) use of FMM techniques
    integer :: fmm
    !!!!!!!!!!!!!! Constants, initialized by ddinit
    !> Maximal degree of used real normalized spherical harmonics
    !!
    !! For example, if FMM is
    !! used, its M2L operation requires computing spherical harmonics of all
    !! degrees up to `pm+pl`. If `force=1` then this parameter might be
    !! increased by 1 because certain implementations of analytical gradients
    !! might rely on spherical harmonics of +1 degree.
    integer :: dmax
    !> Scales of real normalized spherical harmonics of degree up to dmax
    !!
    !! This array has `(dmax+1)**2` entries
    real(kind=rp), allocatable :: vscales(:)
    !> Maximum sqrt of factorial precomputed
    !!
    !! Just like with `dmax` parameter, number of used factorials is either
    !! `2*lmax+1` or `2*(pm+pl)+1` depending on whether FMM is used or not
    integer :: nfact
    !> Array of square roots of factorials of dimension (nfact)
    real(kind=rp), allocatable :: vfact(:)
    !> Coordinates of Lebedev quadrature points of dimension (3, ngrid)
    real(kind=rp), allocatable :: cgrid(:, :)
    !> Weights of Lebedev quadrature points of dimension (ngrid)
    real(kind=rp), allocatable :: wgrid(:)
    !> Maximal degree of spherical harmonics evaluated at Lebedev grid points
    !!
    !! Although we use spherical harmonics of degree up to `dmax`, only
    !! spherical harmonics of degree up to `lmax` and `pl` are evaluated
    !! at Lebedev grid points. In the case `force=1` this degrees might be
    !! increased by 1 depending on implementation of gradients.
    integer :: vgrid_dmax
    !> Number of spherical harmonics evaluated at Lebedev grid points
    integer :: vgrid_nbasis
    !> Values of spherical harmonics at Lebedev grid points
    !!
    !! Dimensions of this array are (vgrid_nbasis, ngrid)
    real(kind=rp), allocatable :: vgrid(:, :)
    !> Weighted values of spherical harmonics at Lebedev grid points
    !!
    !! vwgrid(:, igrid) = vgrid(:, igrid) * wgrid(igrid)
    !! Dimensions of this array are ((vgrid_nbasis, ngrid)
    real(kind=rp), allocatable :: vwgrid(:, :)
    !> Maximum number of neighbours per sphere. Input from a user
    integer :: nngmax
    !> List of intersecting spheres in a CSR format
    integer, allocatable :: inl(:)
    !> List of intersecting spheres in a CSR format
    integer, allocatable :: nl(:)
    !> Characteristic function f. Dimension is (ngrid, npsh)
    real(kind=rp), allocatable :: fi(:, :)
    !> Characteristic function U. Dimension is (ngrid, nsph)
    real(kind=rp), allocatable :: ui(:, :)
    !> Derivative of characteristic function U. Dimension is (3, ngrid, nsph)
    real(kind=rp), allocatable :: zi(:, :, :)
    !> Number of external Lebedev grid points on a molecule surface
    integer :: ncav
    !> Coordinates of external Lebedev grid points of dimension (3, ncav)
    real(kind=rp), allocatable :: ccav(:, :)
    !> Preconditioner for operator R_eps. Dimension is (nbasis, nbasis, nsph)
    real(kind=rp), allocatable :: rx_prc(:, :, :)
    ! Maximum degree of spherical harmonics for M (multipole) expansion
    integer :: pm
    ! Maximum degree of spherical harmonics for L (local) expansion
    integer :: pl
    !! Cluster tree information
    ! Reordering of spheres for better locality
    integer, allocatable :: ind(:)
    ! First and last spheres, corresponding to given node of cluster tree
    integer, allocatable :: cluster(:, :)
    ! Children of each cluster
    integer, allocatable :: children(:, :)
    ! Parent of each cluster
    integer, allocatable :: parent(:)
    ! Center of bounding sphere of each cluster
    real(kind=rp), allocatable :: cnode(:, :)
    ! Radius of bounding sphere of each cluster
    real(kind=rp), allocatable :: rnode(:)
    ! Which leaf node contains only given input sphere
    integer, allocatable :: snode(:)
    ! Number of nodes (subclusters)
    integer :: nclusters
    ! Height of a tree
    integer :: height
    ! Total number of far and near admissible pairs
    integer :: nnfar, nnnear
    ! Number of admissible pairs for each node
    integer, allocatable :: nfar(:), nnear(:)
    ! Arrays of admissible far and near pairs
    integer, allocatable :: far(:), near(:)
    ! Index of first element of array of admissible far and near pairs
    integer, allocatable :: sfar(:), snear(:)
    ! External grid points
    integer :: ngrid_ext, ngrid_ext_near
    integer, allocatable :: ngrid_ext_sph(:), grid_ext_ia(:), grid_ext_ja(:)
    ! Reflection matrices for M2M, M2L and L2L operations
    real(kind=8), allocatable :: m2m_reflect_mat(:, :), m2m_ztrans_mat(:, :)
    real(kind=8), allocatable :: l2l_reflect_mat(:, :), l2l_ztrans_mat(:, :)
    real(kind=8), allocatable :: m2l_reflect_mat(:, :), m2l_ztrans_mat(:, :)
    ! Far-field L2P and near-field M2P
    real(kind=8), allocatable :: l2p_mat(:, :), m2p_mat(:, :)

end type dd_data_type

contains

!> Initialize ddX input with a full set of parameters
!!
!! @param[in] n: Number of atoms
!! @param[in] x: \f$ x \f$ coordinates of atoms. Dimension is `(n)`
!! @param[in] y: \f$ y \f$ coordinates of atoms. Dimension is `(n)`
!! @param[in] z: \f$ z \f$ coordinates of atoms. Dimension is `(n)`
!! @param[in] rvdw: Van-der-Waals radii of atoms. Dimension is `(n)`
!! @param[in] model: Choose model: 1 for COSMO, 2 for PCM and 3 for LPB
!! @param[in] lmax: Maximal degree of modeling spherical harmonics. `lmax` >= 0
!! @param[inout] ngrid: Approximate number of Lebedev grid points on input and
!!      actual number of grid points on exit. `ngrid` > 0
!! @param[in] force: 1 if forces are required and 0 otherwise
!! @param[in] fmm: 1 to use FMM acceleration and 0 otherwise
!! @param[in] pm: Maximal degree of multipole spherical harmonics. `pm` >= 0
!! @param[in] pl: Maximal degree of local spherical harmonics. `pl` >= 0
!! @param[in] iprint: Level of printing debug info
!! @param[in] nngmax: Maximal number of neighbours (intersects) for every
!!      sphere
!! @param[in] eta: regularization parameter
!! @param[in] eps: Relative dielectric permittivity
!! @param[out] dd_data: Object containing all inputs
!! @param[out] info: flag of succesfull exit
!!      = 0: Succesfull exit
!!      < 0: If info=-i then i-th argument had an illegal value
!!      > 0: Allocation of a buffer for the output dd_data failed
!!
!!
!! TODO: intersecting spheres must take into account regularization parameter
subroutine ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, pl, &
        & iprint, nngmax, eta, eps, kappa, dd_data)!, info)
    ! Inputs
    integer, intent(in) :: n, model, lmax, force, fmm, pm, pl, iprint, nngmax
    real(kind=rp) :: x(n), y(n), z(n), rvdw(n), eta, eps, kappa
    ! Output
    type(dd_data_type), intent(out) :: dd_data
    integer :: info
    ! Inouts
    integer, intent(inout) :: ngrid
    ! Local variables
    integer :: istatus, i, ii, inear, jnear, igrid, jgrid, isph, jsph, lnl
    real(kind=rp) :: rho, ctheta, stheta, cphi, sphi
    real(kind=rp), allocatable :: vplm(:), vcos(:), vsin(:)
    real(kind=rp) :: v(3), vv, t, swthr, fac, d2, r2
    double precision, external :: dnrm2
    !!! Check input parameters
    if (n .lt. 0) then
        !write(*, *) "ddinit: wrong value of parameter `n`"
        info = -1
        return
    end if
    dd_data % nsph = n
    allocate(dd_data % csph(3, n), dd_data % rsph(n), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "ddinit: [1] allocation failed!"
        info = 1
        return
    end if
    dd_data % csph(1, :) = x
    dd_data % csph(2, :) = y
    dd_data % csph(3, :) = z
    dd_data % rsph = rvdw
    if ((model .lt. 1) .and. (model .gt. 3)) then
        !write(*, *) "ddinit: wrong value of parameter `model`"
        info = -6
        return
    end if
    dd_data % model = model
    if (lmax .lt. 0) then
        !write(*, *) "ddinit: wrong value of parameter `lmax`"
        info = -7
        return
    end if
    dd_data % lmax = lmax
    dd_data % nbasis = (lmax+1)**2
    dd_data % n = n * dd_data % nbasis
    if (ngrid .lt. 0) then
        !write(*, *) "ddinit: wrong value of parameter `ngrid`"
        info = -8
        return
    end if
    ! Actual value of ngrid will be calculated later
    if ((force .lt. 0) .and. (force .gt. 1)) then
        !write(*, *) "ddinit: wrong value of parameter `force`"
        info = -9
        return
    end if
    dd_data % force = force
    if ((fmm .lt. 0) .and. (fmm .gt. 1)) then
        !write(*, *) "ddinit: wrong value of parameter `fmm`"
        info = -10
        return
    end if
    dd_data % fmm = fmm
    if (iprint .lt. 0) then
        !write(*, *) "ddinit: wrong value of parameter `iprint`"
        info = -13
        return
    end if
    dd_data % iprint = iprint
    if (nngmax .le. 0) then
        !write(*, *) "ddinit: wrong value of parameter `nngmax`"
        info = -14
        return
    end if
    dd_data % nngmax = nngmax
    if ((eta .lt. 0) .or. (eta .gt. 1)) then
        !write(*, *) "ddinit: wrong value of parameter `eta`"
        info = -15
        return
    end if
    dd_data % eta = eta
    if (eps .lt. 0) then
        !write(*, *) "ddinit: wrong value of parameter `eps`"
        info = -16
        return
    end if
    dd_data % eps = eps
    if (kappa .lt. 0) then
        !write(*, *) "ddinit: wrong value of parameter `kappa`"
        info = -16
        return
    end if
    dd_data % kappa = kappa
    ! Set FMM parameters
    if(fmm .eq. 1) then
        if(pm .lt. 0) then
            write(*, *) "ddinit: wrong value of parameter `pm`"
            info = -11
            return
        end if
        if(pl .lt. 0) then
            write(*, *) "ddinit: wrong value of parameter `pl`"
            info = -12
            return
        end if
    end if
    dd_data % pm = pm
    dd_data % pl = pl
    !!! Generate constants and auxiliary arrays
    ! Compute sizes of auxiliary arrays for `fmm=0`
    if (fmm .eq. 0) then
        dd_data % dmax = lmax
        dd_data % vgrid_dmax = lmax
    ! Compute sizes of arrays if `fmm=1`
    else
        ! If forces are required then we need M2P of degree lmax+1 for
        ! near-field analytical gradients
        if (force .eq. 1) then
            dd_data % dmax = max(pm+pl, lmax+1)
        else
            dd_data % dmax = max(pm+pl, lmax)
        end if
        dd_data % vgrid_dmax = max(pl, lmax)
    end if
    dd_data % vgrid_nbasis = (dd_data % vgrid_dmax+1)**2
    dd_data % nfact = max(2*dd_data % dmax+1, 2)
    ! Compute scaling factors of spherical harmonics
    allocate(dd_data % vscales(dd_data % nbasis), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "ddinit: [2] allocation failed!"
        info = 2
        return
    end if
    call ylmscale(dd_data % dmax, dd_data % vscales)
    ! Precompute square roots of factorials
    allocate(dd_data % vfact(dd_data % nfact), stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) "ddinit: [3] allocation failed!"
        stop
    end if
    dd_data % vfact(1) = 1
    do i = 2, dd_data % nfact
        dd_data % vfact(i) = dd_data % vfact(i-1) * sqrt(dble(i-1))
    end do
    ! Get nearest number of Lebedev grid points
    igrid = 0
    inear = 100000
    do i = 1, nllg
        jnear = iabs(ng0(i)-ngrid)
        if (jnear.lt.inear) then
            inear = jnear
            igrid = i
        end if
    end do
    ngrid = ng0(igrid)
    dd_data % ngrid = ngrid
    ! Get weights and coordinates of Lebedev grid points
    allocate(dd_data % cgrid(3, ngrid), dd_data % wgrid(ngrid), stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) "ddinit: [4] allocation failed!"
        stop
    end if
    call llgrid(ngrid, dd_data % wgrid, dd_data % cgrid)
    ! Compute spherical harmonics at grid points
    allocate(dd_data % vgrid(dd_data % vgrid_nbasis, ngrid), &
        & dd_data % vwgrid(dd_data % vgrid_nbasis, ngrid), stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) "ddinit: [5] allocation failed!"
        stop
    end if
    allocate(vplm(dd_data % vgrid_nbasis), vcos(dd_data % vgrid_dmax+1), &
        & vsin(dd_data % vgrid_dmax+1), stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) "ddinit: [6] allocation failed!"
        stop
    end if
    do igrid = 1, ngrid
        call ylmbas(dd_data % cgrid(:, igrid), rho, ctheta, stheta, cphi, &
            & sphi, dd_data % vgrid_dmax, dd_data % vscales, &
            & dd_data % vgrid(:, igrid), vplm, vcos, vsin)
        dd_data % vwgrid(:, igrid) = dd_data % wgrid(igrid) * &
            & dd_data % vgrid(:, igrid)
    end do
    ! Debug printing to compare against ddPCM
    if (iprint.ge.4) then
        call prtsph('facs', dd_data % nbasis, lmax, 1, 0, dd_data % vscales)
        ! TODO: this printing is obviously wrong, but it is the same in ddPCM
        call prtsph('facl', dd_data % nbasis, lmax, 1, 0, dd_data % vscales)
        call prtsph('basis', dd_data % nbasis, lmax, ngrid, 0, dd_data % vgrid)
        call ptcart('grid', ngrid, 3, 0, dd_data % cgrid)
        call ptcart('weights', ngrid, 1, 0, dd_data % wgrid)
    end if
    deallocate(vplm, vcos, vsin, stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) "ddinit: [1] deallocation failed!"
        stop
    end if
    ! Upper bound of switch region. Defines intersection criterion for spheres
    swthr = one + (se+one)*eta/two
    ! Build list of neighbours in CSR format
    allocate(dd_data % inl(n+1), dd_data % nl(n*nngmax), stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) 'ddinit : [7] allocation failed !'
        stop
    end if
    i  = 1
    lnl = 0
    do isph = 1, n
        dd_data % inl(isph) = lnl + 1
        do jsph = 1, dd_data % nsph
            if (isph .ne. jsph) then
                d2 = dnrm2(3, &
                    & dd_data % csph(:, isph)-dd_data % csph(:, jsph), 1)
                r2 = rvdw(isph) + rvdw(jsph)
                ! This check does not take into account eta regularization
                ! TODO: take eta into account
                ! Solution is in the following commented line
                ! r2 = rvdw(isph) + swthr*rvdw(jsph)
                if (d2 .le. r2) then
                    dd_data % nl(i) = jsph
                    i  = i + 1
                    lnl = lnl + 1
                end if
            end if
        end do
    end do
    dd_data % inl(n+1) = lnl+1
    ! Some format data that I will throw away as soon as this code works
    1000 format(t3,'neighbours of sphere ',i6)
    1010 format(t5,12i6)
    ! Debug printing
    if (iprint.ge.4) then
        write(6,*) '   inl:'
        write(6,'(10i6)') dd_data % inl(1:n+1)
        write(6,*)
        do isph = 1, n
            write(6,1000) isph
            write(6,1010) dd_data % nl(dd_data % inl(isph): &
                & dd_data % inl(isph+1)-1)
        end do
        write(6,*)
    end if
    ! Build arrays fi, ui, zi
    allocate(dd_data % fi(ngrid, n), dd_data % ui(ngrid, n), stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) 'ddinit : [8] allocation failed !'
        stop
    end if
    dd_data % fi = zero
    dd_data % ui = zero
    if (force .eq. 1) then
        allocate(dd_data % zi(3, ngrid, n), stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) 'ddinit : [9] allocation failed !'
            stop
        endif
        dd_data % zi = zero
    end if
    do isph = 1, n
        do i = 1, ngrid
            ! Loop over neighbours of i-th sphere
            do ii = dd_data % inl(isph), dd_data % inl(isph+1)-1
                ! Neighbour's index
                jsph = dd_data % nl(ii)
                ! Compute t_n^ij
                v(:) = dd_data % csph(:, isph) - dd_data % csph(:, jsph) + &
                    & rvdw(isph)*dd_data % cgrid(:, i)
                vv = dnrm2(3, v, 1)
                t = vv / dd_data % rsph(jsph)
                ! Accumulate characteristic function \chi(t_n^ij)
                dd_data % fi(i, isph) = dd_data % fi(i, isph) + fsw(t, se, eta)
                ! Check if gradients are required
                if (force .eq. 1) then
                    ! Check if t_n^ij belongs to switch region
                    if ((t .lt. swthr) .and. (t .gt. swthr-eta)) then
                        ! Accumulate gradient of characteristic function \chi
                        fac = dfsw(t, se, eta) / rvdw(jsph)
                        dd_data % zi(:, i, isph) = dd_data % zi(:, i, isph) + &
                            & fac*v/vv
                    end if
                end if
            enddo
            ! Compute characteristic function of a molecular surface ui
            if (dd_data % fi(i, isph) .le. one) then
                dd_data % ui(i, isph) = one - dd_data % fi(i, isph)
            end if
        enddo
    enddo
    ! Debug printing
    if (iprint .ge. 4) then
        call ptcart('fi', ngrid, n, 0, dd_data % fi)
        call ptcart('ui', ngrid, n, 0, dd_data % ui)
    end if
    ! Build cavity array. At first get total count
    dd_data % ncav = 0
    do isph = 1, n
        ! Loop over integration points
        do i = 1, ngrid
            ! Positive contribution from integration point
            if (dd_data % ui(i, isph) .gt. zero) then
                dd_data % ncav = dd_data % ncav + 1
            endif
        enddo
    enddo
    ! Allocate cavity array
    allocate(dd_data % ccav(3, dd_data % ncav) , stat=istatus)
    if (istatus .ne. 0) then
        write(*,*)'ddinit : [10] allocation failed!'
        stop
    endif
    ! Get actual cavity coordinates
    ii = 0
    do isph = 1, n
        ! Loop over integration points
        do i = 1, ngrid
            ! Positive contribution from integration point
            if (dd_data % ui(i, isph) .gt. zero) then
                ! Advance cavity array index
                ii = ii + 1
                ! Store point
                dd_data % ccav(:, ii) = dd_data % csph(:, isph) + &
                    & rvdw(isph)*dd_data % cgrid(:, i)
            endif
        enddo
    enddo
    ! Create preconditioner for PCM
    if (model .eq. 2) then
        allocate(dd_data % rx_prc(dd_data % nbasis, dd_data % nbasis, n), &
            & stat=istatus)
        if (istatus .ne. 0) then
            write(*,*)'ddinit : [11] allocation failed!'
            stop
        endif
        call mkprec(dd_data)
    end if
    ! Again some debug output
    1100  format(t3,i8,3f14.6)
    if (iprint .ge. 4) then
        write(6, *) '   external cavity points:'
        do ii = 1, dd_data % ncav
            write(6,1100) ii, dd_data % ccav(:, ii)
        end do
        write(6, *)
    end if
end subroutine ddinit

!> Deallocate object with corresponding data
!!
!! @param[inout] dd_data: object to deallocate
subroutine ddfree(dd_data)
    ! Input/output
    type(dd_data_type), intent(inout) :: dd_data
    ! Local variables
    integer :: istatus
    deallocate(dd_data % csph, dd_data % rsph, stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) "ddfree: [1] deallocation failed!"
        stop
    end if
    deallocate(dd_data % vscales, dd_data % vfact, stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) "ddfree: [2] deallocation failed!"
        stop
    end if
    deallocate(dd_data % cgrid, dd_data % wgrid, stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) "ddfree: [3] deallocation failed!"
        stop
    end if
    deallocate(dd_data % vgrid, dd_data % vwgrid, stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) "ddfree: [4] deallocation failed!"
        stop
    end if
    deallocate(dd_data % inl, dd_data % nl, stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) "ddfree: [5] deallocation failed!"
        stop
    end if
    deallocate(dd_data % fi, dd_data % ui, stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) "ddfree: [6] deallocation failed!"
        stop
    end if
    if (dd_data % force .eq. 1) then
        deallocate(dd_data % zi, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [7] deallocation failed!"
            stop
        end if
    end if
    deallocate(dd_data % ccav, stat=istatus)
    if (istatus .ne. 0) then
        write(*, *) "ddfree: [8] deallocation failed!"
        stop
    end if
    if (dd_data % model .eq. 2) then
        deallocate(dd_data % rx_prc, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [9] deallocation failed!"
            stop
        end if
    end if
end subroutine ddfree

!> Print array of spherical harmonics
!!
!! Prints (nbasis, ncol) array
!!
!! @param[in] label: Label to print
!! @param[in] nbasis: Number of rows of input x. nbasis >= 0
!! @param[in] lmax: Maximal degree of corresponding spherical harmonics.
!!      (lmax+1)**2 = nbasis
!! @param[in] ncol: Number of columns to print
!! @param[in] icol: This number is only for printing purposes
!! @param[in] x: Actual data to print
subroutine prtsph(label, nbasis, lmax, ncol, icol, x)
    ! Inputs
    character (len=*), intent(in) :: label
    integer, intent(in) :: nbasis, lmax, ncol, icol
    real(kind=rp), intent(in) :: x(nbasis, ncol)
    ! Local variables
    integer :: l, m, ind, noff, nprt, ic, j
    ! Print header:
    if (ncol .eq. 1) then
        write (6,'(3x,a,1x,"(column ",i4")")') label, icol
    else
        write (6,'(3x,a)') label
    endif
    ! Print entries:
    if (ncol .eq. 1) then
        do l = 0, lmax
            ind = l*l + l + 1
            do m = -l, l
                write(6,1000) l, m, x(ind+m, 1)
            end do
        end do
    else
        noff = mod(ncol, 5)
        nprt = max(ncol-noff, 0)
        do ic = 1, nprt, 5
            write(6,1010) (j, j = ic, ic+4)
            do l = 0, lmax
                ind = l*l + l + 1
                do m = -l, l
                    write(6,1020) l, m, x(ind+m, ic:ic+4)
                end do
            end do
        end do
        write (6,1010) (j, j = nprt+1, nprt+noff)
        do l = 0, lmax
            ind = l*l + l + 1
            do m = -l, l
                write(6,1020) l, m, x(ind+m, nprt+1:nprt+noff)
            end do
        end do
    end if
    1000 format(1x,i3,i4,f14.8)
    1010 format(8x,5i14)
    1020 format(1x,i3,i4,5f14.8)
end subroutine prtsph

!> Print array of quadrature points
!!
!! Prints (ngrid, ncol) array
!!
!! @param[in] label: Label to print
!! @param[in] ngrid: Number of rows of input x. ngrid >= 0
!! @param[in] ncol: Number of columns to print
!! @param[in] icol: This number is only for printing purposes
!! @param[in] x: Actual data to print
subroutine ptcart(label, ngrid, ncol, icol, x)
    ! Inputs
    character (len=*), intent(in) :: label
    integer, intent(in) :: ngrid, ncol, icol
    real(kind=rp), intent(in) :: x(ngrid, ncol)
    ! Local variables
    integer :: ig, noff, nprt, ic, j
    ! Print header :
    if (ncol .eq. 1) then
        write (6,'(3x,a,1x,"(column ",i4")")') label, icol
    else
        write (6,'(3x,a)') label
    endif
    ! Print entries :
    if (ncol .eq. 1) then
        do ig = 1, ngrid
            write(6,1000) ig, x(ig, 1)
        enddo
    else
        noff = mod(ncol, 5)
        nprt = max(ncol-noff, 0)
        do ic = 1, nprt, 5
            write(6,1010) (j, j = ic, ic+4)
            do ig = 1, ngrid
                write(6,1020) ig, x(ig, ic:ic+4)
            end do
        end do
        write (6,1010) (j, j = nprt+1, nprt+noff)
        do ig = 1, ngrid
            write(6,1020) ig, x(ig, nprt+1:nprt+noff)
        end do
    end if
    !
    1000 format(1x,i5,f14.8)
    1010 format(6x,5i14)
    1020 format(1x,i5,5f14.8)
    !
end subroutine ptcart

!> Print dd Solution vector
!!
!! @param[in] label : Label to print
!! @param[in] vector: Vector to print

subroutine print_ddvector(dd_data, label, vector)
    implicit none
    type(dd_data_type), intent(in)  :: dd_data
    character(len=*) :: label
    real(kind=rp) :: vector(dd_data % nbasis, dd_data % nsph)
    integer :: isph, lm
  
    write(6,*) label
    do isph = 1, dd_data % nsph
      do lm = 1, dd_data % nbasis
        write(6,'(F15.8)') vector(lm,isph)
      end do
    end do
    return
end subroutine print_ddvector


!> Compute scaling factors of real normalized spherical harmonics
!!
!! Output values of scaling factors of \f$ Y_\ell^m \f$ harmonics are filled
!! only for non-negative \f$ m \f$ since scaling factor of \f$ Y_\ell^{-m} \f$
!! is the same as scaling factor of \f$ Y_\ell^m \f$.
!!
!! @param[in] p: Maximal degree of spherical harmonics. `p` >= 0
!! @param[out] vscales: Array of scaling factors. Dimension is `(p+1)**2`
subroutine ylmscale(p, vscales)
    ! Input
    integer, intent(in) :: p
    ! Output
    real(kind=rp), intent(out) :: vscales((p+1)**2)
    ! Local variables
    real(kind=rp) :: tmp
    integer :: l, ind, m
    do l = 0, p
        ! m = 0
        ind = l*l + l + 1
        vscales(ind) = sqrt(dble(2*l+1)) / sqrt4pi
        tmp = vscales(ind) * sqrt2
        ! m != 0
        do m = 1, l
            tmp = -tmp / sqrt(dble((l-m+1)*(l+m)))
            ! Fill only positive m
            vscales(ind+m) = tmp
        end do
    end do
end subroutine ylmscale

!> Compute arrays of \f$ \cos(m \phi) \f$ and \f$ \sin(m \phi) \f$
!!
!! All values are computed recurrently from input \f$ \cos(\phi) \f$ and \f$
!! \sin(\phi) \f$ without accessing arccos or arcsin functions.
!!
!! @param[in] cphi: \f$ \cos(\phi) \f$. -1 <= `cphi` <= 1
!! @param[in] sphi: \f$ \sin(\phi) \f$. -1 <= `sphi` <= 1
!! @param[in] p: Maximal value of \f$ m \f$, for which to compute \f$ \cos(m
!!      \phi) \f$ and \f$ \sin(m\phi) \f$. `p` >= 0
!! @param[out] vcos: Array of \f$ \cos(m\phi) \f$ for \f$ m=0..p \f$. Dimension
!!      is `(p+1)`
!! @param[out] vsin: Array of \f$ \sin(m\phi) \f$ for \f$ m=0..p \f$. Dimension
!!      is `(p+1)`
subroutine trgev(cphi, sphi, p, vcos, vsin)
    ! Inputs
    real(kind=rp), intent(in) :: cphi, sphi
    integer, intent(in) :: p
    ! Output
    real(kind=rp), intent(out) :: vcos(p+1), vsin(p+1)
    ! Local variables
    integer :: m
    ! Set cos(0) and sin(0)
    vcos(1) = one
    vsin(1) = zero
    ! Define cos(m*phi) and sin(m*phi) recurrently
    do m = 2, p+1
        vcos(m) = vcos(m-1)*cphi - vsin(m-1)*sphi
        vsin(m) = vcos(m-1)*sphi + vsin(m-1)*cphi
    end do
end subroutine trgev

!> Compute associated Legendre polynomials
!!
!! Only polynomials \f$ P_\ell^m (\cos \theta) \f$ with non-negative parameter
!! \f$ m \f$ are computed. Implemented via following recurrent formulas:
!! \f{align}{
!!      &P_0^0(\cos \theta) = 1\\
!!      &P_{m+1}^{m+1}(\cos \theta) = -(2m+1) \sin \theta P_m^m(\cos \theta) \\
!!      &P_{m+1}^m(\cos \theta) = \cos \theta (2m+1) P_m^m(\cos \theta) \\
!!      &P_\ell^m(\cos \theta) = \frac{1}{\ell-m} \left( \cos \theta (2\ell-1)
!!      P_{\ell-1}^m(\cos \theta) - (\ell+m-1)P_{\ell-2}^m(\cos \theta)
!!      \right), \quad \forall \ell \geq m+2.
!! \f}
!!
!! @param[in] ctheta: \f$ \cos(\theta) \f$. -1 <= `ctheta` <= 1
!! @param[in] stheta: \f$ \sin(\theta) \f$. 0 <= `stheta` <= 1
!! @param[in] p: Maximal degree of polynomials to compute. `p` >= 0
!! @param[out] vplm: Values of associated Legendre polynomials. Dimension is
!!      `(p+1)**2`
subroutine polleg(ctheta, stheta, p, vplm)
    ! Inputs
    real(kind=rp), intent(in) :: ctheta, stheta
    integer, intent(in) :: p
    ! Outputs
    real(kind=rp), intent(out) :: vplm((p+1)**2)
    ! Local variables
    integer :: m, ind, l, ind2
    real(kind=rp) :: fact, pmm, pmm1, pmmo, pll, fm, fl
    ! Init aux factors
    fact = one
    pmm = one
    ! This loop goes over non-negative upper index of P_l^m, namely m, and
    ! defines value for all l from m to p. Polynomials P_l^m are defined only
    ! for l >= |m|. Only positive values of m are filled. Here we
    ! define at first P_m^m, then P_{m+1}^m and then all remaining P_l^m for
    ! l from m+2 to p.
    do m = 0, p
        ! index of P_m^m
        ind = (m+1) * (m+1)
        ! Store P_m^m
        vplm(ind) = pmm
        if (m .eq. p) then
            return
        end if
        fm = dble(m)
        ! index of P_{m+1}^m
        ind2 = ind + 2*m + 2
        ! P_{m+1}^m
        pmm1 = ctheta * (2*fm+1) * pmm
        vplm(ind2) = pmm1
        ! Save value P_m^m for recursion
        pmmo = pmm
        ! Fill values of P_l^m
        do l = m+2, p
            ! pmm corresponds to P_{l-2}^m
            ! pmm1 corresponds to P_{l-1}^m
            fl = dble(l)
            ! Compute P_l^m
            pll = ctheta*(2*fl-1)*pmm1 - (fl+fm-1)*pmm
            pll = pll / (fl-fm)
            ! Store P_l^m
            vplm(l*l + l + m + 1) = pll
            ! Save P_{l-1}^m and P_l^m for recursion
            pmm = pmm1
            pmm1 = pll
        end do
        ! Value of P_{m+1}^{m+1}
        pmm = -pmmo * fact * stheta
        fact = fact + two
    end do
end subroutine polleg

!> Compute all spherical harmonics up to a given degree at a given point
!!
!! Spherical harmonics are computed for a point \f$ x / \|x\| \f$. Cartesian
!! coordinate of input `x` is translated into a spherical coordinate \f$ (\rho,
!! \theta, \phi) \f$ that is represented by \f$ \rho, \cos \theta, \sin \theta,
!! \cos \phi \f$ and \f$ \sin \phi \f$. If \f$ \rho=0 \f$ nothing is computed,
!! only zero \f$ \rho \f$ is returned without doing anything else. If \f$
!! \rho>0 \f$ values \f$ \cos \theta \f$ and \f$ \sin \theta \f$ are computed.
!! If \f$ \sin \theta \ne 0 \f$ then \f$ \cos \phi \f$ and \f$ \sin \phi \f$
!! are computed.
!! Auxiliary values of associated Legendre polynomials \f$ P_\ell^m(\theta) \f$
!! are computed along with \f$ \cos (m \phi) \f$ and \f$ \sin(m \phi) \f$.
!!
!! @param[in] x: Target point
!! @param[out] rho: Euclidian length of `x`
!! @param[out] ctheta: \f$ \cos \theta \f$
!! @param[out] stheta: \f$ \sin \theta \f$
!! @param[out] cphi: \f$ \cos \phi \f$
!! @param[out] sphi: \f$ \sin \phi \f$
!! @param[in] p: Maximal degree of spherical harmonics. `p` >= 0
!! @param[in] vscales: Scaling factors of real normalized spherical harmonics.
!!      Dimension is `(p+1)**2`
!! @param[out] vylm: Values of spherical harmonics \f$ Y_\ell^m(x) \f$.
!!      Dimension is `(p+1)**2`
!! @param[out] vplm: Values of associated Legendre polynomials \f$ P_\ell^m(
!!      \theta) \f$. Dimension is `(p+1)**2`
!! @param[out] vcos: Array of alues of \f$ \cos(m\phi) \f$ of a dimension
!!      `(p+1)`
!! @param[out] vsin: array of values of \f$ \sin(m\phi) \f$ of a dimension
!!      `(p+1)`
subroutine ylmbas(x, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
        & vcos, vsin)
    ! Inputs
    real(kind=rp), intent(in) :: x(3)
    integer, intent(in) :: p
    real(kind=rp), intent(in) :: vscales((p+1)**2)
    ! Outputs
    real(kind=rp), intent(out) :: rho, ctheta, stheta, cphi, sphi
    real(kind=rp), intent(out) :: vylm((p+1)**2), vplm((p+1)**2)
    real(kind=rp), intent(out) :: vcos(p+1), vsin(p+1)
    ! Local variables
    integer :: l, m, ind
    real(kind=rp) :: max12, ssq12, tmp
    ! Get rho cos(theta), sin(theta), cos(phi) and sin(phi) from the cartesian
    ! coordinates of x. To support full range of inputs we do it via a scale
    ! and a sum of squares technique.
    ! At first we compute x(1)**2 + x(2)**2
    if (x(1) .eq. zero) then
        max12 = abs(x(2))
        ssq12 = one
    else if (abs(x(2)) .gt. abs(x(1))) then
        max12 = abs(x(2))
        ssq12 = one + (x(1)/x(2))**2
    else
        max12 = abs(x(1))
        ssq12 = one + (x(2)/x(1))**2
    end if
    ! Then we compute rho
    if (x(3) .eq. zero) then
        rho = max12 * sqrt(ssq12)
    else if (abs(x(3)) .gt. max12) then
        rho = one + ssq12 *(max12/x(3))**2
        rho = abs(x(3)) * sqrt(rho)
    else
        rho = ssq12 + (x(3)/max12)**2
        rho = max12 * sqrt(rho)
    end if
    if (rho .eq. zero) then
        return
    end if
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Evalutate cos(m*phi) and sin(m*phi) arrays. Notice that this is
    ! pointless if x(3) = 1, as the only non vanishing terms will be the
    ! ones with m = 0.
    if (stheta .ne. zero) then
        cphi = x(1) / stheta
        sphi = x(2) / stheta
        call trgev(cphi, sphi, p, vcos, vsin)
    else
        cphi = 1
        sphi = 0
        vcos = zero
        vsin = zero
    end if
    ! Normalize ctheta and stheta
    ctheta = x(3) / rho
    stheta = stheta / rho
    ! Evaluate associated Legendre polynomials
    call polleg(ctheta, stheta, p, vplm)
    ! Construct spherical harmonics
    do l = 0, p
        ! Offset of a Y_l^0 harmonic in vplm and vylm arrays
        ind = l**2 + l + 1
        ! m = 0
        vylm(ind) = vscales(ind) * vplm(ind)
        do m = 1, l
            ! only P_l^m for non-negative m is used/defined
            tmp = vplm(ind+m) * vscales(ind+m)
            ! m > 0
            vylm(ind+m) = tmp * vcos(m+1)
            ! m < 0
            vylm(ind-m) = tmp * vsin(m+1)
        end do
    end do
end subroutine ylmbas

!> Switching function
!!
!! This is an implementation of \f$ \chi(t) \f$ with a shift \f$ s \f$:
!! \f[
!!      \chi(t) = \left\{ \begin{array}{ll} 0, & \text{if} \quad x \geq
!!      1 \\ p_\eta(x), & \text{if} \quad 1-\eta < x < 1 \\ 1, & \text{if}
!!      \quad x \leq 1-\eta \end{array} \right.
!! \f]
!! where \f$ x = t - \frac{1+s}{2} \eta \f$ is a shifted coordinate and
!! \f[
!!      p_\eta(x) = \frac{1}{\eta^5} (1-t)^3 (6t^2 + (15\eta-12)t + (10\eta^2
!!      -15\eta+6))
!! \f]
!! is a smoothing polynomial of the 5th degree.
!! In the case shift \f$ s=1 \f$ the switch function \f$ \chi(t) \f$ is 1 for
!! \f$ t \in [0,1] \f$, is 0 for \f$ t \in [1+\eta, \infty) \f$ and varies in
!! \f$ [1, 1+\eta] \f$ which is an external shift.
!! In the case shift \f$ s=-1 \f$ the switch function \f$ \chi(t) \f$ is 1 for
!! \f$ t \in [0,1-\eta] \f$, is 0 for \f$ t \in [1, \infty) \f$ and varies in
!! \f$ [1-eta, 1] \f$ which is an internal shift.
!! In the case shift \f$ s=0 \f$ the switch function \f$ \chi(t) \f$ is 1 for
!! \f$ t \in [0,1-\eta/2] \f$, is 0 for \f$ t \in [1+\eta/2, \infty) \f$ and
!! varies in \f$ [1-eta/2, 1+\eta/2] \f$ which is a centered shift.
!!
!! @param[in] t: real non-negative input value
!! @param[in] s: shift
!! @param[in] eta: regularization parameter \f$ \eta \f$
real(kind=rp) function fsw(t, s, eta)
    ! Inputs
    real(kind=rp), intent(in) :: t, s, eta
    ! Local variables
    real(kind=rp) :: a, b, flow, x
    real(kind=rp), parameter :: f6=6.0d0, f10=10.d0, f12=12.d0, f15=15.d0
    ! Apply shift:
    !   s =  0   =>   t - eta/2  [ CENTERED ]
    !   s =  1   =>   t - eta    [ EXTERIOR ]
    !   s = -1   =>   t          [ INTERIOR ]
    x = t - (s + one)*eta / two
    ! Lower bound of switch region
    flow = one - eta
    ! Define switch function chi
    if (x .ge. one) then
        fsw = zero
    else if (x .le. flow) then
        fsw = one
    else
        a = f15*eta - f12
        b = f10*eta*eta - f15*eta + f6
        fsw = (one-x) ** 3
        fsw = fsw * (f6*x*x+a*x+b) / (eta**5)
    end if
end function fsw

!> Derivative of a switching function
!!
!! This is an implementation of \f$ \chi'(t) \f$ with a shift \f$ s \f$:
!! \f[
!!      \chi'(t) = \left\{ \begin{array}{ll} 0, & \text{if} \quad x \geq
!!      1 \\ p'_\eta(x), & \text{if} \quad 1-\eta < x < 1 \\ 0, & \text{if}
!!      \quad x \leq 1-\eta \end{array} \right.
!! \f]
!! where \f$ x = t - \frac{1+s}{2} \eta \f$ is a shifted coordinate and
!! \f[
!!      p'_\eta(x) = -\frac{30}{\eta^5} (1-t)^2 (t-1+\eta)^2
!! \f]
!! is a derivative of the smoothing polynomial.
!! In the case shift \f$ s=1 \f$ the derivative \f$ \chi'(t) \f$ is 0 for
!! \f$ t \in [0,1] \cup [1+\eta, \infty) \f$ and varies in
!! \f$ [1, 1+\eta] \f$ which is an external shift.
!! In the case shift \f$ s=-1 \f$ the derivative \f$ \chi'(t) \f$ is 0 for
!! \f$ t \in [0,1-\eta] \cup [1, \infty) \f$ and varies in
!! \f$ [1-eta, 1] \f$ which is an internal shift.
!! In the case shift \f$ s=0 \f$ the derivative \f$ \chi'(t) \f$ is 0 for
!! \f$ t \in [0,1-\eta/2] \cup [1+\eta/2, \infty) \f$ and
!! varies in \f$ [1-eta/2, 1+\eta/2] \f$ which is a centered shift.
!!
!! @param[in] t: real non-negative input value
!! @param[in] s: shift
!! @param[in] eta: regularization parameter \f$ \eta \f$
real(kind=rp) function dfsw(t, s, eta)
    ! Inputs
    real(kind=rp), intent(in) :: t, s, eta
    ! Local variables
    real(kind=rp) :: flow, x
    real(kind=rp), parameter :: f30=30.0d0
    ! Apply shift:
    !   s =  0   =>   t - eta/2  [ CENTERED ]
    !   s =  1   =>   t - eta    [ EXTERIOR ]
    !   s = -1   =>   t          [ INTERIOR ]
    x = t - (s + 1.d0)*eta / 2.d0
    ! Lower bound of switch region
    flow = one - eta
    ! Define derivative of switch function chi
    if (x .ge. one) then
        dfsw = zero
    else if (x .le. flow) then
        dfsw = zero
    else
        dfsw = -f30 * (( (one-x)*(x-one+eta) )**2) / (eta**5)
    endif
end function dfsw

!> Integrate against spherical harmonics
!!
!! Integrate by Lebedev spherical quadrature. This function can be simply
!! substituted by a matrix-vector product.
!!
!! TODO: use dgemv. Loop of this cycle can be easily substituted by a dgemm.
!!
!! @param[in] ngrid: Number of Lebedev grid points. `ngrid` > 0
!! @param[in] p: Maximal degree of spherical harmonics. `p` >= 0
!! @param[in] vwgrid: Values of spherical harmonics at Lebedev grid points,
!!      multiplied by weights of grid points. Dimension is ((p+1)**2, ngrid)
!! @param[in] isph: Index of given sphere. Used for debug purpose.
!! @param[in] x: Input values at grid points of the sphere. Dimension is
!!      (ngrid)
!! @param[out] xlm: Output spherical harmonics. Dimension is ((p+1)**2)
subroutine intrhs(iprint, ngrid, p, vwgrid, isph, x, xlm)
    ! Inputs
    integer, intent(in) :: iprint, ngrid, p, isph
    real(kind=rp), intent(in) :: vwgrid((p+1)**2, ngrid)
    real(kind=rp), intent(in) :: x(ngrid)
    ! Output
    real(kind=rp), intent(out) :: xlm((p+1)**2)
    ! Local variables
    integer :: igrid
    ! Initialize
    xlm = zero
    ! Accumulate over integration points
    do igrid = 1, ngrid
        xlm = xlm + vwgrid(:,igrid)*x(igrid)
    end do
    ! Printing (these functions are not implemented yet)
    if (iprint .ge. 5) then
        call ptcart('pot', ngrid, 1, isph, x)
        call prtsph('vlm', (p+1)**2, p, 1, isph, xlm)
    end if
end subroutine intrhs

!> Compute first derivatives of spherical harmonics
!!
!! @param[in] x:
!! @param[out] basloc:
!! @param[out] dbsloc:
!! @param[out] vplm:
!! @param[out] vcos:
!! @param[out] vsin:
!!
!!
!! TODO: rewrite code and fill description. Computing sqrt(one-cthe*cthe)
!! reduces effective range of input double precision values. cthe*cthe for
!! cthe=1d+155 is NaN.
subroutine dbasis(dd_data, x, basloc, dbsloc, vplm, vcos, vsin)
    type(dd_data_type) :: dd_data
    real(kind=rp), dimension(3),        intent(in)    :: x
    real(kind=rp), dimension((dd_data % lmax+1)**2),   intent(inout) :: basloc, vplm
    real(kind=rp), dimension(3,(dd_data % lmax +1)**2), intent(inout) :: dbsloc
    real(kind=rp), dimension(dd_data % lmax+1),   intent(inout) :: vcos, vsin
    integer :: l, m, ind, VC, VS
    real(kind=rp)  :: cthe, sthe, cphi, sphi, plm, fln, pp1, pm1, pp
    real(kind=rp)  :: et(3), ep(3)
    !     get cos(\theta), sin(\theta), cos(\phi) and sin(\phi) from the cartesian
    !     coordinates of x.
    cthe = x(3)
    sthe = sqrt(one - cthe*cthe)
    !     not ( NORTH or SOUTH pole )
    if ( sthe.ne.zero ) then
        cphi = x(1)/sthe
        sphi = x(2)/sthe
        !     NORTH or SOUTH pole
    else
        cphi = one
        sphi = zero
    end if
    !     evaluate the derivatives of theta and phi:
    et(1) = cthe*cphi
    et(2) = cthe*sphi
    et(3) = -sthe
    !     not ( NORTH or SOUTH pole )
    if( sthe.ne.zero ) then
        ep(1) = -sphi/sthe
        ep(2) = cphi/sthe
        ep(3) = zero
        !     NORTH or SOUTH pole
    else
        ep(1) = zero
        ep(2) = one
        ep(3) = zero
    end if
    !     evaluate cos(m*phi) and sin(m*phi) arrays. notice that this is 
    !     pointless if z = 1, as the only non vanishing terms will be the 
    !     ones with m=0.
    !
    !     not ( NORTH or SOUTH pole )
    if ( sthe.ne.zero ) then
        call trgev( cphi, sphi, dd_data % lmax, vcos, vsin )
        !     NORTH or SOUTH pole
    else
        vcos = one
        vsin = zero
    end if
    VC=zero
    VS=cthe
    !     evaluate the generalized legendre polynomials.
    call polleg( cthe, sthe, dd_data % lmax, vplm )
    !
    !     now build the spherical harmonics. we will distinguish m=0,
    !     m>0 and m<0:
    !
    basloc = zero
    dbsloc = zero
    do l = 0, dd_data % lmax
        ind = l*l + l + 1
        ! m = 0
        fln = dd_data % vscales(ind)   
        basloc(ind) = fln*vplm(ind)
        if (l.gt.0) then
            dbsloc(:,ind) = fln*vplm(ind+1)*et(:)
        else
            dbsloc(:,ind) = zero
        end if
        !dir$ simd
        do m = 1, l
            fln = dd_data % vscales(ind+m)
            plm = fln*vplm(ind+m)   
            pp1 = zero
            if (m.lt.l) pp1 = -pt5*vplm(ind+m+1)
            pm1 = pt5*(dble(l+m)*dble(l-m+1)*vplm(ind+m-1))
            pp  = pp1 + pm1   
            !
            !         m > 0
            !         -----
            !
            basloc(ind+m) = plm*vcos(m+1)
            !          
            !         not ( NORTH or SOUTH pole )
            if ( sthe.ne.zero ) then
                ! 
                dbsloc(:,ind+m) = -fln*pp*vcos(m+1)*et(:) - dble(m)*plm*vsin(m+1)*ep(:)
                !
                !            
                !         NORTH or SOUTH pole
            else
                !                  
                dbsloc(:,ind+m) = -fln*pp*vcos(m+1)*et(:) - fln*pp*ep(:)*VC
                !
                !
            endif
            !
            !         m < 0
            !         -----
            !
            basloc(ind-m) = plm*vsin(m+1)
            !          
            !         not ( NORTH or SOUTH pole )
            if ( sthe.ne.zero ) then
                ! 
                dbsloc(:,ind-m) = -fln*pp*vsin(m+1)*et(:) + dble(m)*plm*vcos(m+1)*ep(:)
                !
                !         NORTH or SOUTH pole
            else
                !                  
                dbsloc(:,ind-m) = -fln*pp*vsin(m+1)*et(:) - fln*pp*ep(:)*VS
                !            
            endif
            !  
        enddo
    enddo
end subroutine dbasis

! Purpose : compute
!
!                               l
!             sum   4pi/(2l+1) t  * Y_l^m( s ) * sigma_l^m
!             l,m                                           
!
!           which is need to compute action of COSMO matrix L.
!------------------------------------------------------------------------------------------------
!
!!> TODO
real(kind=rp) function intmlp( dd_data, t, sigma, basloc )
!  
      implicit none
      type(dd_data_type) :: dd_data
      real(kind=rp), intent(in) :: t
      real(kind=rp), dimension((dd_data % lmax+1)**2), intent(in) :: sigma, basloc
!
      integer :: l, ind
      real(kind=rp)  :: tt, ss, fac
!
!------------------------------------------------------------------------------------------------
!
!     initialize t^l
      tt = one
!
!     initialize
      ss = zero
!
!     loop over l
      do l = 0, dd_data % lmax
!      
        ind = l*l + l + 1
!
!       update factor 4pi / (2l+1) * t^l
        fac = tt / dd_data % vscales(ind)**2
!
!       contract over l,m and accumulate
        ss = ss + fac * dot_product( basloc(ind-l:ind+l), &
                                     sigma( ind-l:ind+l)   )
!
!       update t^l
        tt = tt*t
!        
      enddo
!      
!     redirect
      intmlp = ss
!
!
end function intmlp

! Purpose : weigh potential at cavity points by characteristic function "ui"
!------------------------------------------------------------------------------------------------
!> TODO
subroutine wghpot( dd_data, phi, g )
!
      implicit none
!
    type(dd_data_type) :: dd_data
    real(kind=rp), dimension(dd_data % ncav), intent(in)  :: phi
    real(kind=rp), dimension(dd_data % ngrid, dd_data % nsph), intent(out) :: g
!
    integer isph, ig, ic
!
!------------------------------------------------------------------------------------------------
!
!   initialize
    ic = 0 ; g(:,:)=0.d0
!      
!   loop over spheres
    do isph = 1, dd_data % nsph
!
!   loop over points
      do ig = 1, dd_data % ngrid
!
!       nonzero contribution from point
        if ( dd_data % ui(ig,isph).ne.zero ) then
!
!         advance cavity point counter
          ic = ic + 1
!            
!         weigh by (negative) characteristic function
          g(ig,isph) = -dd_data % ui(ig,isph) * phi(ic)
        endif
!          
       enddo
   enddo
end subroutine wghpot

! Purpose : compute H-norm
!------------------------------------------------------------------------------------------------
!> TODO
subroutine hsnorm( dd_data, u, unorm )
!          
      implicit none
      type(dd_data_type) :: dd_data
      real(kind=rp), dimension((dd_data % lmax+1)**2), intent(in)    :: u
      real(kind=rp),                    intent(inout) :: unorm
!
      integer :: l, m, ind
      real(kind=rp)  :: fac
!
!------------------------------------------------------------------------------------------------
!
!     initialize
      unorm = zero
!      
!     loop over l
      do l = 0, dd_data % lmax
!      
!       first index associated to l
        ind = l*l + l + 1
!
!       scaling factor
        fac = one/(one + dble(l))
!
!       loop over m
        do m = -l, l
!
!         accumulate
          unorm = unorm + fac*u(ind+m)*u(ind+m)
!          
        enddo
      enddo
!
!     the much neglected square root
      unorm = sqrt(unorm)
!
      return
!
!
end subroutine hsnorm

! compute the h^-1/2 norm of the increment on each sphere, then take the
! rms value.
!-------------------------------------------------------------------------------
!
!> TODO
real(kind=rp) function hnorm(dd_data, x)
    type(dd_data_type), intent(in) :: dd_data
      real(kind=rp),  dimension(dd_data % nbasis, dd_data % nsph), intent(in) :: x
!
      integer                                     :: isph, istatus
      real(kind=rp)                                      :: vrms, vmax
      real(kind=rp), allocatable                         :: u(:)
!
!-------------------------------------------------------------------------------
!
!     allocate workspace
      allocate( u(dd_data % nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'hnorm: allocation failed !'
        stop
      endif
!
!     loop over spheres
      do isph = 1, dd_data % nsph
!
!       compute norm contribution
        call hsnorm(dd_data, x(:,isph), u(isph))
      enddo
!
!     compute rms of norms
      call rmsvec( dd_data % nsph, u, vrms, vmax )
!
!     return value
      hnorm = vrms
!
!     deallocate workspace
      deallocate( u , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'hnorm: deallocation failed !'
        stop
      endif
!
!
end function hnorm

!------------------------------------------------------------------------------
! Purpose : compute root-mean-square and max norm
!------------------------------------------------------------------------------
subroutine rmsvec( n, v, vrms, vmax )
!
      implicit none
      integer,               intent(in)    :: n
      real(kind=rp),  dimension(n), intent(in)    :: v
      real(kind=rp),                intent(inout) :: vrms, vmax
!
      integer :: i
      real(kind=rp), parameter :: zero=0.0d0
!      
!------------------------------------------------------------------------------
!      
!     initialize
      vrms = zero
      vmax = zero
!
!     loop over entries
      do i = 1,n
!
!       max norm
        vmax = max(vmax,abs(v(i)))
!
!       rms norm
        vrms = vrms + v(i)*v(i)
!        
      enddo
!
!     the much neglected square root
      vrms = sqrt(vrms/dble(n))
!      
      return
!      
!      
endsubroutine rmsvec

!-----------------------------------------------------------------------------------
! Purpose : compute
!
!   v_l^m = v_l^m +
!
!               4 pi           l
!     sum  sum  ---- ( t_n^ji )  Y_l^m( s_n^ji ) W_n^ji [ \xi_j ]_n
!      j    n   2l+1
!
! which is related to the action of the adjont COSMO matrix L^* in the following
! way. Set
!
!   [ \xi_j ]_n = sum  Y_l^m( s_n ) [ s_j ]_l^m
!                 l,m
!
! then
!
!   v_l^m = -   sum    (L^*)_ij s_j
!             j \ne i 
!
! The auxiliary quantity [ \xi_j ]_l^m needs to be computed explicitly.
!-----------------------------------------------------------------------------------
!
!> TODO
subroutine adjrhs( dd_data, isph, xi, vlm, basloc, vplm, vcos, vsin )
!
      implicit none
      type(dd_data_type) :: dd_data
      integer,                       intent(in)    :: isph
      real(kind=rp), dimension(dd_data % ngrid, dd_data % nsph), intent(in)    :: xi
      real(kind=rp), dimension((dd_data % lmax+1)**2),     intent(inout) :: vlm
      real(kind=rp), dimension((dd_data % lmax+1)**2),     intent(inout) :: basloc, vplm
      real(kind=rp), dimension(dd_data % lmax+1),     intent(inout) :: vcos, vsin
!
      integer :: ij, jsph, ig, l, ind, m
      real(kind=rp)  :: vji(3), vvji, tji, sji(3), xji, oji, fac, ffac, t
!      
!-----------------------------------------------------------------------------------
!
!     loop over neighbors of i-sphere
      do ij = dd_data % inl(isph), dd_data % inl(isph+1)-1
!
!       j-sphere is neighbor
        jsph = dd_data % nl(ij)
!
!       loop over integration points
        do ig = 1, dd_data % ngrid
!        
!         compute t_n^ji = | r_j + \rho_j s_n - r_i | / \rho_i
          vji  = dd_data % csph(:,jsph) + dd_data % rsph(jsph)* &
              & dd_data % cgrid(:,ig) - dd_data % csph(:,isph)
          vvji = sqrt(dot_product(vji,vji))
          tji  = vvji/dd_data % rsph(isph)
!
!         point is INSIDE i-sphere (+ transition layer)
!         ---------------------------------------------
          if ( tji.lt.( one + (se+one)/two*dd_data % eta ) ) then
!                  
!           compute s_n^ji
            sji = vji/vvji
!
!           compute \chi( t_n^ji )
            xji = fsw( tji, se, dd_data % eta )
!
!           compute W_n^ji
            if ( dd_data % fi(ig,jsph).gt.one ) then
!                    
              oji = xji/ dd_data % fi(ig,jsph)
!              
            else
!                    
              oji = xji
!              
            endif
!            
!           compute Y_l^m( s_n^ji )
            !call ylmbas( sji, basloc, vplm, vcos, vsin )
!            
!           initialize ( t_n^ji )^l
            t = one
!            
!           compute w_n * xi(n,j) * W_n^ji
            fac = dd_data % wgrid(ig) * xi(ig,jsph) * oji
!            
!           loop over l
            do l = 0, dd_data % lmax
!            
              ind  = l*l + l + 1
!
!             compute 4pi / (2l+1) * ( t_n^ji )^l * w_n * xi(n,j) * W_n^ji
              ffac = fac*t/ dd_data % vscales(ind)**2
!
!             loop over m
              do m = -l,l
!              
                vlm(ind+m) = vlm(ind+m) + ffac*basloc(ind+m)
!                
              enddo
!
!             update ( t_n^ji )^l
              t = t*tji
!              
            enddo
!            
          endif
        enddo
      enddo
!
!
end subroutine adjrhs

!------------------------------------------------------------------------
! Purpose : compute
!
!   \Phi( n ) =
!     
!                       4 pi           l
!     sum  W_n^ij  sum  ---- ( t_n^ij )  Y_l^m( s_n^ij ) [ \sigma_j ]_l^m
!      j           l,m  2l+1
!
! which is related to the action of the COSMO matrix L in the following
! way :
!
!   -   sum    L_ij \sigma_j = sum  w_n Y_l^m( s_n ) \Phi( n ) 
!     j \ne i                   n
!
! This second step is performed by routine "intrhs".
!------------------------------------------------------------------------
!
!> TODO
subroutine calcv( dd_data, first, isph, pot, sigma, basloc, vplm, vcos, vsin )
!
    type(dd_data_type) :: dd_data
      logical,                        intent(in)    :: first
      integer,                        intent(in)    :: isph
      real(kind=rp), dimension((dd_data % lmax+1)**2, dd_data % nsph), intent(in)    :: sigma
      real(kind=rp), dimension(dd_data % ngrid),       intent(inout) :: pot
      real(kind=rp), dimension((dd_data % lmax+1)**2),      intent(inout) :: basloc
      real(kind=rp), dimension((dd_data % lmax+1)**2),      intent(inout) :: vplm
      real(kind=rp), dimension(dd_data % lmax+1),      intent(inout) :: vcos
      real(kind=rp), dimension(dd_data % lmax+1),      intent(inout) :: vsin
!
      integer :: its, ij, jsph
      real(kind=rp)  :: vij(3), sij(3)
      real(kind=rp)  :: vvij, tij, xij, oij, stslm, stslm2, stslm3, &
          & thigh, rho, ctheta, stheta, cphi, sphi
!
!------------------------------------------------------------------------
      thigh = one + pt5*(se + one)*dd_data % eta
!
!     initialize
      pot(:) = zero
!
!     if 1st iteration of Jacobi method, then done!
      if ( first )  return
!
!     loop over grid points
      do its = 1, dd_data % ngrid
!
!       contribution from integration point present
        if ( dd_data % ui(its,isph).lt.one ) then
!
!         loop over neighbors of i-sphere
          do ij = dd_data % inl(isph), dd_data % inl(isph+1)-1
!
!           neighbor is j-sphere
            jsph = dd_data % nl(ij)
!            
!           compute t_n^ij = | r_i + \rho_i s_n - r_j | / \rho_j
            vij  = dd_data % csph(:,isph) + dd_data % rsph(isph)* &
                & dd_data % cgrid(:,its) - dd_data % csph(:,jsph)
            vvij = sqrt( dot_product( vij, vij ) )
            tij  = vvij / dd_data % rsph(jsph) 
!
!           point is INSIDE j-sphere
!           ------------------------
            if ( tij.lt.thigh .and. tij.gt.zero) then
!
!             compute s_n^ij = ( r_i + \rho_i s_n - r_j ) / | ... |
              sij = vij / vvij
!            
!             compute \chi( t_n^ij )
              xij = fsw( tij, se, dd_data % eta )
!
!             compute W_n^ij
              if ( dd_data % fi(its,isph).gt.one ) then
!
                oij = xij / dd_data % fi(its,isph)
!
              else
!
                oij = xij
!
              endif
!
!             compute Y_l^m( s_n^ij )
              !call ylmbas( sij, basloc, vplm, vcos, vsin )
              call ylmbas(sij, rho, ctheta, stheta, cphi, sphi, dd_data % lmax, &
                  & dd_data % vscales, basloc, vplm, vcos, vsin)
!                    
!             accumulate over j, l, m
              pot(its) = pot(its) + oij * intmlp( dd_data, tij, sigma(:,jsph), basloc )
!              
            endif
          end do
        end if
      end do
!      
!      
!      
end subroutine calcv
!------------------------------------------------------------------------
!
!------------------------------------------------------------------------
! Purpose : compute
!
! \xi(n,i) = 
!
!  sum w_n U_n^i Y_l^m(s_n) [S_i]_l^m
!  l,m
!
!------------------------------------------------------------------------
subroutine ddmkxi( dd_data, s, xi)
!
    type(dd_data_type) :: dd_data
       real(kind=rp), dimension((dd_data % lmax+1)**2, dd_data % nsph), intent(in)    :: s
       real(kind=rp), dimension(dd_data % ncav),      intent(inout) :: xi
!
       integer :: its, isph, ii
!
       ii = 0
       do isph = 1, dd_data % nsph
         do its = 1, dd_data % ngrid
           if (dd_data % ui(its,isph) .gt. zero) then
             ii     = ii + 1
             xi(ii) = dd_data % ui(its,isph)*dot_product(dd_data % vwgrid(:,its),s(:,isph))
           end if
         end do
       end do
!
       return
end subroutine ddmkxi
!
!------------------------------------------------------------------------
! Purpose : compute
!
! \zeta(n,i) = 
!
!  1/2 f(\eps) sum w_n U_n^i Y_l^m(s_n) [S_i]_l^m
!              l,m
!
!------------------------------------------------------------------------
subroutine ddmkzeta( dd_data, s, zeta)
!
    type(dd_data_type) :: dd_data
       real(kind=rp), dimension((dd_data % lmax+1)**2, dd_data % nsph), intent(in)    :: s
       real(kind=rp), dimension(dd_data % ncav),      intent(inout) :: zeta
!
       integer :: its, isph, ii
!
       ii = 0
       do isph = 1, dd_data % nsph
         do its = 1, dd_data % ngrid
           if (dd_data % ui(its,isph) .gt. zero) then
             ii     = ii + 1
             zeta(ii) = dd_data % ui(its,isph)* &
                 & dot_product(dd_data % vwgrid(:,its),s(:,isph))
           end if
         end do
       end do
!
       zeta = pt5*((dd_data % eps-one)/dd_data % eps)*zeta
       return
end subroutine ddmkzeta

!> Compute preconditioner
!!
!! assemble the diagonal blocks of the reps matrix
!! then invert them to build the preconditioner
subroutine mkprec(dd_data)
    ! Inouts
    type(dd_data_type), intent(inout) :: dd_data
    integer :: isph, lm, ind, l1, m1, ind1, its, istatus
    real(kind=rp)  :: f, f1
    integer, allocatable :: ipiv(:)
    real(kind=rp),  allocatable :: work(:)
    external :: dgetrf, dgetri
    ! Allocation of temporaries
    allocate(ipiv(dd_data % nbasis),work(dd_data % nbasis**2),stat=istatus)
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

!> Compute multipole coefficients for a particle of a unit charge
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{\ell+1} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^{-l-1} \f$.
!! This function is not needed by ddX, but it is used for testing purposes.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] r: Radius of spherical harmonics
!! @param[in] p: Maximal degral of multipole basis functions
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$. Dimension
!!      is `(p+1)**2`
!! @param[out] m: Multipole coefficients. Dimension is `(p+1)**2`
subroutine fmm_p2m(c, r, p, vscales, m)
    ! Inputs
    real(kind=rp), intent(in) :: c(3), r, vscales((p+1)**2)
    integer, intent(in) :: p
    ! Output
    real(kind=rp), intent(out) :: m((p+1)**2)
    ! Local variables
    real(kind=rp) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=rp) :: vylm((p+1)**2), vplm((p+1)**2), t, tmp, rcoef
    integer :: n, k, ind
    ! Get radius and values of spherical harmonics
    call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
        & vcos, vsin)
    ! Harmonics are available only if rho > 0
    if (rho .ne. 0) then
        rcoef = rho / r
        t = 1 / r
        do n = 0, p
            ind = n*n + n + 1
            m(ind-n:ind+n) = t * vylm(ind-n:ind+n)
            t = t * rcoef
        end do
    ! Naive case of rho = 0
    else
        m(1) = 1 / r
        m(2:) = 0
    end if
end subroutine fmm_p2m

!> Compute potential, induced by multipole spherical harmonics
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{\ell+1} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^{-l-1} \f$.
!! This function is not needed by ddX, but it is used for testing purposes.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of multipole basis functions
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$. Dimension
!!      is `(p+1)**2`
!! @param[in] m: Multipole coefficients. Dimension is `(p+1)**2`
!! @param[out] v: Value of induced potential
subroutine fmm_m2p(c, r, p, vscales, m, v)
    ! Inputs
    real(kind=8), intent(in) :: c(3), r, vscales((p+1)*(p+1)), m((p+1)*(p+1))
    integer, intent(in) :: p
    ! Output
    real(kind=8), intent(out) :: v
    ! Local variables
    real(kind=8) :: tmp, tmp1, tmp2
    real(kind=8) :: rho, t, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(kind=8) :: vylm((p+1)**2), vplm((p+1)**2), rcoef
    integer :: n, k, ind
    ! Get radius and values of spherical harmonics
    call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
        & vcos, vsin)
    ! Exit in case of singularity
    if (rho .eq. 0) then
        return
    end if
    ! Compute actual potential
    rcoef = r / rho
    t = 1
    v = 0
    do n = 0, p
        t = t * rcoef
        ind = n*n + n + 1
        tmp = dot_product(vylm(ind-n:ind+n), m(ind-n:ind+n))
        ! Here vscales(ind)**2 is 4*pi/sqrt(2n+1)
        v = v + t*tmp/vscales(ind)**2
    end do
end subroutine fmm_m2p

end module dd_core

