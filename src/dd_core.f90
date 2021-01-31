!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/dd_core.f90
!! Core routines and parameters of entire ddX software
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2021-01-31

!> Core routines and parameters of ddX software
module dd_core
implicit none

!> Kind of double precision
integer, parameter :: dp = kind(1d0)

!! Compile-time constants
real(dp), parameter :: zero = 0d0, one = 1d0, two = 2d0, three = 3d0
real(dp), parameter :: four = 4d0, pt5 = 5d-1
real(dp), parameter :: sqrt2 = sqrt(two)
real(dp), parameter :: pi4 = atan(one)
real(dp), parameter :: pi = four * pi4
real(dp), parameter :: fourpi = four * pi
real(dp), parameter :: twopi = two * pi
real(dp), parameter :: sqrt4pi = four * sqrt(pi4)
real(dp), parameter :: machine_eps = epsilon(zero)
!> Number of supported Lebedev grids
integer, parameter :: nllg = 32
!> Number of grid points of each Lebedev grid
integer, parameter :: ng0(nllg) = (/ 6, 14, 26, 38, 50, 74, 86, 110, 146, &
    & 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, &
    & 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810 /)

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
    real(dp), allocatable :: csph(:, :)
    !> Array of radii of atoms of dimension (nsph)
    real(dp), allocatable :: rsph(:)
    !> Dielectric permittivity inside
    real(dp) :: epsin
    !> Dielectric permittivity outside
    real(dp) :: epsout
    !> Relative dielectric permittivity
    real(dp) :: eps
    !> Debye H\"{u}ckel parameter
    real(dp) :: kappa
    !> Shift of characteristic function
    !!
    !! -1 for interior, 0 for centered and 1 for outer regularization. ddCOSMO
    !! and ddLPB use value -1 while ddPCM uses value 0 by default.
    real(dp) :: se
    !> Regularization parameter
    real(dp) :: eta
    !> Relative error threshold for iterative solvers
    real(dp) :: tol
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
    !> Total number of used real spherical harmonics and a size of `vscales`
    integer :: nscales
    !> Scales of real normalized spherical harmonics of degree up to dmax
    !!
    !! This array has `nscales` entries
    real(dp), allocatable :: vscales(:)
    !> Maximum sqrt of factorial precomputed
    !!
    !! Just like with `dmax` parameter, number of used factorials is either
    !! `2*lmax+1` or `2*(pm+pl)+1` depending on whether FMM is used or not
    integer :: nfact
    !> Array of square roots of factorials of dimension (nfact)
    real(dp), allocatable :: vfact(:)
    !> Coordinates of Lebedev quadrature points of dimension (3, ngrid)
    real(dp), allocatable :: cgrid(:, :)
    !> Weights of Lebedev quadrature points of dimension (ngrid)
    real(dp), allocatable :: wgrid(:)
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
    real(dp), allocatable :: vgrid(:, :)
    !> Weighted values of spherical harmonics at Lebedev grid points
    !!
    !! vwgrid(:, igrid) = vgrid(:, igrid) * wgrid(igrid)
    !! Dimensions of this array are ((vgrid_nbasis, ngrid)
    real(dp), allocatable :: vwgrid(:, :)
    !> Maximum number of neighbours per sphere. Input from a user
    integer :: nngmax
    !> List of intersecting spheres in a CSR format
    integer, allocatable :: inl(:)
    !> List of intersecting spheres in a CSR format
    integer, allocatable :: nl(:)
    !> Characteristic function f. Dimension is (ngrid, npsh)
    real(dp), allocatable :: fi(:, :)
    !> Characteristic function U. Dimension is (ngrid, nsph)
    real(dp), allocatable :: ui(:, :)
    !> Derivative of characteristic function U. Dimension is (3, ngrid, nsph)
    real(dp), allocatable :: zi(:, :, :)
    !> Number of external Lebedev grid points on a molecule surface
    integer :: ncav
    !> Coordinates of external Lebedev grid points of dimension (3, ncav)
    real(dp), allocatable :: ccav(:, :)
    !> Preconditioner for operator R_eps. Dimension is (nbasis, nbasis, nsph)
    real(dp), allocatable :: rx_prc(:, :, :)
    !> Maximum degree of spherical harmonics for M (multipole) expansion
    integer :: pm
    !> Maximum degree of spherical harmonics for L (local) expansion
    integer :: pl
    !! Cluster tree information
    !> Reordering of spheres for better locality of dimension (nsph)
    integer, allocatable :: order(:)
    !> Number of clusters
    integer :: nclusters
    !> The first and the last spheres of each node. Dimension is (2, nclusters)
    integer, allocatable :: cluster(:, :)
    !> Children of each cluster. Dimension is (2, nclusters)
    integer, allocatable :: children(:, :)
    !> Parent of each cluster. Dimension is (nclusters)
    integer, allocatable :: parent(:)
    !> Center of bounding sphere of each cluster. Dimension is (3, nclusters)
    real(dp), allocatable :: cnode(:, :)
    !> Radius of bounding sphere of each cluster. Dimension is (nclusters)
    real(dp), allocatable :: rnode(:)
    !> Which leaf node contains only given input sphere. Dimension is (nsph)
    integer, allocatable :: snode(:)
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
    !! Reflection matrices for M2M, M2L and L2L operations
    !> Size of each M2M reflection matrix
    integer :: m2m_reflect_mat_size
    !> Array of M2M reflection matrices
    real(dp), allocatable :: m2m_reflect_mat(:, :)
    !> Size of each M2M OZ translation matrices
    integer :: m2m_ztranslate_mat_size
    !> Array of M2M OZ translation matrices
    real(dp), allocatable :: m2m_ztranslate_mat(:, :)
    !> Size of each L2L reflection matrix
    integer :: l2l_reflect_mat_size
    !> Array of L2L reflection matrices
    real(dp), allocatable :: l2l_reflect_mat(:, :)
    !> Size of each L2L OZ translation matrices
    integer :: l2l_ztranslate_mat_size
    !> Array of L2L OZ translation matrices
    real(dp), allocatable :: l2l_ztranslate_mat(:, :)
    !> Size of each M2L reflection matrix
    integer :: m2l_reflect_mat_size
    !> Array of M2L reflection matrices
    real(dp), allocatable :: m2l_reflect_mat(:, :)
    !> Size of each M2L OZ translation matrices
    integer :: m2l_ztranslate_mat_size
    !> Array of M2L OZ translation matrices
    real(dp), allocatable :: m2l_ztranslate_mat(:, :)
    !> Far-field L2P matrices
    real(dp), allocatable :: l2p_mat(:, :)
    !> Near-field M2P matrices
    real(dp), allocatable :: m2p_mat(:, :)

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
!!      actual number of grid points on exit. `ngrid` >= 0
!! @param[in] force: 1 if forces are required and 0 otherwise
!! @param[in] fmm: 1 to use FMM acceleration and 0 otherwise
!! @param[in] pm: Maximal degree of multipole spherical harmonics. `pm` >= 0
!! @param[in] pl: Maximal degree of local spherical harmonics. `pl` >= 0
!! @param[in] iprint: Level of printing debug info
!! @param[in] nngmax: Maximal number of neighbours (intersects) for every
!!      sphere
!! @param[in] se: Shift of characteristic function. -1 for interior, 0 for
!!      centered and 1 for outer regularization
!! @param[in] eta: Regularization parameter
!! @param[in] eps: Relative dielectric permittivity
!! @param[in] kappa: Debye H\"{u}ckel parameter
!! @param[out] dd_data: Object containing all inputs
!! @param[out] info: flag of succesfull exit
!!      = 0: Succesfull exit
!!      < 0: If info=-i then i-th argument had an illegal value
!!      > 0: Allocation of a buffer for the output dd_data failed
subroutine ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, pl, &
        & iprint, nngmax, se, eta, eps, kappa, dd_data, info)
    ! Inputs
    integer, intent(in) :: n, model, lmax, force, fmm, pm, pl, iprint, nngmax
    real(dp), intent(in):: x(n), y(n), z(n), rvdw(n), se, eta, eps, kappa
    ! Output
    type(dd_data_type), intent(out) :: dd_data
    integer, intent(out) :: info
    ! Inouts
    integer, intent(inout) :: ngrid
    ! Local variables
    integer :: istatus, i, ii, inear, jnear, igrid, jgrid, isph, jsph, lnl
    real(dp) :: rho, ctheta, stheta, cphi, sphi
    real(dp), allocatable :: vplm(:), vcos(:), vsin(:)
    real(dp) :: v(3), vv, t, swthr, fac, d2, r2
    double precision, external :: dnrm2
    ! Reset info
    info = 0
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
    ! No checks for NaNs and Infs in input coordinates
    dd_data % csph(1, :) = x
    dd_data % csph(2, :) = y
    dd_data % csph(3, :) = z
    dd_data % rsph = rvdw
    if ((model .lt. 1) .or. (model .gt. 3)) then
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
    if ((force .lt. 0) .or. (force .gt. 1)) then
        !write(*, *) "ddinit: wrong value of parameter `force`"
        info = -9
        return
    end if
    dd_data % force = force
    if ((fmm .lt. 0) .or. (fmm .gt. 1)) then
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
    if ((se .lt. -1) .or. (se .gt. 1)) then
        !write(*, *) "ddinit: wrong value of parameter `se`"
        info = -15
        return
    end if
    dd_data % se = se
    if ((eta .lt. 0) .or. (eta .gt. 1)) then
        !write(*, *) "ddinit: wrong value of parameter `eta`"
        info = -16
        return
    end if
    dd_data % eta = eta
    if (eps .lt. 0) then
        !write(*, *) "ddinit: wrong value of parameter `eps`"
        info = -17
        return
    end if
    dd_data % eps = eps
    if (kappa .lt. 0) then
        !write(*, *) "ddinit: wrong value of parameter `kappa`"
        info = -18
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
        dd_data % pm = pm
        dd_data % pl = pl
    else
        ! These values are ignored if fmm flag is 0
        dd_data % pm = -1
        dd_data % pl = -1
    end if
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
    dd_data % vgrid_nbasis = (dd_data % vgrid_dmax+1) ** 2
    dd_data % nfact = max(2*dd_data % dmax+1, 2)
    dd_data % nscales = (dd_data % dmax+1) ** 2
    ! Compute scaling factors of spherical harmonics
    allocate(dd_data % vscales(dd_data % nscales), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "ddinit: [2] allocation failed!"
        info = 2
        return
    end if
    call ylmscale(dd_data % dmax, dd_data % vscales)
    ! Precompute square roots of factorials
    allocate(dd_data % vfact(dd_data % nfact), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "ddinit: [3] allocation failed!"
        info = 3
        return
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
        if (jnear .lt. inear) then
            inear = jnear
            igrid = i
        end if
    end do
    ! Update inout parameter `ngrid` also
    ngrid = ng0(igrid)
    dd_data % ngrid = ngrid
    ! Get weights and coordinates of Lebedev grid points
    allocate(dd_data % cgrid(3, ngrid), dd_data % wgrid(ngrid), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "ddinit: [4] allocation failed!"
        info = 4
        return
    end if
    call llgrid(ngrid, dd_data % wgrid, dd_data % cgrid)
    ! Compute non-weighted and weighted spherical harmonics at grid points
    allocate(dd_data % vgrid(dd_data % vgrid_nbasis, ngrid), &
        & dd_data % vwgrid(dd_data % vgrid_nbasis, ngrid), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "ddinit: [5] allocation failed!"
        info = 5
        return
    end if
    allocate(vplm(dd_data % vgrid_nbasis), vcos(dd_data % vgrid_dmax+1), &
        & vsin(dd_data % vgrid_dmax+1), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) "ddinit: [6] allocation failed!"
        info = 6
        return
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
        !write(*, *) "ddinit: [7] deallocation failed!"
        info = 7
        return
    end if
    ! Upper bound of switch region. Defines intersection criterion for spheres
    swthr = one + (se+one)*eta/two
    ! Build list of neighbours in CSR format
    allocate(dd_data % inl(n+1), dd_data % nl(n*nngmax), stat=istatus)
    if (istatus .ne. 0) then
        !write(*, *) 'ddinit : [8] allocation failed !'
        info = 8
        return
    end if
    i  = 1
    lnl = 0
    do isph = 1, n
        dd_data % inl(isph) = lnl + 1
        do jsph = 1, dd_data % nsph
            if (isph .ne. jsph) then
                d2 = dnrm2(3, &
                    & dd_data % csph(:, isph)-dd_data % csph(:, jsph), 1)
                ! Take regularization parameter into account with respect to
                ! shift se. It is described properly by the upper bound of a
                ! switch region `swthr`.
                r2 = rvdw(isph) + swthr*rvdw(jsph)
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
        !write(*, *) 'ddinit : [9] allocation failed !'
        info = 9
        return
    end if
    dd_data % fi = zero
    dd_data % ui = zero
    if (force .eq. 1) then
        allocate(dd_data % zi(3, ngrid, n), stat=istatus)
        if (istatus .ne. 0) then
            !write(*, *) 'ddinit : [10] allocation failed !'
            info = 10
            return
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
        !write(*,*)'ddinit : [11] allocation failed!'
        info = 11
        return
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
            !write(*,*)'ddinit : [12] allocation failed!'
            info = 12
            return
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
    !! Prepare FMM structures if needed
    if (fmm .eq. 1) then
        ! Allocate space for a cluster tree
        allocate(dd_data % order(n), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [13] allocation failed!'
            info = 13
            return
        endif
        dd_data % nclusters = 2*n - 1
        allocate(dd_data % cluster(2, dd_data % nclusters), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [14] allocation failed!'
            info = 14
            return
        endif
        allocate(dd_data % children(2, dd_data % nclusters), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [15] allocation failed!'
            info = 15
            return
        endif
        allocate(dd_data % parent(dd_data % nclusters), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [16] allocation failed!'
            info = 16
            return
        endif
        allocate(dd_data % cnode(3, dd_data % nclusters), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [17] allocation failed!'
            info = 17
            return
        endif
        allocate(dd_data % rnode(dd_data % nclusters), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [18] allocation failed!'
            info = 18
            return
        endif
        allocate(dd_data % snode(dd_data % nsph), stat=istatus)
        if (istatus .ne. 0) then
            !write(*,*)'ddinit : [19] allocation failed!'
            info = 19
            return
        endif
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
    if (allocated(dd_data % csph)) then
        deallocate(dd_data % csph, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [1] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(dd_data % rsph)) then
        deallocate(dd_data % rsph, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [2] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(dd_data % vscales)) then
        deallocate(dd_data % vscales, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [3] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(dd_data % vfact)) then
        deallocate(dd_data % vfact, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [4] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(dd_data % cgrid)) then
        deallocate(dd_data % cgrid, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [5] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(dd_data % wgrid)) then
        deallocate(dd_data % wgrid, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [6] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(dd_data % vgrid)) then
        deallocate(dd_data % vgrid, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [7] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(dd_data % vwgrid)) then
        deallocate(dd_data % vwgrid, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [8] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(dd_data % inl)) then
        deallocate(dd_data % inl, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [9] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(dd_data % nl)) then
        deallocate(dd_data % nl, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [10] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(dd_data % fi)) then
        deallocate(dd_data % fi, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [11] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(dd_data % ui)) then
        deallocate(dd_data % ui, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [12] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(dd_data % zi)) then
        deallocate(dd_data % zi, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [13] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(dd_data % ccav)) then
        deallocate(dd_data % ccav, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [14] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(dd_data % rx_prc)) then
        deallocate(dd_data % rx_prc, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [15] deallocation failed!"
            stop 1
        end if
    end if
    if (allocated(dd_data % order)) then
        deallocate(dd_data % order, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [16] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(dd_data % cluster)) then
        deallocate(dd_data % cluster, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [17] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(dd_data % children)) then
        deallocate(dd_data % children, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [18] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(dd_data % parent)) then
        deallocate(dd_data % parent, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [19] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(dd_data % cnode)) then
        deallocate(dd_data % cnode, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [20] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(dd_data % rnode)) then
        deallocate(dd_data % rnode, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [21] deallocation failed!"
            stop 1
        endif
    end if
    if (allocated(dd_data % snode)) then
        deallocate(dd_data % snode, stat=istatus)
        if (istatus .ne. 0) then
            write(*, *) "ddfree: [22] deallocation failed!"
            stop 1
        endif
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
    real(dp), intent(in) :: x(nbasis, ncol)
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
    real(dp), intent(in) :: x(ngrid, ncol)
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
    real(dp) :: vector(dd_data % nbasis, dd_data % nsph)
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
!! @param[out] vscales: Array of scaling factors. Dimension is `nscales`
subroutine ylmscale(p, vscales)
    ! Input
    integer, intent(in) :: p
    ! Output
    real(dp), intent(out) :: vscales((p+1)**2)
    ! Local variables
    real(dp) :: tmp
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
    real(dp), intent(in) :: cphi, sphi
    integer, intent(in) :: p
    ! Output
    real(dp), intent(out) :: vcos(p+1), vsin(p+1)
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
    real(dp), intent(in) :: ctheta, stheta
    integer, intent(in) :: p
    ! Outputs
    real(dp), intent(out) :: vplm((p+1)**2)
    ! Local variables
    integer :: m, ind, l, ind2
    real(dp) :: fact, pmm, pmm1, pmmo, pll, fm, fl
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

!> Convert input cartesian coordinate into spherical coordinate
!!
!! Output coordinate \f$ (\rho, \theta, \phi) \f$ is presented by \f$ (\rho,
!! \cos \theta, \sin \theta, \cos \phi, \sin\phi) \f$.
!!
!! @param[in] x: Cartesian coordinate
!! @param[out] rho: \f$ \rho \f$
!! @param[out] ctheta: \f$ \cos \theta \f$
!! @param[out] stheta: \f$ \sin \theta \f$
!! @param[out] cphi: \f$ \cos \phi \f$
!! @param[out] sphi: \f$ \sin \phi \f$
subroutine carttosph(x, rho, ctheta, stheta, cphi, sphi)
    ! Input
    real(dp), intent(in) :: x(3)
    ! Output
    real(dp), intent(out) :: rho, ctheta, stheta, cphi, sphi
    ! Local variables
    real(dp) :: max12, ssq12
    ! Check x(1:2) = 0
    if ((x(1) .eq. zero) .and. (x(2) .eq. zero)) then
        rho = abs(x(3))
        ctheta = sign(one, x(3))
        stheta = zero
        cphi = one
        sphi = zero
        return
    end if
    ! In other situations use sum-of-scaled-squares technique
    ! Get norm of x(1:2) and cphi with sphi outputs
    if (abs(x(2)) .gt. abs(x(1))) then
        max12 = abs(x(2))
        ssq12 = one + (x(1)/x(2))**2
    else
        max12 = abs(x(1))
        ssq12 = one + (x(2)/x(1))**2
    end if
    stheta = max12 * sqrt(ssq12)
    cphi = x(1) / stheta
    sphi = x(2) / stheta
    ! Then compute rho, ctheta and stheta outputs
    if (abs(x(3)) .gt. max12) then
        rho = one + ssq12*(max12/x(3))**2
        rho = abs(x(3)) * sqrt(rho)
        stheta = stheta / rho
        ctheta = x(3) / rho
    else
        rho = ssq12 + (x(3)/max12)**2
        rho = max12 * sqrt(rho)
        stheta = stheta / rho
        ctheta = x(3) / rho
    end if
end subroutine carttosph

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
!! @param[out] ctheta: \f$ -1 \leq \cos \theta \leq 1\f$
!! @param[out] stheta: \f$ 0 \leq \sin \theta \leq 1\f$
!! @param[out] cphi: \f$ -1 \leq \cos \phi \leq 1\f$
!! @param[out] sphi: \f$ -1 \leq \sin \phi \leq 1\f$
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
    real(dp), intent(in) :: x(3)
    integer, intent(in) :: p
    real(dp), intent(in) :: vscales((p+1)**2)
    ! Outputs
    real(dp), intent(out) :: rho, ctheta, stheta, cphi, sphi
    real(dp), intent(out) :: vylm((p+1)**2), vplm((p+1)**2)
    real(dp), intent(out) :: vcos(p+1), vsin(p+1)
    ! Local variables
    integer :: l, m, ind
    real(dp) :: max12, ssq12, tmp
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
    ! In case x=0 just exit without setting any other variable
    if (rho .eq. zero) then
        return
    end if
    ! Length of a vector x(1:2)
    stheta = max12 * sqrt(ssq12)
    ! Case x(1:2) != 0
    if (stheta .ne. zero) then
        ! Evaluate cos(m*phi) and sin(m*phi) arrays
        cphi = x(1) / stheta
        sphi = x(2) / stheta
        call trgev(cphi, sphi, p, vcos, vsin)
        ! Normalize ctheta and stheta
        ctheta = x(3) / rho
        stheta = stheta / rho
        ! Evaluate associated Legendre polynomials
        call polleg(ctheta, stheta, p, vplm)
        ! Construct spherical harmonics
        do l = 0, p
            ! Offset of a Y_l^0 harmonic in vplm and vylm arrays
            ind = l**2 + l + 1
            ! m = 0 implicitly uses `vcos(1) = 1`
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
    ! Case of x(1:2) = 0 and x(3) != 0
    else
        ! Set spherical coordinates
        cphi = one
        sphi = zero
        ctheta = sign(one, x(3))
        stheta = zero
        ! Set output arrays vcos and vsin
        vcos = one
        vsin = zero
        ! Evaluate spherical harmonics. P_l^m = 0 for m > 0. In the case m = 0
        ! it depends if l is odd or even. Additionally, vcos = one and vsin =
        ! zero for all elements
        vylm = zero
        vplm = zero
        do l = 0, p, 2
            ind = l**2 + l + 1
            ! only case m = 0
            vplm(ind) = one
            vylm(ind) = vscales(ind)
        end do
        do l = 1, p, 2
            ind = l**2 + l + 1
            ! only case m = 0
            vplm(ind) = ctheta
            vylm(ind) = ctheta * vscales(ind)
        end do
    end if
end subroutine ylmbas

!> Switching function
!!
!! This is an implementation of \f$ \chi(t) \f$ with a shift \f$ se \f$:
!! \f[
!!      \chi(t) = \left\{ \begin{array}{ll} 0, & \text{if} \quad x \geq
!!      1 \\ p_\eta(x), & \text{if} \quad 1-\eta < x < 1 \\ 1, & \text{if}
!!      \quad x \leq 1-\eta \end{array} \right.
!! \f]
!! where \f$ x = t - \frac{1+se}{2} \eta \f$ is a shifted coordinate and
!! \f[
!!      p_\eta(x) = \frac{1}{\eta^5} (1-t)^3 (6t^2 + (15\eta-12)t + (10\eta^2
!!      -15\eta+6))
!! \f]
!! is a smoothing polynomial of the 5th degree.
!! In the case shift \f$ se=1 \f$ the switch function \f$ \chi(t) \f$ is 1 for
!! \f$ t \in [0,1] \f$, is 0 for \f$ t \in [1+\eta, \infty) \f$ and varies in
!! \f$ [1, 1+\eta] \f$ which is an external shift.
!! In the case shift \f$ se=-1 \f$ the switch function \f$ \chi(t) \f$ is 1 for
!! \f$ t \in [0,1-\eta] \f$, is 0 for \f$ t \in [1, \infty) \f$ and varies in
!! \f$ [1-\eta, 1] \f$ which is an internal shift.
!! In the case shift \f$ se=0 \f$ the switch function \f$ \chi(t) \f$ is 1 for
!! \f$ t \in [0,1-\eta/2] \f$, is 0 for \f$ t \in [1+\eta/2, \infty) \f$ and
!! varies in \f$ [1-\eta/2, 1+\eta/2] \f$ which is a centered shift.
!!
!! @param[in] t: real non-negative input value
!! @param[in] se: shift
!! @param[in] eta: regularization parameter \f$ \eta \f$
real(dp) function fsw(t, se, eta)
    ! Inputs
    real(dp), intent(in) :: t, se, eta
    ! Local variables
    real(dp) :: a, b, flow, x
    real(dp), parameter :: f6=6.0d0, f10=10.d0, f12=12.d0, f15=15.d0
    ! Apply shift:
    !   se =  0   =>   t - eta/2  [ CENTERED ]
    !   se =  1   =>   t - eta    [ EXTERIOR ]
    !   se = -1   =>   t          [ INTERIOR ]
    x = t - (se + one)*eta / two
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
!! This is an implementation of \f$ \chi'(t) \f$ with a shift \f$ se \f$:
!! \f[
!!      \chi'(t) = \left\{ \begin{array}{ll} 0, & \text{if} \quad x \geq
!!      1 \\ p'_\eta(x), & \text{if} \quad 1-\eta < x < 1 \\ 0, & \text{if}
!!      \quad x \leq 1-\eta \end{array} \right.
!! \f]
!! where \f$ x = t - \frac{1+se}{2} \eta \f$ is a shifted coordinate and
!! \f[
!!      p'_\eta(x) = -\frac{30}{\eta^5} (1-t)^2 (t-1+\eta)^2
!! \f]
!! is a derivative of the smoothing polynomial.
!! In the case shift \f$ se=1 \f$ the derivative \f$ \chi'(t) \f$ is 0 for
!! \f$ t \in [0,1] \cup [1+\eta, \infty) \f$ and varies in
!! \f$ [1, 1+\eta] \f$ which is an external shift.
!! In the case shift \f$ se=-1 \f$ the derivative \f$ \chi'(t) \f$ is 0 for
!! \f$ t \in [0,1-\eta] \cup [1, \infty) \f$ and varies in
!! \f$ [1-\eta, 1] \f$ which is an internal shift.
!! In the case shift \f$ se=0 \f$ the derivative \f$ \chi'(t) \f$ is 0 for
!! \f$ t \in [0,1-\eta/2] \cup [1+\eta/2, \infty) \f$ and
!! varies in \f$ [1-\eta/2, 1+\eta/2] \f$ which is a centered shift.
!!
!! @param[in] t: real non-negative input value
!! @param[in] se: shift
!! @param[in] eta: regularization parameter \f$ \eta \f$
real(dp) function dfsw(t, se, eta)
    ! Inputs
    real(dp), intent(in) :: t, se, eta
    ! Local variables
    real(dp) :: flow, x
    real(dp), parameter :: f30=30.0d0
    ! Apply shift:
    !   s =  0   =>   t - eta/2  [ CENTERED ]
    !   s =  1   =>   t - eta    [ EXTERIOR ]
    !   s = -1   =>   t          [ INTERIOR ]
    x = t - (se + 1.d0)*eta / 2.d0
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
    real(dp), intent(in) :: vwgrid((p+1)**2, ngrid)
    real(dp), intent(in) :: x(ngrid)
    ! Output
    real(dp), intent(out) :: xlm((p+1)**2)
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

!> Integrate against spherical harmonics by dgemv
!!
!! Integrate by Lebedev spherical quadrature by BLAS dgemv.
!!
!! @param[in] ngrid: Number of Lebedev grid points. `ngrid` > 0
!! @param[in] p: Maximal degree of spherical harmonics. `p` >= 0
!! @param[in] vwgrid: Values of spherical harmonics at Lebedev grid points,
!!      multiplied by weights of grid points. Dimension is ((p+1)**2, ngrid)
!! @param[in] isph: Index of given sphere. Used for debug purpose.
!! @param[in] x: Input values at grid points of the sphere. Dimension is
!!      (ngrid)
!! @param[out] xlm: Output spherical harmonics. Dimension is ((p+1)**2)
subroutine intrhs2(iprint, ngrid, p, vwgrid, isph, x, xlm)
    ! Inputs
    integer, intent(in) :: iprint, ngrid, p, isph
    real(dp), intent(in) :: vwgrid((p+1)**2, ngrid)
    real(dp), intent(in) :: x(ngrid)
    ! Output
    real(dp), intent(out) :: xlm((p+1)**2)
    ! Local variables
    integer :: nbasis
    nbasis = (p+1)**2
    ! Integrate by dgemv
    call dgemv('N', nbasis, ngrid, one, vwgrid, nbasis, x, 1, zero, xlm, 1)
    ! Printing (these functions are not implemented yet)
    if (iprint .ge. 5) then
        call ptcart('pot', ngrid, 1, isph, x)
        call prtsph('vlm', (p+1)**2, p, 1, isph, xlm)
    end if
end subroutine intrhs2

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
    real(dp), dimension(3),        intent(in)    :: x
    real(dp), dimension((dd_data % lmax+1)**2),   intent(inout) :: basloc, vplm
    real(dp), dimension(3,(dd_data % lmax +1)**2), intent(inout) :: dbsloc
    real(dp), dimension(dd_data % lmax+1),   intent(inout) :: vcos, vsin
    integer :: l, m, ind, VC, VS
    real(dp)  :: cthe, sthe, cphi, sphi, plm, fln, pp1, pm1, pp
    real(dp)  :: et(3), ep(3)
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
real(dp) function intmlp( dd_data, t, sigma, basloc )
!  
      implicit none
      type(dd_data_type) :: dd_data
      real(dp), intent(in) :: t
      real(dp), dimension((dd_data % lmax+1)**2), intent(in) :: sigma, basloc
!
      integer :: l, ind
      real(dp)  :: tt, ss, fac
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
      real(dp), dimension(dd_data % ncav),       intent(in)  :: phi
      real(dp), dimension(dd_data % ngrid, dd_data % nsph), intent(out) :: g
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
      real(dp), dimension((dd_data % lmax+1)**2), intent(in)    :: u
      real(dp),                    intent(inout) :: unorm
!
      integer :: l, m, ind
      real(dp)  :: fac
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
real(dp) function hnorm(dd_data, x)
    type(dd_data_type), intent(in) :: dd_data
      real(dp),  dimension(dd_data % nbasis, dd_data % nsph), intent(in) :: x
!
      integer                                     :: isph, istatus
      real(dp)                                      :: vrms, vmax
      real(dp), allocatable                         :: u(:)
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
      real(dp),  dimension(n), intent(in)    :: v
      real(dp),                intent(inout) :: vrms, vmax
!
      integer :: i
      real(dp), parameter :: zero=0.0d0
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
      real(dp), dimension(dd_data % ngrid, dd_data % nsph), intent(in)    :: xi
      real(dp), dimension((dd_data % lmax+1)**2),     intent(inout) :: vlm
      real(dp), dimension((dd_data % lmax+1)**2),     intent(inout) :: basloc, vplm
      real(dp), dimension(dd_data % lmax+1),     intent(inout) :: vcos, vsin
!
      integer :: ij, jsph, ig, l, ind, m
      real(dp)  :: vji(3), vvji, tji, sji(3), xji, oji, fac, ffac, t
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
          if ( tji.lt.( one + (dd_data % se+one)/two*dd_data % eta ) ) then
!                  
!           compute s_n^ji
            sji = vji/vvji
!
!           compute \chi( t_n^ji )
            xji = fsw( tji, dd_data % se, dd_data % eta )
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
      real(dp), dimension((dd_data % lmax+1)**2, dd_data % nsph), intent(in)    :: sigma
      real(dp), dimension(dd_data % ngrid),       intent(inout) :: pot
      real(dp), dimension((dd_data % lmax+1)**2),      intent(inout) :: basloc
      real(dp), dimension((dd_data % lmax+1)**2),      intent(inout) :: vplm
      real(dp), dimension(dd_data % lmax+1),      intent(inout) :: vcos
      real(dp), dimension(dd_data % lmax+1),      intent(inout) :: vsin
!
      integer :: its, ij, jsph
      real(dp)  :: vij(3), sij(3)
      real(dp)  :: vvij, tij, xij, oij, stslm, stslm2, stslm3, &
          & thigh, rho, ctheta, stheta, cphi, sphi
!
!------------------------------------------------------------------------
      thigh = one + pt5*(dd_data % se + one)*dd_data % eta
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
              xij = fsw( tij, dd_data % se, dd_data % eta )
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
       real(dp), dimension(dd_data % nbasis, dd_data % nsph), intent(in)    :: s
       real(dp), dimension(dd_data % ncav),      intent(inout) :: xi
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
       real(dp), dimension(dd_data % nbasis, dd_data % nsph), intent(in)    :: s
       real(dp), dimension(dd_data % ncav),      intent(inout) :: zeta
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
    real(dp)  :: f, f1
    integer, allocatable :: ipiv(:)
    real(dp),  allocatable :: work(:)
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

!> Accumulate a multipole expansion induced by a particle of a given charge
!!
!! Computes the following sums:
!! \f[
!!      \forall \ell=0, \ldots, p, \quad \forall m=-\ell, \ldots, \ell : \quad
!!      M_\ell^m = \beta M_\ell^m + \frac{q \|c\|^\ell}{r^{\ell+1}}
!!      Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of output harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \beta \f$ is a scaling factor, \f$ q \f$ and \f$
!! c \f$ are a charge and coordinates of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{\ell+1} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^{-\ell-1} \f$.
!!
!! @param[in] c: Radius-vector from the particle to the center of harmonics
!! @param[in] src_q: Charge of the source particle
!! @param[in] dst_r: Radius of output multipole spherical harmonics
!! @param[in] p: Maximal degree of output multipole spherical harmonics
!! @param[in] vscales: Normalization constants for spherical harmonics
!! @param[in] beta: Scaling factor for `dst_m`
!! @param[inout] dst_m: Multipole coefficients
!!
!! @sa fmm_m2p
subroutine fmm_p2m(c, src_q, dst_r, p, vscales, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales((p+1)**2), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)**2)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(dp) :: vylm((p+1)**2), vplm((p+1)**2), t, rcoef
    integer :: n, ind
    ! Get radius and values of spherical harmonics
    call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
        & vcos, vsin)
    ! Harmonics are available only if rho > 0
    if (rho .ne. zero) then
        rcoef = rho / dst_r
        t = src_q / dst_r
        ! Ignore input `m` in case of zero scaling factor
        if (beta .eq. zero) then
            do n = 0, p
                ind = n*n + n + 1
                dst_m(ind-n:ind+n) = t * vylm(ind-n:ind+n)
                t = t * rcoef
            end do
        ! Update `m` otherwise
        else
            do n = 0, p
                ind = n*n + n + 1
                dst_m(ind-n:ind+n) = beta*dst_m(ind-n:ind+n) + &
                    & t*vylm(ind-n:ind+n)
                t = t * rcoef
            end do
        end if
    ! Naive case of rho = 0
    else
        ! Ignore input `m` in case of zero scaling factor
        if (beta .eq. zero) then
            dst_m(1) = src_q / dst_r / sqrt4pi
            dst_m(2:) = zero
        ! Update `m` otherwise
        else
            dst_m(1) = beta*dst_m(1) + src_q/dst_r/sqrt4pi
            dst_m(2:) = beta * dst_m(2:)
        end if
    end if
end subroutine fmm_p2m

!> Accumulate potential, induced by multipole spherical harmonics
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{r}{\|c\|} \right)^{\ell+1} \sum_{m=-\ell}^\ell
!!      M_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of multipole basis functions
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$. Dimension
!!      is `(p+1)**2`
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole coefficients. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `v`
!! @param[inout] v: Value of induced potential
subroutine fmm_m2p(c, r, p, vscales, alpha, src_m, beta, v)
    ! Inputs
    real(dp), intent(in) :: c(3), r, vscales((p+1)*(p+1)), alpha, &
        & src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: v
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(dp) :: vylm((p+1)**2), vplm((p+1)**2), rcoef, t, tmp
    integer :: n, ind
    ! Scale output
    if (beta .eq. zero) then
        v = zero
    else
        v = beta * v
    end if
    ! In case of zero alpha nothing else is required no matter what is the
    ! value of the induced potential
    if (alpha .eq. zero) then
        return
    end if
    ! Get radius and values of spherical harmonics
    call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
        & vcos, vsin)
    ! In case of a singularity (rho=zero) induced potential is infinite and is
    ! not taken into account.
    if (rho .eq. zero) then
        return
    end if
    ! Compute the actual induced potential
    rcoef = r / rho
    t = alpha
    do n = 0, p
        t = t * rcoef
        ind = n*n + n + 1
        tmp = dot_product(vylm(ind-n:ind+n), src_m(ind-n:ind+n))
        ! Here vscales(ind)**2 is 4*pi/sqrt(2n+1)
        v = v + t*tmp/vscales(ind)**2
    end do
end subroutine fmm_m2p

!> Accumulate a local expansion induced by a particle of a given charge
!!
!! Computes the following sums:
!! \f[
!!      \forall \ell=0, \ldots, p, \quad \forall m=-\ell, \ldots, \ell : \quad
!!      M_\ell^m = \beta M_\ell^m + \frac{q \|c\|^\ell}{r^{\ell+1}}
!!      Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ M \f$ is a vector of coefficients of output harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \beta \f$ is a scaling factor, \f$ q \f$ and \f$
!! c \f$ are a charge and coordinates of a particle.
!!
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{\ell+1} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^{-\ell-1} \f$.
!!
!! @param[in] c: Radius-vector from the particle to the center of harmonics
!! @param[in] src_q: Charge of the source particle
!! @param[in] dst_r: Radius of output local spherical harmonics
!! @param[in] p: Maximal degree of output local spherical harmonics
!! @param[in] vscales: Normalization constants for spherical harmonics
!! @param[in] beta: Scaling factor for `dst_l`
!! @param[inout] dst_l: Local coefficients
!!
!! @sa fmm_l2p
subroutine fmm_p2l(c, src_q, dst_r, p, vscales, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_q, dst_r, vscales((p+1)**2), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)**2)
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(dp) :: vylm((p+1)**2), vplm((p+1)**2), t, rcoef
    integer :: n, ind
    ! Get radius and values of spherical harmonics
    call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
        & vcos, vsin)
    ! Local expansion represents potential from the outside of the sphere of
    ! the given radius `dst_r`. In case `rho=zero` source particle is inside of
    ! the sphere no matter what radius `r` is, so we just ignore such a case.
    if (rho .eq. zero) then
        if (beta .eq. zero) then
            dst_l = zero
        else
            dst_l = beta * dst_l
        end if
        return
    end if
    rcoef = dst_r / rho
    t = src_q / rho
    ! Ignore input `m` in case of zero scaling factor
    if (beta .eq. zero) then
        do n = 0, p
            ind = n*n + n + 1
            dst_l(ind-n:ind+n) = t * vylm(ind-n:ind+n)
            t = t * rcoef
        end do
    ! Update `m` otherwise
    else
        do n = 0, p
            ind = n*n + n + 1
            dst_l(ind-n:ind+n) = beta*dst_l(ind-n:ind+n) + t*vylm(ind-n:ind+n)
            t = t * rcoef
        end do
    end if
end subroutine fmm_p2l

!> Accumulate potential, induced by local spherical harmonics
!!
!! Computes the following sum:
!! \f[
!!      v = \beta v + \alpha \sum_{\ell=0}^p \frac{4\pi}{\sqrt{2\ell+1}}
!!      \left( \frac{\|c\|}{r} \right)^\ell \sum_{m=-\ell}^\ell
!!      L_\ell^m Y_\ell^m \left( \frac{c}{\|c\|} \right),
!! \f]
!! where \f$ L \f$ is a vector of coefficients of input harmonics of
!! a degree up to \f$ p \f$ inclusively with a convergence radius \f$ r \f$
!! located at the origin, \f$ \alpha \f$ and \f$ \beta \f$ are scaling factors
!! and \f$ c \f$ is a location of a particle.
!! Based on normalized real spherical harmonics \f$ Y_\ell^m \f$, scaled by \f$
!! r^{-\ell} \f$. It means corresponding coefficients are simply scaled by an
!! additional factor \f$ r^\ell \f$.
!!
!! @param[in] c: Coordinates of a particle (relative to center of harmonics)
!! @param[in] src_r: Radius of spherical harmonics
!! @param[in] p: Maximal degree of local basis functions
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$. Dimension
!!      is `(p+1)**2`
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Local coefficients. Dimension is `(p+1)**2`
!! @param[in] beta: Scalar multiplier for `dst_v`
!! @param[inout] dst_v: Value of induced potential
subroutine fmm_l2p(c, src_r, p, vscales, alpha, src_l, beta, dst_v)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, vscales((p+1)*(p+1)), alpha, &
        & src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_v
    ! Local variables
    real(dp) :: tmp, tmp1, tmp2
    real(dp) :: rho, t, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1)
    real(dp) :: vplm((p+1)**2), vylm((p+1)**2), rcoef
    integer :: n, k, ind
    ! Scale output
    if (beta .eq. zero) then
        dst_v = zero
    else
        dst_v = beta * dst_v
    end if
    ! In case of zero alpha nothing else is required no matter what is the
    ! value of the induced potential
    if (alpha .eq. zero) then
        return
    end if
    ! Get radius and values of spherical harmonics
    call ylmbas(c, rho, ctheta, stheta, cphi, sphi, p, vscales, vylm, vplm, &
        & vcos, vsin)
    ! In case `rho=zero` values of spherical harmonics make no sense and only
    ! spherical harmonic of degree 0 shall be taken into account
    if (rho .eq. zero) then
        dst_v = dst_v + alpha*src_l(1)/vscales(1)
    ! Compute the actual induced potential otherwise
    else
        rcoef = rho / src_r
        t = alpha
        do n = 0, p
            ind = n*n + n + 1
            tmp = dot_product(vylm(ind-n:ind+n), src_l(ind-n:ind+n))
            dst_v = dst_v + t*tmp/(vscales(ind)**2)
            t = t * rcoef
        end do
    end if
end subroutine fmm_l2p

!> Transform coefficients of spherical harmonics to a new cartesion system
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$ R \f$
!! is a '(p+1)**2'-by-'(p+1)**2' matrix that transforms coefficients of the
!! source spherical harmonics to the destination spherical harmonics while
!! preserving values of any function on a unit sphere defined as a linear
!! combination of the source spherical harmonics.
!!
!! Input 3-by-3 matrix `r1` must be an orthogonal matrix \f$ R_1 \f$ of a
!! transform of new cartesion coordinates \f$ (\widetilde{y}, \widetilde{z},
!! \widetilde{x}) \f$ into initial cartesian coordinates \f$ (y, z, x) \f$.
!! This is due to the following equalities:
!! \f{align}{
!!      Y_1^{-1} (\theta, \phi) &= \sqrt{\frac{3}{4\pi}} \sin \theta \sin \phi
!!      = \sqrt{\frac{3}{4\pi}} y, \\ Y_1^0 (\theta, \phi) &=
!!      \sqrt{\frac{3}{4\pi}} \cos \theta = \sqrt{\frac{3}{4\pi}} z, \\
!!      Y_1^1 (\theta, \phi) &= \sqrt{\frac{3}{4\pi}} \sin \theta \cos \phi =
!!      \sqrt{\frac{3}{4\pi}} x.
!! \f}
!! So, to find a column-vector \f$ \widetilde{c} \f$ of coefficients of
!! spherical harmonics \f$ Y_1^{-1}, Y_1^0 \f$ and \f$ Y_1^1 \f$ in a new
!! system of coordinates \f$ (\widetilde{y}, \widetilde{z}, \widetilde{x}) \f$
!! the following system needs to be solved:
!! \f[
!!      \widetilde{c}^\top \cdot \begin{bmatrix} Y_1^{-1}
!!      (\widetilde{\theta}, \widetilde{\phi}) \\ Y_1^0 (\widetilde{\theta},
!!      \widetilde{\phi}) \\ Y_1^1 (\widetilde{\theta}, \widetilde{\phi})
!!      \end{bmatrix} = c ^\top \cdot \begin{bmatrix} Y_1^{-1} (\theta, \phi)
!!      \\ Y_1^0 (\theta, \phi) \\ Y_1^1 (\theta, \phi) \end{bmatrix}.
!! \f]
!! The solution has the following obvious form:
!! \f[
!!      \widetilde{c} = R_1^\top c.
!! \f]
!!
!! Translation of spherical harmonics of all other degrees is computed
!! recursively as is described in the following source:
!!      @cite ir-realharms-1996
!!      @cite ir-realharms-1998
!!
!!
!! @param[in] p: Maximum degree of spherical harmonics
!! @param[in] r1: Transformation from new to old cartesian coordinates
!! @param[in] alpha: Scalar multipler for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[inout] dst: Coefficients of transformed spherical harmonics
!!
!! @sa fmm_sph_transform_get_mat, fmm_sph_transform_use_mat
subroutine fmm_sph_transform(p, r1, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: r1(-1:1, -1:1), alpha, src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)*(p+1))
    ! Local variables
    real(dp) :: tmp((p+1)*(p+1))
    ! Take care of beta param
    if (beta .eq. zero) then
        call fmm_sph_transform_out(p, r1, alpha*src, dst)
    else
        call fmm_sph_transform_out(p, r1, alpha*src, tmp)
        dst = beta*dst + tmp
    end if
end subroutine fmm_sph_transform

!> Transform spherical harmonics to a new cartesion system of coordinates
!!
!! This function implements @ref fmm_sph_transform with predefined values of
!! parameters \p alpha=one and \p beta=zero.
!! 
!!
!! @param[in] p: Maximum degree of spherical harmonics
!! @param[in] r1: Transformation from new to old cartesian coordinates
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[out] dst: Coefficients of transformed spherical harmonics
!!
!! @sa fmm_sph_transform_get_mat, fmm_sph_transform_use_mat
subroutine fmm_sph_transform_out(p, r1, src, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: r1(-1:1, -1:1), src((p+1)*(p+1))
    ! Output
    real(dp), intent(out) :: dst((p+1)*(p+1))
    ! Local variables
    real(dp) :: u, v, w, r_prev(-p:p, -p:p), r(-p:p, -p:p)
    real(dp) :: scal_uvw_m(-p:p), scal_u_n(-p:p), scal_v_n(-p:p)
    real(dp) :: scal_w_n(-p:p)
    integer :: l, m, n, ind
    ! Compute rotations/reflections
    ! l = 0
    dst(1) = src(1)
    if (p .eq. 0) then
        return
    end if
    ! l = 1
    do m = -1, 1
        dst(3+m) = 0
        do n = -1, 1
            dst(3+m) = dst(3+m) + r1(n, m)*src(3+n)
        end do
    end do
    if (p .eq. 1) then
        return
    end if
    ! l = 2
    r(2, 2) = (r1(1, 1)*r1(1, 1) - r1(1, -1)*r1(1, -1) - &
        & r1(-1, 1)*r1(-1, 1) + r1(-1, -1)*r1(-1, -1)) / 2
    r(2, 1) = r1(1, 1)*r1(1, 0) - r1(-1, 1)*r1(-1, 0)
    r(2, 0) = sqrt(dble(3)) / 2 * (r1(1, 0)*r1(1, 0) - r1(-1, 0)*r1(-1, 0))
    r(2, -1) = r1(1, -1)*r1(1, 0) - r1(-1, -1)*r1(-1, 0)
    r(2, -2) = r1(1, 1)*r1(1, -1) - r1(-1, 1)*r1(-1, -1)
    r(1, 2) = r1(1, 1)*r1(0, 1) - r1(1, -1)*r1(0, -1)
    r(1, 1) = r1(1, 1)*r1(0, 0) + r1(1, 0)*r1(0, 1)
    r(1, 0) = sqrt(dble(3)) * r1(1, 0) * r1(0, 0)
    r(1, -1) = r1(1, -1)*r1(0, 0) + r1(1, 0)*r1(0, -1)
    r(1, -2) = r1(1, 1)*r1(0, -1) + r1(1, -1)*r1(0, 1)
    r(0, 2) = sqrt(dble(3)) / 2 * (r1(0, 1)*r1(0, 1) - r1(0, -1)*r1(0, -1))
    r(0, 1) = sqrt(dble(3)) * r1(0, 1) * r1(0, 0)
    r(0, 0) = (3*r1(0, 0)*r1(0, 0)-1) / 2
    r(0, -1) = sqrt(dble(3)) * r1(0, -1) * r1(0, 0)
    r(0, -2) = sqrt(dble(3)) * r1(0, 1) * r1(0, -1)
    r(-1, 2) = r1(-1, 1)*r1(0, 1) - r1(-1, -1)*r1(0, -1)
    r(-1, 1) = r1(-1, 1)*r1(0, 0) + r1(-1, 0)*r1(0, 1)
    r(-1, 0) = sqrt(dble(3)) * r1(-1, 0) * r1(0, 0)
    r(-1, -1) = r1(-1, -1)*r1(0, 0) + r1(0, -1)*r1(-1, 0)
    r(-1, -2) = r1(-1, 1)*r1(0, -1) + r1(-1, -1)*r1(0, 1)
    r(-2, 2) = r1(1, 1)*r1(-1, 1) - r1(1, -1)*r1(-1, -1)
    r(-2, 1) = r1(1, 1)*r1(-1, 0) + r1(1, 0)*r1(-1, 1)
    r(-2, 0) = sqrt(dble(3)) * r1(1, 0) * r1(-1, 0)
    r(-2, -1) = r1(1, -1)*r1(-1, 0) + r1(1, 0)*r1(-1, -1)
    r(-2, -2) = r1(1, 1)*r1(-1, -1) + r1(1, -1)*r1(-1, 1)
    do m = -2, 2
        dst(7+m) = 0
        do n = -2, 2
            dst(7+m) = dst(7+m) + src(7+n)*r(n, m)
        end do
    end do
    ! l > 2
    do l = 3, p
        ! Prepare previous rotation matrix
        r_prev(1-l:l-1, 1-l:l-1) = r(1-l:l-1, 1-l:l-1)
        ! Prepare scalar factors
        scal_uvw_m(0) = l
        do m = 1, l-1
        u = sqrt(dble(l*l-m*m))
        scal_uvw_m(m) = u
        scal_uvw_m(-m) = u
        end do
        u = 2 * l
        u = sqrt(dble(u*(u-1)))
        scal_uvw_m(l) = u
        scal_uvw_m(-l) = u
        scal_u_n(0) = l
        scal_v_n(0) = -sqrt(dble(l*(l-1))) / sqrt2
        scal_w_n(0) = 0
        do n = 1, l-2
        u = sqrt(dble(l*l-n*n))
        scal_u_n(n) = u
        scal_u_n(-n) = u
        v = l + n
        v = sqrt(v*(v-1)) / 2
        scal_v_n(n) = v
        scal_v_n(-n) = v
        w = l - n
        w = -sqrt(w*(w-1)) / 2
        scal_w_n(n) = w
        scal_w_n(-n) = w
        end do
        u = sqrt(dble(2*l-1))
        scal_u_n(l-1) = u
        scal_u_n(1-l) = u
        scal_u_n(l) = 0
        scal_u_n(-l) = 0
        v = sqrt(dble((2*l-1)*(2*l-2))) / 2
        scal_v_n(l-1) = v
        scal_v_n(1-l) = v
        v = sqrt(dble(2*l*(2*l-1))) / 2
        scal_v_n(l) = v
        scal_v_n(-l) = v
        scal_w_n(l-1) = 0
        scal_w_n(l) = 0
        scal_w_n(-l) = 0
        scal_w_n(1-l) = 0
        ind = l*l + l + 1
        ! m = l, n = l
        v = r1(1, 1)*r_prev(l-1, l-1) - r1(1, -1)*r_prev(l-1, 1-l) - &
            & r1(-1, 1)*r_prev(1-l, l-1) + r1(-1, -1)*r_prev(1-l, 1-l)
        r(l, l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind+l) = src(ind+l) * r(l, l)
        ! m = l, n = -l
        v = r1(1, 1)*r_prev(1-l, l-1) - r1(1, -1)*r_prev(1-l, 1-l) + &
            & r1(-1, 1)*r_prev(l-1, l-1) - r1(-1, -1)*r_prev(l-1, 1-l)
        r(-l, l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind-l)*r(-l, l)
        ! m = l, n = l-1
        u = r1(0, 1)*r_prev(l-1, l-1) - r1(0, -1)*r_prev(l-1, 1-l)
        v = r1(1, 1)*r_prev(l-2, l-1) - r1(1, -1)*r_prev(l-2, 1-l) - &
            & r1(-1, 1)*r_prev(2-l, l-1) + r1(-1, -1)*r_prev(2-l, 1-l)
        r(l-1, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind+l-1)*r(l-1, l)
        ! m = l, n = 1-l
        u = r1(0, 1)*r_prev(1-l, l-1) - r1(0, -1)*r_prev(1-l, 1-l)
        v = r1(1, 1)*r_prev(2-l, l-1) - r1(1, -1)*r_prev(2-l, 1-l) + &
            & r1(-1, 1)*r_prev(l-2, l-1) - r1(-1, -1)*r_prev(l-2, 1-l)
        r(1-l, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind+1-l)*r(1-l, l)
        ! m = l, n = 1
        u = r1(0, 1)*r_prev(1, l-1) - r1(0, -1)*r_prev(1, 1-l)
        v = r1(1, 1)*r_prev(0, l-1) - r1(1, -1)*r_prev(0, 1-l)
        w = r1(1, 1)*r_prev(2, l-1) - r1(1, -1)*r_prev(2, 1-l) + &
            & r1(-1, 1)*r_prev(-2, l-1) - r1(-1, -1)*r_prev(-2, 1-l)
        r(1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(1, l) = r(1, l) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind+1)*r(1, l)
        ! m = l, n = -1
        u = r1(0, 1)*r_prev(-1, l-1) - r1(0, -1)*r_prev(-1, 1-l)
        v = r1(-1, 1)*r_prev(0, l-1) - r1(-1, -1)*r_prev(0, 1-l)
        w = r1(1, 1)*r_prev(-2, l-1) - r1(1, -1)*r_prev(-2, 1-l) - &
            & r1(-1, 1)*r_prev(2, l-1) + r1(-1, -1)*r_prev(2, 1-l)
        r(-1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(-1, l) = r(-1, l) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind-1)*r(-1, l)
        ! m = l, n = 0
        u = r1(0, 1)*r_prev(0, l-1) - r1(0, -1)*r_prev(0, 1-l)
        v = r1(1, 1)*r_prev(1, l-1) - r1(1, -1)*r_prev(1, 1-l) + &
            & r1(-1, 1)*r_prev(-1, l-1) - r1(-1, -1)*r_prev(-1, 1-l)
        r(0, l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind)*r(0, l)
        ! m = l, n = 2-l..l-2, n != -1,0,1
        do n = 2, l-2
            u = r1(0, 1)*r_prev(n, l-1) - r1(0, -1)*r_prev(n, 1-l)
            v = r1(1, 1)*r_prev(n-1, l-1) - r1(1, -1)*r_prev(n-1, 1-l) - &
                & r1(-1, 1)*r_prev(1-n, l-1) + r1(-1, -1)*r_prev(1-n, 1-l)
            w = r1(1, 1)*r_prev(n+1, l-1) - r1(1, -1)*r_prev(n+1, 1-l) + &
                & r1(-1, 1)*r_prev(-n-1, l-1) - r1(-1, -1)*r_prev(-n-1, 1-l)
            r(n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(n, l) = r(n, l) / scal_uvw_m(l)
            dst(ind+l) = dst(ind+l) + src(ind+n)*r(n, l)
            u = r1(0, 1)*r_prev(-n, l-1) - r1(0, -1)*r_prev(-n, 1-l)
            v = r1(1, 1)*r_prev(1-n, l-1) - r1(1, -1)*r_prev(1-n, 1-l) + &
                & r1(-1, 1)*r_prev(n-1, l-1) - r1(-1, -1)*r_prev(n-1, 1-l)
            w = r1(1, 1)*r_prev(-n-1, l-1) - r1(1, -1)*r_prev(-n-1, 1-l) - &
                & r1(-1, 1)*r_prev(n+1, l-1) + r1(-1, -1)*r_prev(n+1, 1-l)
            r(-n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(-n, l) = r(-n, l) / scal_uvw_m(l)
            dst(ind+l) = dst(ind+l) + src(ind-n)*r(-n, l)
        end do
        ! m = -l, n = l
        v = r1(1, 1)*r_prev(l-1, 1-l) + r1(1, -1)*r_prev(l-1, l-1) - &
            & r1(-1, 1)*r_prev(1-l, 1-l) - r1(-1, -1)*r_prev(1-l, l-1)
        r(l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind-l) = src(ind+l)*r(l, -l)
        ! m = -l, n = -l
        v = r1(1, 1)*r_prev(1-l, 1-l) + r1(1, -1)*r_prev(1-l, l-1) + &
            & r1(-1, 1)*r_prev(l-1, 1-l) + r1(-1, -1)*r_prev(l-1, l-1)
        r(-l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind-l)*r(-l, -l)
        ! m = -l, n = l-1
        u = r1(0, 1)*r_prev(l-1, 1-l) + r1(0, -1)*r_prev(l-1, l-1)
        v = r1(1, 1)*r_prev(l-2, 1-l) + r1(1, -1)*r_prev(l-2, l-1) - &
            & r1(-1, 1)*r_prev(2-l, 1-l) - r1(-1, -1)*r_prev(2-l, l-1)
        r(l-1, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind+l-1)*r(l-1, -l)
        ! m = -l, n = 1-l
        u = r1(0, 1)*r_prev(1-l, 1-l) + r1(0, -1)*r_prev(1-l, l-1)
        v = r1(1, 1)*r_prev(2-l, 1-l) + r1(1, -1)*r_prev(2-l, l-1) + &
            & r1(-1, 1)*r_prev(l-2, 1-l) + r1(-1, -1)*r_prev(l-2, l-1)
        r(1-l, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind+1-l)*r(1-l, -l)
        ! m = -l, n = 1
        u = r1(0, 1)*r_prev(1, 1-l) + r1(0, -1)*r_prev(1, l-1)
        v = r1(1, 1)*r_prev(0, 1-l) + r1(1, -1)*r_prev(0, l-1)
        w = r1(1, 1)*r_prev(2, 1-l) + r1(1, -1)*r_prev(2, l-1) + &
            & r1(-1, 1)*r_prev(-2, 1-l) + r1(-1, -1)*r_prev(-2, l-1)
        r(1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(1, -l) = r(1, -l) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind+1)*r(1, -l)
        ! m = -l, n = -1
        u = r1(0, 1)*r_prev(-1, 1-l) + r1(0, -1)*r_prev(-1, l-1)
        v = r1(-1, 1)*r_prev(0, 1-l) + r1(-1, -1)*r_prev(0, l-1)
        w = r1(1, 1)*r_prev(-2, 1-l) + r1(1, -1)*r_prev(-2, l-1) - &
            & r1(-1, 1)*r_prev(2, 1-l) - r1(-1, -1)*r_prev(2, l-1)
        r(-1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(-1, -l) = r(-1, -l) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind-1)*r(-1, -l)
        ! m = -l, n = 0
        u = r1(0, 1)*r_prev(0, 1-l) + r1(0, -1)*r_prev(0, l-1)
        v = r1(1, 1)*r_prev(1, 1-l) + r1(1, -1)*r_prev(1, l-1) + &
            & r1(-1, 1)*r_prev(-1, 1-l) + r1(-1, -1)*r_prev(-1, l-1)
        r(0, -l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind)*r(0, -l)
        ! m = -l, n = 2-l..l-2, n != -1,0,1
        do n = 2, l-2
            u = r1(0, 1)*r_prev(n, 1-l) + r1(0, -1)*r_prev(n, l-1)
            v = r1(1, 1)*r_prev(n-1, 1-l) + r1(1, -1)*r_prev(n-1, l-1) - &
                & r1(-1, 1)*r_prev(1-n, 1-l) - r1(-1, -1)*r_prev(1-n, l-1)
            w = r1(1, 1)*r_prev(n+1, 1-l) + r1(1, -1)*r_prev(n+1, l-1) + &
                & r1(-1, 1)*r_prev(-n-1, 1-l) + r1(-1, -1)*r_prev(-n-1, l-1)
            r(n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(n, -l) = r(n, -l) / scal_uvw_m(l)
            dst(ind-l) = dst(ind-l) + src(ind+n)*r(n, -l)
            u = r1(0, 1)*r_prev(-n, 1-l) + r1(0, -1)*r_prev(-n, l-1)
            v = r1(1, 1)*r_prev(1-n, 1-l) + r1(1, -1)*r_prev(1-n, l-1) + &
                & r1(-1, 1)*r_prev(n-1, 1-l) + r1(-1, -1)*r_prev(n-1, l-1)
            w = r1(1, 1)*r_prev(-n-1, 1-l) + r1(1, -1)*r_prev(-n-1, l-1) - &
                & r1(-1, 1)*r_prev(n+1, 1-l) - r1(-1, -1)*r_prev(n+1, l-1)
            r(-n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(-n, -l) = r(-n, -l) / scal_uvw_m(l)
            dst(ind-l) = dst(ind-l) + src(ind-n)*r(-n, -l)
        end do
        ! Now deal with m=1-l..l-1
        do m = 1-l, l-1
            ! n = l
            v = r1(1, 0)*r_prev(l-1, m) - r1(-1, 0)*r_prev(1-l, m)
            r(l, m) = v * scal_v_n(l) / scal_uvw_m(m)
            dst(ind+m) = src(ind+l) * r(l, m)
            ! n = -l
            v = r1(1, 0)*r_prev(1-l, m) + r1(-1, 0)*r_prev(l-1, m)
            r(-l, m) = v * scal_v_n(l) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind-l)*r(-l, m)
            ! n = l-1
            u = r1(0, 0) * r_prev(l-1, m)
            v = r1(1, 0)*r_prev(l-2, m) - r1(-1, 0)*r_prev(2-l, m)
            r(l-1, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
            r(l-1, m) = r(l-1, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind+l-1)*r(l-1, m)
            ! n = 1-l
            u = r1(0, 0) * r_prev(1-l, m)
            v = r1(1, 0)*r_prev(2-l, m) + r1(-1, 0)*r_prev(l-2, m)
            r(1-l, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
            r(1-l, m) = r(1-l, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind+1-l)*r(1-l, m)
            ! n = 0
            u = r1(0, 0) * r_prev(0, m)
            v = r1(1, 0)*r_prev(1, m) + r1(-1, 0)*r_prev(-1, m)
            r(0, m) = u*scal_u_n(0) + v*scal_v_n(0)
            r(0, m) = r(0, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind)*r(0, m)
            ! n = 1
            u = r1(0, 0) * r_prev(1, m)
            v = r1(1, 0) * r_prev(0, m)
            w = r1(1, 0)*r_prev(2, m) + r1(-1, 0)*r_prev(-2, m)
            r(1, m) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(1, m) = r(1, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind+1)*r(1, m)
            ! n = -1
            u = r1(0, 0) * r_prev(-1, m)
            v = r1(-1, 0) * r_prev(0, m)
            w = r1(1, 0)*r_prev(-2, m) - r1(-1, 0)*r_prev(2, m)
            r(-1, m) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(-1, m) = r(-1, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind-1)*r(-1, m)
            ! n = 2-l..l-2, n != -1,0,1
            do n = 2, l-2
                u = r1(0, 0) * r_prev(n, m)
                v = r1(1, 0)*r_prev(n-1, m) - r1(-1, 0)*r_prev(1-n, m)
                w = r1(1, 0)*r_prev(n+1, m) + r1(-1, 0)*r_prev(-1-n, m)
                r(n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(n, m) = r(n, m) / scal_uvw_m(m)
                dst(ind+m) = dst(ind+m) + src(ind+n)*r(n, m)
                u = r1(0, 0) * r_prev(-n, m)
                v = r1(1, 0)*r_prev(1-n, m) + r1(-1, 0)*r_prev(n-1, m)
                w = r1(1, 0)*r_prev(-n-1, m) - r1(-1, 0)*r_prev(n+1, m)
                r(-n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(-n, m) = r(-n, m) / scal_uvw_m(m)
                dst(ind+m) = dst(ind+m) + src(ind-n)*r(-n, m)
            end do
        end do
    end do
end subroutine fmm_sph_transform_out

!> Save matrix of transformation of spherical harmonics
!!
!! This procedure is based on @ref fmm_sph_transform_out but stores
!! corresponding transformation matrix in a data-efficient sparse way. For a
!! full description of the input parameter `r1` please take a look at @ref
!! fmm_sph_transform.
!!
!!
!! @param[in] p: Maximum degree of spherical harmonics
!! @param[in] r1: Transformation from new to old cartesian coordinates
!! @param[out] mat: Transformation matrix in a sparse storage
!!
!! @sa fmm_sph_transform, fmm_sph_transform_use_mat
subroutine fmm_sph_transform_get_mat(p, r1, mat)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: r1(-1:1, -1:1)
    ! Output
    real(dp), intent(out) :: mat((p+1)*(2*p+1)*(2*p+3)/3)
    ! Local variables
    real(dp) :: u, v, w, r_prev(-p:p, -p:p), r(-p:p, -p:p)
    real(dp) :: scal_uvw_m(-p:p), scal_u_n(-p:p), scal_v_n(-p:p)
    real(dp) :: scal_w_n(-p:p)
    integer :: l, m, n, ind
    ! Compute rotations/reflections
    ! l = 0
    !dst(1) = src(1)
    mat(1) = 1
    if (p .eq. 0) then
        return
    end if
    ! l = 1
    do m = -1, 1
        !dst(3+m) = 0
        ind = 3*m + 6
        do n = -1, 1
            !dst(3+m) = dst(3+m) + r1(n, m)*src(3+n)
            mat(ind+n) = r1(n, m)
        end do
    end do
    if (p .eq. 1) then
        return
    end if
    ! l = 2
    r(2, 2) = (r1(1, 1)*r1(1, 1) - r1(1, -1)*r1(1, -1) - &
        & r1(-1, 1)*r1(-1, 1) + r1(-1, -1)*r1(-1, -1)) / 2
    r(2, 1) = r1(1, 1)*r1(1, 0) - r1(-1, 1)*r1(-1, 0)
    r(2, 0) = sqrt(dble(3)) / 2 * (r1(1, 0)*r1(1, 0) - r1(-1, 0)*r1(-1, 0))
    r(2, -1) = r1(1, -1)*r1(1, 0) - r1(-1, -1)*r1(-1, 0)
    r(2, -2) = r1(1, 1)*r1(1, -1) - r1(-1, 1)*r1(-1, -1)
    r(1, 2) = r1(1, 1)*r1(0, 1) - r1(1, -1)*r1(0, -1)
    r(1, 1) = r1(1, 1)*r1(0, 0) + r1(1, 0)*r1(0, 1)
    r(1, 0) = sqrt(dble(3)) * r1(1, 0) * r1(0, 0)
    r(1, -1) = r1(1, -1)*r1(0, 0) + r1(1, 0)*r1(0, -1)
    r(1, -2) = r1(1, 1)*r1(0, -1) + r1(1, -1)*r1(0, 1)
    r(0, 2) = sqrt(dble(3)) / 2 * (r1(0, 1)*r1(0, 1) - r1(0, -1)*r1(0, -1))
    r(0, 1) = sqrt(dble(3)) * r1(0, 1) * r1(0, 0)
    r(0, 0) = (3*r1(0, 0)*r1(0, 0)-1) / 2
    r(0, -1) = sqrt(dble(3)) * r1(0, -1) * r1(0, 0)
    r(0, -2) = sqrt(dble(3)) * r1(0, 1) * r1(0, -1)
    r(-1, 2) = r1(-1, 1)*r1(0, 1) - r1(-1, -1)*r1(0, -1)
    r(-1, 1) = r1(-1, 1)*r1(0, 0) + r1(-1, 0)*r1(0, 1)
    r(-1, 0) = sqrt(dble(3)) * r1(-1, 0) * r1(0, 0)
    r(-1, -1) = r1(-1, -1)*r1(0, 0) + r1(0, -1)*r1(-1, 0)
    r(-1, -2) = r1(-1, 1)*r1(0, -1) + r1(-1, -1)*r1(0, 1)
    r(-2, 2) = r1(1, 1)*r1(-1, 1) - r1(1, -1)*r1(-1, -1)
    r(-2, 1) = r1(1, 1)*r1(-1, 0) + r1(1, 0)*r1(-1, 1)
    r(-2, 0) = sqrt(dble(3)) * r1(1, 0) * r1(-1, 0)
    r(-2, -1) = r1(1, -1)*r1(-1, 0) + r1(1, 0)*r1(-1, -1)
    r(-2, -2) = r1(1, 1)*r1(-1, -1) + r1(1, -1)*r1(-1, 1)
    do m = -2, 2
        !dst(7+m) = 0
        ind = 5*m + 23
        do n = -2, 2
            !dst(7+m) = dst(7+m) + src(7+n)*r(n, m)
            mat(ind+n) = r(n, m)
        end do
    end do
    ! l > 2
    do l = 3, p
        ! Prepare previous rotation matrix
        r_prev(1-l:l-1, 1-l:l-1) = r(1-l:l-1, 1-l:l-1)
        ! Prepare scalar factors
        scal_uvw_m(0) = l
        do m = 1, l-1
        u = sqrt(dble(l*l-m*m))
        scal_uvw_m(m) = u
        scal_uvw_m(-m) = u
        end do
        u = 2 * l
        u = sqrt(dble(u*(u-1)))
        scal_uvw_m(l) = u
        scal_uvw_m(-l) = u
        scal_u_n(0) = l
        scal_v_n(0) = -sqrt(dble(l*(l-1))) / sqrt2
        scal_w_n(0) = 0
        do n = 1, l-2
        u = sqrt(dble(l*l-n*n))
        scal_u_n(n) = u
        scal_u_n(-n) = u
        v = l + n
        v = sqrt(v*(v-1)) / 2
        scal_v_n(n) = v
        scal_v_n(-n) = v
        w = l - n
        w = -sqrt(w*(w-1)) / 2
        scal_w_n(n) = w
        scal_w_n(-n) = w
        end do
        u = sqrt(dble(2*l-1))
        scal_u_n(l-1) = u
        scal_u_n(1-l) = u
        scal_u_n(l) = 0
        scal_u_n(-l) = 0
        v = sqrt(dble((2*l-1)*(2*l-2))) / 2
        scal_v_n(l-1) = v
        scal_v_n(1-l) = v
        v = sqrt(dble(2*l*(2*l-1))) / 2
        scal_v_n(l) = v
        scal_v_n(-l) = v
        scal_w_n(l-1) = 0
        scal_w_n(l) = 0
        scal_w_n(-l) = 0
        scal_w_n(1-l) = 0
        ind = l*l + l + 1
        ! m = l, n = l
        v = r1(1, 1)*r_prev(l-1, l-1) - r1(1, -1)*r_prev(l-1, 1-l) - &
            & r1(-1, 1)*r_prev(1-l, l-1) + r1(-1, -1)*r_prev(1-l, 1-l)
        r(l, l) = v * scal_v_n(l) / scal_uvw_m(l)
        !dst(ind+l) = src(ind+l) * r(l, l)
        ! m = l, n = -l
        v = r1(1, 1)*r_prev(1-l, l-1) - r1(1, -1)*r_prev(1-l, 1-l) + &
            & r1(-1, 1)*r_prev(l-1, l-1) - r1(-1, -1)*r_prev(l-1, 1-l)
        r(-l, l) = v * scal_v_n(l) / scal_uvw_m(l)
        !dst(ind+l) = dst(ind+l) + src(ind-l)*r(-l, l)
        ! m = l, n = l-1
        u = r1(0, 1)*r_prev(l-1, l-1) - r1(0, -1)*r_prev(l-1, 1-l)
        v = r1(1, 1)*r_prev(l-2, l-1) - r1(1, -1)*r_prev(l-2, 1-l) - &
            & r1(-1, 1)*r_prev(2-l, l-1) + r1(-1, -1)*r_prev(2-l, 1-l)
        r(l-1, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        !dst(ind+l) = dst(ind+l) + src(ind+l-1)*r(l-1, l)
        ! m = l, n = 1-l
        u = r1(0, 1)*r_prev(1-l, l-1) - r1(0, -1)*r_prev(1-l, 1-l)
        v = r1(1, 1)*r_prev(2-l, l-1) - r1(1, -1)*r_prev(2-l, 1-l) + &
            & r1(-1, 1)*r_prev(l-2, l-1) - r1(-1, -1)*r_prev(l-2, 1-l)
        r(1-l, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        !dst(ind+l) = dst(ind+l) + src(ind+1-l)*r(1-l, l)
        ! m = l, n = 1
        u = r1(0, 1)*r_prev(1, l-1) - r1(0, -1)*r_prev(1, 1-l)
        v = r1(1, 1)*r_prev(0, l-1) - r1(1, -1)*r_prev(0, 1-l)
        w = r1(1, 1)*r_prev(2, l-1) - r1(1, -1)*r_prev(2, 1-l) + &
            & r1(-1, 1)*r_prev(-2, l-1) - r1(-1, -1)*r_prev(-2, 1-l)
        r(1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(1, l) = r(1, l) / scal_uvw_m(l)
        !dst(ind+l) = dst(ind+l) + src(ind+1)*r(1, l)
        ! m = l, n = -1
        u = r1(0, 1)*r_prev(-1, l-1) - r1(0, -1)*r_prev(-1, 1-l)
        v = r1(-1, 1)*r_prev(0, l-1) - r1(-1, -1)*r_prev(0, 1-l)
        w = r1(1, 1)*r_prev(-2, l-1) - r1(1, -1)*r_prev(-2, 1-l) - &
            & r1(-1, 1)*r_prev(2, l-1) + r1(-1, -1)*r_prev(2, 1-l)
        r(-1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(-1, l) = r(-1, l) / scal_uvw_m(l)
        !dst(ind+l) = dst(ind+l) + src(ind-1)*r(-1, l)
        ! m = l, n = 0
        u = r1(0, 1)*r_prev(0, l-1) - r1(0, -1)*r_prev(0, 1-l)
        v = r1(1, 1)*r_prev(1, l-1) - r1(1, -1)*r_prev(1, 1-l) + &
            & r1(-1, 1)*r_prev(-1, l-1) - r1(-1, -1)*r_prev(-1, 1-l)
        r(0, l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
        !dst(ind+l) = dst(ind+l) + src(ind)*r(0, l)
        ! m = l, n = 2-l..l-2, n != -1,0,1
        do n = 2, l-2
            u = r1(0, 1)*r_prev(n, l-1) - r1(0, -1)*r_prev(n, 1-l)
            v = r1(1, 1)*r_prev(n-1, l-1) - r1(1, -1)*r_prev(n-1, 1-l) - &
                & r1(-1, 1)*r_prev(1-n, l-1) + r1(-1, -1)*r_prev(1-n, 1-l)
            w = r1(1, 1)*r_prev(n+1, l-1) - r1(1, -1)*r_prev(n+1, 1-l) + &
                & r1(-1, 1)*r_prev(-n-1, l-1) - r1(-1, -1)*r_prev(-n-1, 1-l)
            r(n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(n, l) = r(n, l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind+n)*r(n, l)
            u = r1(0, 1)*r_prev(-n, l-1) - r1(0, -1)*r_prev(-n, 1-l)
            v = r1(1, 1)*r_prev(1-n, l-1) - r1(1, -1)*r_prev(1-n, 1-l) + &
                & r1(-1, 1)*r_prev(n-1, l-1) - r1(-1, -1)*r_prev(n-1, 1-l)
            w = r1(1, 1)*r_prev(-n-1, l-1) - r1(1, -1)*r_prev(-n-1, 1-l) - &
                & r1(-1, 1)*r_prev(n+1, l-1) + r1(-1, -1)*r_prev(n+1, 1-l)
            r(-n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(-n, l) = r(-n, l) / scal_uvw_m(l)
            !dst(ind+l) = dst(ind+l) + src(ind-n)*r(-n, l)
        end do
        ! m = -l, n = l
        v = r1(1, 1)*r_prev(l-1, 1-l) + r1(1, -1)*r_prev(l-1, l-1) - &
            & r1(-1, 1)*r_prev(1-l, 1-l) - r1(-1, -1)*r_prev(1-l, l-1)
        r(l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
        !dst(ind-l) = src(ind+l)*r(l, -l)
        ! m = -l, n = -l
        v = r1(1, 1)*r_prev(1-l, 1-l) + r1(1, -1)*r_prev(1-l, l-1) + &
            & r1(-1, 1)*r_prev(l-1, 1-l) + r1(-1, -1)*r_prev(l-1, l-1)
        r(-l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
        !dst(ind-l) = dst(ind-l) + src(ind-l)*r(-l, -l)
        ! m = -l, n = l-1
        u = r1(0, 1)*r_prev(l-1, 1-l) + r1(0, -1)*r_prev(l-1, l-1)
        v = r1(1, 1)*r_prev(l-2, 1-l) + r1(1, -1)*r_prev(l-2, l-1) - &
            & r1(-1, 1)*r_prev(2-l, 1-l) - r1(-1, -1)*r_prev(2-l, l-1)
        r(l-1, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        !dst(ind-l) = dst(ind-l) + src(ind+l-1)*r(l-1, -l)
        ! m = -l, n = 1-l
        u = r1(0, 1)*r_prev(1-l, 1-l) + r1(0, -1)*r_prev(1-l, l-1)
        v = r1(1, 1)*r_prev(2-l, 1-l) + r1(1, -1)*r_prev(2-l, l-1) + &
            & r1(-1, 1)*r_prev(l-2, 1-l) + r1(-1, -1)*r_prev(l-2, l-1)
        r(1-l, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        !dst(ind-l) = dst(ind-l) + src(ind+1-l)*r(1-l, -l)
        ! m = -l, n = 1
        u = r1(0, 1)*r_prev(1, 1-l) + r1(0, -1)*r_prev(1, l-1)
        v = r1(1, 1)*r_prev(0, 1-l) + r1(1, -1)*r_prev(0, l-1)
        w = r1(1, 1)*r_prev(2, 1-l) + r1(1, -1)*r_prev(2, l-1) + &
            & r1(-1, 1)*r_prev(-2, 1-l) + r1(-1, -1)*r_prev(-2, l-1)
        r(1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(1, -l) = r(1, -l) / scal_uvw_m(l)
        !dst(ind-l) = dst(ind-l) + src(ind+1)*r(1, -l)
        ! m = -l, n = -1
        u = r1(0, 1)*r_prev(-1, 1-l) + r1(0, -1)*r_prev(-1, l-1)
        v = r1(-1, 1)*r_prev(0, 1-l) + r1(-1, -1)*r_prev(0, l-1)
        w = r1(1, 1)*r_prev(-2, 1-l) + r1(1, -1)*r_prev(-2, l-1) - &
            & r1(-1, 1)*r_prev(2, 1-l) - r1(-1, -1)*r_prev(2, l-1)
        r(-1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(-1, -l) = r(-1, -l) / scal_uvw_m(l)
        !dst(ind-l) = dst(ind-l) + src(ind-1)*r(-1, -l)
        ! m = -l, n = 0
        u = r1(0, 1)*r_prev(0, 1-l) + r1(0, -1)*r_prev(0, l-1)
        v = r1(1, 1)*r_prev(1, 1-l) + r1(1, -1)*r_prev(1, l-1) + &
            & r1(-1, 1)*r_prev(-1, 1-l) + r1(-1, -1)*r_prev(-1, l-1)
        r(0, -l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
        !dst(ind-l) = dst(ind-l) + src(ind)*r(0, -l)
        ! m = -l, n = 2-l..l-2, n != -1,0,1
        do n = 2, l-2
            u = r1(0, 1)*r_prev(n, 1-l) + r1(0, -1)*r_prev(n, l-1)
            v = r1(1, 1)*r_prev(n-1, 1-l) + r1(1, -1)*r_prev(n-1, l-1) - &
                & r1(-1, 1)*r_prev(1-n, 1-l) - r1(-1, -1)*r_prev(1-n, l-1)
            w = r1(1, 1)*r_prev(n+1, 1-l) + r1(1, -1)*r_prev(n+1, l-1) + &
                & r1(-1, 1)*r_prev(-n-1, 1-l) + r1(-1, -1)*r_prev(-n-1, l-1)
            r(n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(n, -l) = r(n, -l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind+n)*r(n, -l)
            u = r1(0, 1)*r_prev(-n, 1-l) + r1(0, -1)*r_prev(-n, l-1)
            v = r1(1, 1)*r_prev(1-n, 1-l) + r1(1, -1)*r_prev(1-n, l-1) + &
                & r1(-1, 1)*r_prev(n-1, 1-l) + r1(-1, -1)*r_prev(n-1, l-1)
            w = r1(1, 1)*r_prev(-n-1, 1-l) + r1(1, -1)*r_prev(-n-1, l-1) - &
                & r1(-1, 1)*r_prev(n+1, 1-l) - r1(-1, -1)*r_prev(n+1, l-1)
            r(-n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(-n, -l) = r(-n, -l) / scal_uvw_m(l)
            !dst(ind-l) = dst(ind-l) + src(ind-n)*r(-n, -l)
        end do
        ! Now deal with m=1-l..l-1
        do m = 1-l, l-1
            ! n = l
            v = r1(1, 0)*r_prev(l-1, m) - r1(-1, 0)*r_prev(1-l, m)
            r(l, m) = v * scal_v_n(l) / scal_uvw_m(m)
            !dst(ind+m) = src(ind+l) * r(l, m)
            ! n = -l
            v = r1(1, 0)*r_prev(1-l, m) + r1(-1, 0)*r_prev(l-1, m)
            r(-l, m) = v * scal_v_n(l) / scal_uvw_m(m)
            !dst(ind+m) = dst(ind+m) + src(ind-l)*r(-l, m)
            ! n = l-1
            u = r1(0, 0) * r_prev(l-1, m)
            v = r1(1, 0)*r_prev(l-2, m) - r1(-1, 0)*r_prev(2-l, m)
            r(l-1, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
            r(l-1, m) = r(l-1, m) / scal_uvw_m(m)
            !dst(ind+m) = dst(ind+m) + src(ind+l-1)*r(l-1, m)
            ! n = 1-l
            u = r1(0, 0) * r_prev(1-l, m)
            v = r1(1, 0)*r_prev(2-l, m) + r1(-1, 0)*r_prev(l-2, m)
            r(1-l, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
            r(1-l, m) = r(1-l, m) / scal_uvw_m(m)
            !dst(ind+m) = dst(ind+m) + src(ind+1-l)*r(1-l, m)
            ! n = 0
            u = r1(0, 0) * r_prev(0, m)
            v = r1(1, 0)*r_prev(1, m) + r1(-1, 0)*r_prev(-1, m)
            r(0, m) = u*scal_u_n(0) + v*scal_v_n(0)
            r(0, m) = r(0, m) / scal_uvw_m(m)
            !dst(ind+m) = dst(ind+m) + src(ind)*r(0, m)
            ! n = 1
            u = r1(0, 0) * r_prev(1, m)
            v = r1(1, 0) * r_prev(0, m)
            w = r1(1, 0)*r_prev(2, m) + r1(-1, 0)*r_prev(-2, m)
            r(1, m) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(1, m) = r(1, m) / scal_uvw_m(m)
            !dst(ind+m) = dst(ind+m) + src(ind+1)*r(1, m)
            ! n = -1
            u = r1(0, 0) * r_prev(-1, m)
            v = r1(-1, 0) * r_prev(0, m)
            w = r1(1, 0)*r_prev(-2, m) - r1(-1, 0)*r_prev(2, m)
            r(-1, m) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(-1, m) = r(-1, m) / scal_uvw_m(m)
            !dst(ind+m) = dst(ind+m) + src(ind-1)*r(-1, m)
            ! n = 2-l..l-2, n != -1,0,1
            do n = 2, l-2
                u = r1(0, 0) * r_prev(n, m)
                v = r1(1, 0)*r_prev(n-1, m) - r1(-1, 0)*r_prev(1-n, m)
                w = r1(1, 0)*r_prev(n+1, m) + r1(-1, 0)*r_prev(-1-n, m)
                r(n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(n, m) = r(n, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind+n)*r(n, m)
                u = r1(0, 0) * r_prev(-n, m)
                v = r1(1, 0)*r_prev(1-n, m) + r1(-1, 0)*r_prev(n-1, m)
                w = r1(1, 0)*r_prev(-n-1, m) - r1(-1, 0)*r_prev(n+1, m)
                r(-n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(-n, m) = r(-n, m) / scal_uvw_m(m)
                !dst(ind+m) = dst(ind+m) + src(ind-n)*r(-n, m)
            end do
        end do
        do m = -l, l
            ind = 2*l*(l+1)*(2*l+1)/3 + l + 1
            do n = -l, l
                mat(ind+m*(2*l+1)+n) = r(n, m)
            end do
        end do
    end do
end subroutine fmm_sph_transform_get_mat

!> Apply matrix of transformation of spherical harmonics
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha \mathrm{mat} \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! \mathrm{mat} \f$ is a matrix of transformation of spherical harmonics
!! precomputed by @ref fmm_sph_transform_get_mat.
!!
!!
!! @param[in] p: Maximum degree of spherical harmonics
!! @param[in] mat: Transformation matrix in a sparse storage
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multiplier for `dst`
!! @param[inout] dst: Coefficients of transformed spherical harmonics
!!
!! @sa fmm_sph_transform, fmm_sph_transform_get_mat
subroutine fmm_sph_transform_use_mat(p, mat, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: mat((p+1)*(2*p+1)*(2*p+3)/3), alpha, &
        & src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)*(p+1))
    ! Local variables
    integer :: l, m, n, ind, indl
    ! Init dst
    if (beta .eq. zero) then
        dst = zero
    else
        dst = beta * dst
    end if
    ! l = 0
    dst(1) = dst(1) + alpha*src(1)
    if(p .eq. 0) then
        return
    end if
    do l = 1, p
        ! magical value for the offset to the current reflection matrix
        ind = 2*l*(l+1)*(2*l+1)/3 + l + 1
        ! offset for current spherical harmonics
        indl = l*l + l + 1
        do m = -l, l
            do n = -l, l
                dst(indl+m) = dst(indl+m) + &
                    & alpha*mat(ind+m*(2*l+1)+n)*src(indl+n)
            end do
        end do
    end do
end subroutine fmm_sph_transform_use_mat

!> Rotate spherical harmonics around OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! R \f$ is a matrix of rotation of coordinates around OZ axis on angle \f$
!! \phi \f$, presented by \f$ \cos(m \phi) \f$ and \f$ \sin(m \phi) \f$.
!!
!!
!! @param[in] p: Maximal order of spherical harmonics
!! @param[in] vcos: Vector \f$ \{ \cos(m \phi) \}_{m=0}^p \f$
!! @param[in] vsin: Vector \f$ \{ \sin(m \phi) \}_{m=0}^p \f$
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[out] dst: Coefficients of rotated spherical harmonics
subroutine fmm_sph_rotate_oz(p, vcos, vsin, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: vcos(p+1), vsin(p+1), alpha, src((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(out) :: dst((p+1)*(p+1))
    ! Local variables
    integer :: l, m, ind
    ! Init output
    if (beta .eq. zero) then
        dst = zero
    else
        dst = beta * dst
    end if
    ! l = 0
    dst(1) = dst(1) + alpha*src(1)
    ! l > 0
    do l = 1, p
        ind = l*l + l + 1
        ! m = 0
        dst(ind) = dst(ind) + alpha*src(ind)
        ! m != 0
        do m = 1, l
            ! m > 0
            dst(ind+m) = dst(ind+m) + &
                & alpha*(src(ind+m)*vcos(1+m) - src(ind-m)*vsin(1+m))
            ! m < 0
            dst(ind-m) = dst(ind-m) + &
                & alpha*(src(ind+m)*vsin(1+m) + src(ind-m)*vcos(1+m))
        end do
    end do
end subroutine fmm_sph_rotate_oz

!> Transform spherical harmonics in the OXZ plane
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha R \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics corresponding to a new cartesion system of coordinates, \f$
!! \mathrm{src} \f$ is a vector of coefficients of input spherical harmonics
!! corresponding to the standard cartesian system of coordinates and \f$
!! R \f$ is a matrix of transformation of coordinates in the OXZ plane (y
!! coordinate remains the same) presented by 2-by-2 matrix `r1xz`.
!!
!! Based on @ref fmm_sph_transform
!! by assuming `r1(-1, 0) = r1(-1, 1) = r1(0, -1) = r1(1, -1) = 0` and
!! `r1(-1, -1) = 1`, which corresponds to the following transformation matrix:
!! \f[
!!      R_1 = \begin{bmatrix} 1 & 0 & 0 \\ 0 & a & b \\ 0 & c & d
!!      \end{bmatrix},
!! \f]
!! where unkown elements represent input `r1xz` 2x2 array:
!! \f[
!!      R_1^{xz} = \begin{bmatrix} a & b \\ c & d \end{bmatrix}.
!! \f]
!!
!!
!! @param[in] p: Maximal order of spherical harmonics
!! @param[in] r1xz: Transformation from new to old coordinates in the OXZ plane
!! @param[in] alpha: Scalar multiplier for `src`
!! @param[in] src: Coefficients of initial spherical harmonics
!! @param[in] beta: Scalar multipler for `dst`
!! @param[out] dst: Coefficients of transformed spherical harmonics
subroutine fmm_sph_transform_oxz(p, r1xz, alpha, src, beta, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: r1xz(2, 2), alpha, src((p+1)**2), beta
    ! Output
    real(dp), intent(inout) :: dst((p+1)**2)
    ! Local variables
    real(dp) :: tmp((p+1)**2)
    ! Call corresponding routine with built-in alpha=one and beta=zero
    if (beta .eq. zero) then
        call fmm_sph_transform_oxz_out(p, r1xz, alpha*src, dst)
    else
        call fmm_sph_transform_oxz_out(p, r1xz, src, tmp)
        dst = beta*dst + alpha*tmp
    end if
end subroutine fmm_sph_transform_oxz

!> Transform spherical harmonics in the OXZ plane
!!
!! This function implements @ref fmm_sph_transform_oxz with predefined values
!! of parameters \p alpha=one and \p beta=zero.
!! 
!!
!! @param[in] p: maximum order of spherical harmonics
!! @param[in] r1xz: 2D transformation matrix in the OXZ plane
!! @param[in] src: coefficients of initial spherical harmonics
!! @param[out] dst: coefficients of rotated spherical harmonics
subroutine fmm_sph_transform_oxz_out(p, r1xz, src, dst)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: r1xz(0:1, 0:1), src((p+1)**2)
    ! Output
    real(dp), intent(out) :: dst((p+1)*(p+1))
    ! Local variables
    real(dp) :: u, v, w, r_prev(-p:p, -p:p), r(-p:p, -p:p)
    real(dp) :: scal_uvw_m(-p:p), scal_u_n(-p:p), scal_v_n(-p:p)
    real(dp) :: scal_w_n(-p:p)
    integer :: l, m, n, ind
    ! Compute rotations/reflections
    ! l = 0
    dst(1) = src(1)
    if (p .eq. 0) then
        return
    end if
    ! l = 1
    dst(2) = src(2)
    dst(3) = src(3)*r1xz(0, 0) + src(4)*r1xz(1, 0)
    dst(4) = src(3)*r1xz(0, 1) + src(4)*r1xz(1, 1)
    if (p .eq. 1) then
        return
    end if
    ! l = 2
    r(2, 2) = (r1xz(1, 1)*r1xz(1, 1) + 1) / 2
    r(2, 1) = r1xz(1, 1)*r1xz(1, 0)
    r(2, 0) = sqrt(dble(3)) / 2 * r1xz(1, 0) * r1xz(1, 0)
    r(2, -1) = 0
    r(2, -2) = 0
    r(1, 2) = r1xz(1, 1)*r1xz(0, 1)
    r(1, 1) = r1xz(1, 1)*r1xz(0, 0) + r1xz(1, 0)*r1xz(0, 1)
    r(1, 0) = sqrt(dble(3)) * r1xz(1, 0) * r1xz(0, 0)
    r(1, -1) = 0
    r(1, -2) = 0
    r(0, 2) = sqrt(dble(3)) / 2 * r1xz(0, 1) * r1xz(0, 1)
    r(0, 1) = sqrt(dble(3)) * r1xz(0, 1) * r1xz(0, 0)
    r(0, 0) = (3*r1xz(0, 0)*r1xz(0, 0)-1) / 2
    r(0, -1) = 0
    r(0, -2) = 0
    r(-1, 2) = 0
    r(-1, 1) = 0
    r(-1, 0) = 0
    r(-1, -1) = r1xz(0, 0)
    r(-1, -2) = r1xz(0, 1)
    r(-2, 2) = 0
    r(-2, 1) = 0
    r(-2, 0) = 0
    r(-2, -1) = r1xz(1, 0)
    r(-2, -2) = r1xz(1, 1)
    dst(9) = src(9)*r(2, 2) + src(8)*r(1, 2) + src(7)*r(0, 2)
    dst(8) = src(9)*r(2, 1) + src(8)*r(1, 1) + src(7)*r(0, 1)
    dst(7) = src(9)*r(2, 0) + src(8)*r(1, 0) + src(7)*r(0, 0)
    dst(6) = src(6)*r(-1, -1) + src(5)*r(-2, -1)
    dst(5) = src(6)*r(-1, -2) + src(5)*r(-2, -2)
    ! l > 2
    do l = 3, p
        ! Prepare previous rotation matrix
        r_prev(1-l:l-1, 1-l:l-1) = r(1-l:l-1, 1-l:l-1)
        ! Prepare scalar factors
        scal_uvw_m(0) = l
        do m = 1, l-1
        u = sqrt(dble(l*l-m*m))
        scal_uvw_m(m) = u
        scal_uvw_m(-m) = u
        end do
        u = 2 * l
        u = sqrt(dble(u*(u-1)))
        scal_uvw_m(l) = u
        scal_uvw_m(-l) = u
        scal_u_n(0) = l
        scal_v_n(0) = -sqrt(dble(l*(l-1))) / sqrt2
        scal_w_n(0) = 0
        do n = 1, l-2
        u = sqrt(dble(l*l-n*n))
        scal_u_n(n) = u
        scal_u_n(-n) = u
        v = l + n
        v = sqrt(v*(v-1)) / 2
        scal_v_n(n) = v
        scal_v_n(-n) = v
        w = l - n
        w = -sqrt(w*(w-1)) / 2
        scal_w_n(n) = w
        scal_w_n(-n) = w
        end do
        u = sqrt(dble(2*l-1))
        scal_u_n(l-1) = u
        scal_u_n(1-l) = u
        scal_u_n(l) = 0
        scal_u_n(-l) = 0
        v = sqrt(dble((2*l-1)*(2*l-2))) / 2
        scal_v_n(l-1) = v
        scal_v_n(1-l) = v
        v = sqrt(dble(2*l*(2*l-1))) / 2
        scal_v_n(l) = v
        scal_v_n(-l) = v
        scal_w_n(l-1) = 0
        scal_w_n(l) = 0
        scal_w_n(-l) = 0
        scal_w_n(1-l) = 0
        ind = l*l + l + 1
        ! m = l, n = l
        v = r1xz(1, 1)*r_prev(l-1, l-1) + r_prev(1-l, 1-l)
        r(l, l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind+l) = src(ind+l) * r(l, l)
        ! m = l, n = -l
        v = r1xz(1, 1)*r_prev(1-l, l-1) - r_prev(l-1, 1-l)
        r(-l, l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind-l)*r(-l, l)
        ! m = l, n = l-1
        u = r1xz(0, 1) * r_prev(l-1, l-1)
        v = r1xz(1, 1)*r_prev(l-2, l-1) + r_prev(2-l, 1-l)
        r(l-1, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind+l-1)*r(l-1, l)
        ! m = l, n = 1-l
        u = r1xz(0, 1) * r_prev(1-l, l-1)
        v = r1xz(1, 1)*r_prev(2-l, l-1) - r_prev(l-2, 1-l)
        r(1-l, l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind+1-l)*r(1-l, l)
        ! m = l, n = 1
        u = r1xz(0, 1) * r_prev(1, l-1)
        v = r1xz(1, 1) * r_prev(0, l-1)
        w = r1xz(1, 1)*r_prev(2, l-1) - r_prev(-2, 1-l)
        r(1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(1, l) = r(1, l) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind+1)*r(1, l)
        ! m = l, n = -1
        u = r1xz(0, 1) * r_prev(-1, l-1)
        v = r_prev(0, 1-l)
        w = r1xz(1, 1)*r_prev(-2, l-1) + r_prev(2, 1-l)
        r(-1, l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(-1, l) = r(-1, l) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind-1)*r(-1, l)
        ! m = l, n = 0
        u = r1xz(0, 1) * r_prev(0, l-1)
        v = r1xz(1, 1)*r_prev(1, l-1) - r_prev(-1, 1-l)
        r(0, l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
        dst(ind+l) = dst(ind+l) + src(ind)*r(0, l)
        ! m = l, n = 2-l..l-2, n != -1,0,1
        do n = 2, l-2
            u = r1xz(0, 1) * r_prev(n, l-1)
            v = r1xz(1, 1)*r_prev(n-1, l-1) + r_prev(1-n, 1-l)
            w = r1xz(1, 1)*r_prev(n+1, l-1) - r_prev(-n-1, 1-l)
            r(n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(n, l) = r(n, l) / scal_uvw_m(l)
            dst(ind+l) = dst(ind+l) + src(ind+n)*r(n, l)
            u = r1xz(0, 1) * r_prev(-n, l-1)
            v = r1xz(1, 1)*r_prev(1-n, l-1) - r_prev(n-1, 1-l)
            w = r1xz(1, 1)*r_prev(-n-1, l-1) + r_prev(n+1, 1-l)
            r(-n, l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(-n, l) = r(-n, l) / scal_uvw_m(l)
            dst(ind+l) = dst(ind+l) + src(ind-n)*r(-n, l)
        end do
        ! m = -l, n = l
        v = r1xz(1, 1)*r_prev(l-1, 1-l) - r_prev(1-l, l-1)
        r(l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind-l) = src(ind+l)*r(l, -l)
        ! m = -l, n = -l
        v = r1xz(1, 1)*r_prev(1-l, 1-l) + r_prev(l-1, l-1)
        r(-l, -l) = v * scal_v_n(l) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind-l)*r(-l, -l)
        ! m = -l, n = l-1
        u = r1xz(0, 1) * r_prev(l-1, 1-l)
        v = r1xz(1, 1)*r_prev(l-2, 1-l) - r_prev(2-l, l-1)
        r(l-1, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind+l-1)*r(l-1, -l)
        ! m = -l, n = 1-l
        u = r1xz(0, 1) * r_prev(1-l, 1-l)
        v = r1xz(1, 1)*r_prev(2-l, 1-l) + r_prev(l-2, l-1)
        r(1-l, -l) = (u*scal_u_n(l-1)+v*scal_v_n(l-1)) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind+1-l)*r(1-l, -l)
        ! m = -l, n = 1
        u = r1xz(0, 1) * r_prev(1, 1-l)
        v = r1xz(1, 1) * r_prev(0, 1-l)
        w = r1xz(1, 1)*r_prev(2, 1-l) + r_prev(-2, l-1)
        r(1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(1, -l) = r(1, -l) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind+1)*r(1, -l)
        ! m = -l, n = -1
        u = r1xz(0, 1) * r_prev(-1, 1-l)
        v = r_prev(0, l-1)
        w = r1xz(1, 1)*r_prev(-2, 1-l) - r_prev(2, l-1)
        r(-1, -l) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
        r(-1, -l) = r(-1, -l) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind-1)*r(-1, -l)
        ! m = -l, n = 0
        u = r1xz(0, 1) * r_prev(0, 1-l)
        v = r1xz(1, 1)*r_prev(1, 1-l) + r_prev(-1, l-1)
        r(0, -l) = (u*scal_u_n(0) + v*scal_v_n(0)) / scal_uvw_m(l)
        dst(ind-l) = dst(ind-l) + src(ind)*r(0, -l)
        ! m = -l, n = 2-l..l-2, n != -1,0,1
        do n = 2, l-2
            u = r1xz(0, 1) * r_prev(n, 1-l)
            v = r1xz(1, 1)*r_prev(n-1, 1-l) - r_prev(1-n, l-1)
            w = r1xz(1, 1)*r_prev(n+1, 1-l) + r_prev(-n-1, l-1)
            r(n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(n, -l) = r(n, -l) / scal_uvw_m(l)
            dst(ind-l) = dst(ind-l) + src(ind+n)*r(n, -l)
            u = r1xz(0, 1) * r_prev(-n, 1-l)
            v = r1xz(1, 1)*r_prev(1-n, 1-l) + r_prev(n-1, l-1)
            w = r1xz(1, 1)*r_prev(-n-1, 1-l) - r_prev(n+1, l-1)
            r(-n, -l) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
            r(-n, -l) = r(-n, -l) / scal_uvw_m(l)
            dst(ind-l) = dst(ind-l) + src(ind-n)*r(-n, -l)
        end do
        ! Now deal with m=1-l..l-1
        do m = 1-l, l-1
            ! n = l
            v = r1xz(1, 0) * r_prev(l-1, m)
            r(l, m) = v * scal_v_n(l) / scal_uvw_m(m)
            dst(ind+m) = src(ind+l) * r(l, m)
            ! n = -l
            v = r1xz(1, 0) * r_prev(1-l, m)
            r(-l, m) = v * scal_v_n(l) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind-l)*r(-l, m)
            ! n = l-1
            u = r1xz(0, 0) * r_prev(l-1, m)
            v = r1xz(1, 0) * r_prev(l-2, m)
            r(l-1, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
            r(l-1, m) = r(l-1, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind+l-1)*r(l-1, m)
            ! n = 1-l
            u = r1xz(0, 0) * r_prev(1-l, m)
            v = r1xz(1, 0) * r_prev(2-l, m)
            r(1-l, m) = u*scal_u_n(l-1) + v*scal_v_n(l-1)
            r(1-l, m) = r(1-l, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind+1-l)*r(1-l, m)
            ! n = 0
            u = r1xz(0, 0) * r_prev(0, m)
            v = r1xz(1, 0) * r_prev(1, m)
            r(0, m) = u*scal_u_n(0) + v*scal_v_n(0)
            r(0, m) = r(0, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind)*r(0, m)
            ! n = 1
            u = r1xz(0, 0) * r_prev(1, m)
            v = r1xz(1, 0) * r_prev(0, m)
            w = r1xz(1, 0) * r_prev(2, m)
            r(1, m) = u*scal_u_n(1) + sqrt2*v*scal_v_n(1) + w*scal_w_n(1)
            r(1, m) = r(1, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind+1)*r(1, m)
            ! n = -1
            u = r1xz(0, 0) * r_prev(-1, m)
            w = r1xz(1, 0) * r_prev(-2, m)
            r(-1, m) = u*scal_u_n(1) + w*scal_w_n(1)
            r(-1, m) = r(-1, m) / scal_uvw_m(m)
            dst(ind+m) = dst(ind+m) + src(ind-1)*r(-1, m)
            ! n = 2-l..l-2, n != -1,0,1
            do n = 2, l-2
                u = r1xz(0, 0) * r_prev(n, m)
                v = r1xz(1, 0) * r_prev(n-1, m)
                w = r1xz(1, 0) * r_prev(n+1, m)
                r(n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(n, m) = r(n, m) / scal_uvw_m(m)
                dst(ind+m) = dst(ind+m) + src(ind+n)*r(n, m)
                u = r1xz(0, 0) * r_prev(-n, m)
                v = r1xz(1, 0) * r_prev(1-n, m)
                w = r1xz(1, 0) * r_prev(-n-1, m)
                r(-n, m) = u*scal_u_n(n) + v*scal_v_n(n) + w*scal_w_n(n)
                r(-n, m) = r(-n, m) / scal_uvw_m(m)
                dst(ind+m) = dst(ind+m) + src(ind-n)*r(-n, m)
            end do
        end do
    end do
end subroutine fmm_sph_transform_oxz_out

!> Direct M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multipler for `alpha`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_m2m_ztranslate(z, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, indj, indn, indjn
    ! Init output
    if (beta .eq. zero) then
        dst_m = zero
    else
        dst_m = beta * dst_m
    end if
    ! If harmonics have different centers
    if (z .ne. 0) then
        r1 = src_r / dst_r
        r2 = z / dst_r
        pow_r1(1) = r1
        pow_r2(1) = one
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = alpha * vscales(indj) * vfact(j-k+1) * vfact(j+k+1)
                do n = 0, j-k
                    indn = n*n + n + 1
                    indjn = (j-n)**2 + (j-n) + 1
                    tmp2 = tmp1 * pow_r1(j-n+1) * pow_r2(n+1) / &
                        & vscales(indjn) / vfact(n+1) / vfact(n+1) / &
                        & vfact(j-n-k+1) / vfact(j-n+k+1)
                    if (k .eq. 0) then
                        dst_m(indj) = dst_m(indj) + tmp2*src_m(indjn)
                    else
                        dst_m(indj+k) = dst_m(indj+k) + tmp2*src_m(indjn+k)
                        dst_m(indj-k) = dst_m(indj-k) + tmp2*src_m(indjn-k)
                    end if
                end do
            end do
        end do
    ! If harmonics are located at the same point
    else
        r1 = src_r / dst_r
        tmp1 = alpha * r1
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_m(k) = dst_m(k) + src_m(k)*tmp1
            end do
            tmp1 = tmp1 * r1
        end do
    end if
end subroutine fmm_m2m_ztranslate

!> Adjoint M2M translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M^\top \f$ is an adjoint matrix of
!! multipole-to-multipole translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multipler for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_m2m_ztranslate_adj(z, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, indj, indn, indjn
    ! Init output
    if (beta .eq. zero) then
        dst_m = zero
    else
        dst_m = beta * dst_m
    end if
    ! If harmonics have different centers
    if (z .ne. 0) then
        r1 = dst_r / src_r
        r2 = -z / src_r
        pow_r1(1) = r1
        pow_r2(1) = one
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = alpha * vscales(indj) * vfact(j-k+1) * vfact(j+k+1)
                do n = 0, j-k
                    indn = n*n + n + 1
                    indjn = (j-n)**2 + (j-n) + 1
                    tmp2 = tmp1 * pow_r1(j-n+1) * pow_r2(n+1) / &
                        & vscales(indjn) / vfact(n+1) / vfact(n+1) / &
                        & vfact(j-n-k+1) / vfact(j-n+k+1)
                    if (k .eq. 0) then
                        dst_m(indjn) = dst_m(indjn) + tmp2*src_m(indj)
                    else
                        dst_m(indjn+k) = dst_m(indjn+k) + tmp2*src_m(indj+k)
                        dst_m(indjn-k) = dst_m(indjn-k) + tmp2*src_m(indj-k)
                    end if
                end do
            end do
        end do
    ! If harmonics are located at the same point
    else
        r1 = dst_r / src_r
        tmp1 = alpha * r1
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_m(k) = dst_m(k) + src_m(k)*tmp1
            end do
            tmp1 = tmp1 * r1
        end do
    end if
end subroutine fmm_m2m_ztranslate_adj

!> Save matrix of M2M translation along OZ axis
!!
!! In a case input `z` is zero no translation matrix is computed, as @ref
!! fmm_m2m_scale shall be used in this case without any precomputed matrices
!!
!!
!! @param[in] z: the OZ coordinate of the radius-vector from new to old centers
!!      of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$
!! @param[in] vfact: Square roots of factorials
!! @param[in] mat: Translation matrix for spherical harmonics
subroutine fmm_m2m_ztranslate_get_mat(z, src_r, dst_r, p, vscales, vfact, mat)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), vfact(2*p+1)
    integer, intent(in) :: p
    ! Output
    real(dp), intent(out) :: mat((p+1)*(p+2)*(p+3)/6)
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, indj, indn, indjn, indmat
    if (z .eq. zero) then
        ! Fill matrix with NaNs
        tmp1 = zero
        mat = tmp1 / tmp1
        return
    end if
    r1 = src_r / dst_r
    r2 = z / dst_r
    pow_r1(1) = r1
    pow_r2(1) = one
    do j = 2, p+1
        pow_r1(j) = pow_r1(j-1) * r1
        pow_r2(j) = pow_r2(j-1) * r2
    end do
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        do k = 0, j
            tmp1 = vscales(indj) * vfact(j-k+1) * vfact(j+k+1)
            do n = 0, j-k
                indn = n*n + n + 1
                indjn = (j-n)**2 + (j-n) + 1
                mat(indmat) = tmp1 * pow_r1(j-n+1) * pow_r2(n+1) / &
                    & vscales(indjn) / vfact(n+1) / vfact(n+1) / &
                    & vfact(j-n-k+1) / vfact(j-n+k+1)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2m_ztranslate_get_mat

!> Apply matrix of M2M translation along OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is an matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!! This function shall not be used in case if centers of spherical harmonics
!! are identical. The OZ translation matrix is not computed by @ref
!! fmm_m2m_ztranslate_get_mat in this case as @ref fmm_m2m_scale is intended to
!! treat the case.
!!
!!
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] mat: The OZ translation matrix for spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion in new harmonics
subroutine fmm_m2m_ztranslate_use_mat(p, mat, alpha, src_m, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: mat((p+1)*(p+2)*(p+3)/6), alpha, &
        & src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    integer :: indmat, j, indj, n, indn, indjn, k
    real(dp) :: tmp1, tmp2
    ! Cycle over matrix elements stored in a sparse way
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        ! Init output properly
        if (beta .eq. zero) then
            dst_m(indj-j:indj+j) = zero
        else
            dst_m(indj-j:indj+j) = beta * dst_m(indj-j:indj+j)
        end if
        ! k = 0
        tmp1 = zero
        do n = 0, j
            indn = n*n + n + 1
            indjn = (j-n)**2 + (j-n) + 1
            tmp1 = tmp1 + mat(indmat)*src_m(indjn)
            indmat = indmat + 1
        end do
        dst_m(indj) = dst_m(indj) + alpha*tmp1
        ! k > 0
        do k = 1, j
            tmp1 = zero
            tmp2 = zero
            do n = 0, j-k
                indn = n*n + n + 1
                indjn = (j-n)**2 + (j-n) + 1
                tmp1 = tmp1 + mat(indmat)*src_m(indjn+k)
                tmp2 = tmp2 + mat(indmat)*src_m(indjn-k)
                indmat = indmat + 1
            end do
            dst_m(indj+k) = dst_m(indj+k) + alpha*tmp1
            dst_m(indj-k) = dst_m(indj-k) + alpha*tmp2
        end do
    end do
end subroutine fmm_m2m_ztranslate_use_mat

!> Apply adjoint matrix of M2M translation along OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M^\top \f$ is an adjoint matrix of a
!! multipole-to-multipole translation over OZ axis.
!!
!! This function shall not be used in case if centers of spherical harmonics
!! are identical. The OZ translation matrix is not computed by @ref
!! fmm_m2m_ztranslate_get_mat in this case as @ref fmm_m2m_scale is intended to
!! treat the case.
!!
!!
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] mat: The OZ translation matrix for spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion in new harmonics
!!
!!
!! TODO: improve performance
subroutine fmm_m2m_ztranslate_use_mat_adj(p, mat, alpha, src_m, beta, dst_m)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: mat((p+1)*(p+2)*(p+3)/6), alpha, &
        & src_m((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    integer :: indmat, j, indj, n, indn, indjn, k
    real(dp) :: tmp1, tmp2
    ! Init output properly
    if (beta .eq. zero) then
        dst_m = zero
    else
        dst_m = beta * dst_m
    end if
    ! Cycle over matrix elements stored in a sparse way
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        ! k = 0
        do n = 0, j
            indn = n*n + n + 1
            indjn = (j-n)**2 + (j-n) + 1
            dst_m(indjn) = dst_m(indjn) + alpha*mat(indmat)*src_m(indj)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            do n = 0, j-k
                indn = n*n + n + 1
                indjn = (j-n)**2 + (j-n) + 1
                dst_m(indjn+k) = dst_m(indjn+k) + &
                    & alpha*mat(indmat)*src_m(indj+k)
                dst_m(indjn-k) = dst_m(indjn-k) + &
                    & alpha*mat(indmat)*src_m(indj-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2m_ztranslate_use_mat_adj

!> Scale M2M, when spherical harmonics are centered in the same point
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of a multipole-to-multipole
!! translation when harmonics are centered in the same point.
!!
!!
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion in new harmonics
subroutine fmm_m2m_scale(src_r, dst_r, p, alpha, src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: src_r, dst_r, alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, tmp1
    integer :: j, k, indj
    ! Init output
    if (beta .eq. zero) then
        dst_m = zero
    else
        dst_m = beta * dst_m
    end if
    r1 = src_r / dst_r
    tmp1 = alpha * r1
    do j = 0, p
        indj = j*j + j + 1
        do k = indj-j, indj+j
            dst_m(k) = dst_m(k) + src_m(k)*tmp1
        end do
        tmp1 = tmp1 * r1
    end do
end subroutine fmm_m2m_scale

!> Adjoint scale M2M, when spherical harmonics are centered in the same point
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M^\top \f$ is an adjoint matrix of a
!! multipole-to-multipole translation when harmonics are centered in the same
!! point.
!!
!!
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion in new harmonics
subroutine fmm_m2m_scale_adj(src_r, dst_r, p, alpha, src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: src_r, dst_r, alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, tmp1
    integer :: j, k, indj
    ! Init output
    if (beta .eq. zero) then
        dst_m = zero
    else
        dst_m = beta * dst_m
    end if
    r1 = dst_r / src_r
    tmp1 = alpha * r1
    do j = 0, p
        indj = j*j + j + 1
        do k = indj-j, indj+j
            dst_m(k) = dst_m(k) + src_m(k)*tmp1
        end do
        tmp1 = tmp1 * r1
    end do
end subroutine fmm_m2m_scale_adj

!> Direct M2M translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M \f$ is a matrix of a multipole-to-multipole
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_m2m_rotation(c, src_r, dst_r, p, vscales, vfact, alpha, src_m, &
        & beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1), &
        & vmsin(p+1), tmp_m((p+1)*(p+1)), r1xz(2, 2), tmp_m2((p+1)*(p+1))
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. zero) then
        call fmm_m2m_ztranslate(c(3), src_r, dst_r, p, vscales, vfact, alpha, &
            & src_m, beta, dst_m)
        return
    end if
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis
    vmsin = -vsin
    call fmm_sph_rotate_oz(p, vcos, vmsin, alpha, src_m, zero, tmp_m)
    ! Perform rotation in the OXZ plane
    r1xz(1, 1) = ctheta
    r1xz(1, 2) = -stheta
    r1xz(2, 1) = stheta
    r1xz(2, 2) = ctheta
    call fmm_sph_transform_oxz(p, r1xz, one, tmp_m, zero, tmp_m2)
    ! OZ translation
    call fmm_m2m_ztranslate(rho, src_r, dst_r, p, vscales, vfact, one, &
        & tmp_m2, zero, tmp_m)
    ! Backward rotation in the OXZ plane
    r1xz(1, 2) = stheta
    r1xz(2, 1) = -stheta
    call fmm_sph_transform_oxz(p, r1xz, one, tmp_m, zero, tmp_m2)
    ! Backward rotation around OZ axis
    call fmm_sph_rotate_oz(p, vcos, vsin, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_rotation

!> Adjoint M2M translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha M_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ M_M^\top \f$ is aa adjoint matrix of a
!! multipole-to-multipole translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] src_m: expansion in old harmonics
!! @param[inout] dst_m: expansion in new harmonics
subroutine fmm_m2m_rotation_adj(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1), &
        & vmsin(p+1), tmp_m((p+1)*(p+1)), r1xz(2, 2), tmp_m2((p+1)*(p+1))
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. 0) then
        call fmm_m2m_ztranslate_adj(c(3), src_r, dst_r, p, vscales, vfact, &
            & alpha, src_m, beta, dst_m)
        return
    end if
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis
    vmsin = -vsin
    call fmm_sph_rotate_oz(p, vcos, vmsin, alpha, src_m, zero, tmp_m)
    ! Perform rotation in the OXZ plane
    r1xz(1, 1) = ctheta
    r1xz(1, 2) = -stheta
    r1xz(2, 1) = stheta
    r1xz(2, 2) = ctheta
    call fmm_sph_transform_oxz(p, r1xz, one, tmp_m, zero, tmp_m2)
    ! OZ translation
    call fmm_m2m_ztranslate_adj(rho, src_r, dst_r, p, vscales, vfact, one, &
        & tmp_m2, zero, tmp_m)
    ! Backward rotation in the OXZ plane
    r1xz(1, 2) = stheta
    r1xz(2, 1) = -stheta
    call fmm_sph_transform_oxz(p, r1xz, one, tmp_m, zero, tmp_m2)
    ! Backward rotation around OZ axis
    call fmm_sph_rotate_oz(p, vcos, vsin, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_rotation_adj

!> Compute reflection that aligns given vector along OZ axis
!!
!! Computes stable reflection, just as in Householder QR algorithm. This
!! procedure is for testing purposes here, as it is already inplace in
!! performance-oriented implementations.
!!
!! @param[in] c: Input vector to align along OZ axis
!! @param[out] z: OZ coordinate of the reflected vector `c`
!! @param[out] mat: Reflection matrix from new (y,z,x) to old (y,z,x)
subroutine coord_reflect_get_mat(c, z, mat)
    ! Input
    real(dp), intent(in) :: c(3)
    ! Output
    real(dp), intent(out) :: z, mat(3, 3)
    ! Local variables
    real(dp) :: max123, ssq123, c1(3), c1_norm
    integer :: m, n
    ! Compute rho (norm of the input vector c)
    max123 = abs(c(1))
    ssq123 = one
    if (c(2) .ne. zero) then
        if (abs(c(2)) .gt. max123) then
            ssq123 = one + ssq123*(max123/c(2))**2
            max123 = abs(c(2))
        else
            ssq123 = ssq123 + (c(2)/max123)**2
        end if
    end if
    if (c(3) .ne. zero) then
        if (abs(c(3)) .gt. max123) then
            ssq123 = one + ssq123*(max123/c(3))**2
            max123 = abs(c(3))
        else
            ssq123 = ssq123 + (c(3)/max123)**2
        end if
    end if
    z = max123 * sqrt(ssq123)
    ! Reorder (x,y,z) -> (y,z,x) and set proper reflection normal vector c1
    if (c(3) .ge. zero) z = -z
    c1(1) = c(2)
    c1(2) = c(3) - z ! -z and c(3) have the same sign
    c1(3) = c(1)
    ! Normalize vector c1. We know for sure c1(2) is maximum and is non-zero
    max123 = abs(c1(2))
    ssq123 = one + (c1(1)/c1(2))**2 + (c1(3)/c1(2))**2
    c1_norm = max123 * sqrt(ssq123)
    c1 = c1 / c1_norm
    mat = zero
    do m = 1, 3
        mat(m, m) = 1
        do n = 1, 3
            mat(n, m) = mat(n, m) - two*c1(n)*c1(m)
        end do
    end do
end subroutine coord_reflect_get_mat

!> Direct M2M translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: expansion in new harmonics
subroutine fmm_m2m_reflection(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: max123, ssq123, rho, c1(3), c1_norm, r1(3, 3), &
        & tmp_m((p+1)*(p+1)), tmp_m2((p+1)*(p+1))
    integer :: m, n
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        call fmm_m2m_ztranslate(c(3), src_r, dst_r, p, vscales, vfact, alpha, &
            & src_m, beta, dst_m)
        return
    end if
    ! Compute rho, c(1:2) != zero already
    if (abs(c(1)) .gt. abs(c(2))) then
        max123 = abs(c(1))
        ssq123 = one + (c(2)/c(1))**2
    else
        max123 = abs(c(2))
        ssq123 = one + (c(1)/c(2))**2
    end if
    if (abs(c(3)) .gt. max123) then
        rho = one + ssq123*(max123/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq123 + (c(3)/max123)**2
        rho = max123 * sqrt(rho)
    end if
    ! Reorder (x,y,z) -> (y,z,x) and set proper reflection normal vector c1
    if (c(3) .ge. zero) rho = -rho
    c1(1) = c(2)
    c1(2) = c(3) - rho ! -rho and c(3) have the same sign
    c1(3) = c(1)
    ! Normalize vector c1. We know for sure c1(2) is maximum and is non-zero
    max123 = abs(c1(2))
    ssq123 = one + (c1(1)/c1(2))**2 + (c1(3)/c1(2))**2
    c1_norm = max123 * sqrt(ssq123)
    c1 = c1 / c1_norm
    r1 = 0
    do m = 1, 3
        r1(m, m) = 1
        do n = 1, 3
            r1(n, m) = r1(n, m) - two*c1(n)*c1(m)
        end do
    end do
    call fmm_sph_transform(p, r1, alpha, src_m, zero, tmp_m)
    call fmm_m2m_ztranslate(rho, src_r, dst_r, p, vscales, vfact, one, tmp_m, &
        & zero, tmp_m2)
    call fmm_sph_transform(p, r1, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_reflection

!> Adjoint M2M translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] src_m: expansion in old harmonics
!! @param[inout] dst_m: expansion in new harmonics
subroutine fmm_m2m_reflection_adj(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: max123, ssq123, rho, c1(3), c1_norm, r1(3, 3), &
        & tmp_m((p+1)*(p+1)), tmp_m2((p+1)*(p+1))
    integer :: m, n
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        call fmm_m2m_ztranslate_adj(c(3), src_r, dst_r, p, vscales, vfact, &
            & alpha, src_m, beta, dst_m)
        return
    end if
    ! Compute rho, c(1:2) != zero already
    if (abs(c(1)) .gt. abs(c(2))) then
        max123 = abs(c(1))
        ssq123 = one + (c(2)/c(1))**2
    else
        max123 = abs(c(2))
        ssq123 = one + (c(1)/c(2))**2
    end if
    if (abs(c(3)) .gt. max123) then
        rho = one + ssq123*(max123/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq123 + (c(3)/max123)**2
        rho = max123 * sqrt(rho)
    end if
    ! Reorder (x,y,z) -> (y,z,x) and set proper reflection normal vector c1
    if (c(3) .ge. zero) rho = -rho
    c1(1) = c(2)
    c1(2) = c(3) - rho ! -rho and c(3) have the same sign
    c1(3) = c(1)
    ! Normalize vector c1. We know for sure c1(2) is maximum and is non-zero
    max123 = abs(c1(2))
    ssq123 = one + (c1(1)/c1(2))**2 + (c1(3)/c1(2))**2
    c1_norm = max123 * sqrt(ssq123)
    c1 = c1 / c1_norm
    r1 = 0
    do m = 1, 3
        r1(m, m) = 1
        do n = 1, 3
            r1(n, m) = r1(n, m) - two*c1(n)*c1(m)
        end do
    end do
    call fmm_sph_transform(p, r1, alpha, src_m, zero, tmp_m)
    call fmm_m2m_ztranslate_adj(rho, src_r, dst_r, p, vscales, vfact, one, &
        & tmp_m, zero, tmp_m2)
    call fmm_sph_transform(p, r1, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_reflection_adj

!> M2M translation by 2 reflections and 1 translation
!!
!! Slightly shorter implementation of @ref fmm_m2m_reflection
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion in new harmonics
subroutine fmm_m2m_reflection2(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: rho, r1(3, 3), tmp_m((p+1)*(p+1)), tmp_m2((p+1)*(p+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        call fmm_m2m_ztranslate(c(3), src_r, dst_r, p, vscales, vfact, alpha, &
            & src_m, beta, dst_m)
        return
    end if
    call coord_reflect_get_mat(c, rho, r1)
    call fmm_sph_transform(p, r1, alpha, src_m, zero, tmp_m)
    call fmm_m2m_ztranslate(rho, src_r, dst_r, p, vscales, vfact, one, tmp_m, &
        & zero, tmp_m2)
    call fmm_sph_transform(p, r1, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_reflection2

!> Adjoint M2M translation by 2 reflections and 1 translation
!!
!! Slightly shorter implementation of @ref fmm_m2m_reflection
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion in new harmonics
subroutine fmm_m2m_reflection2_adj(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: rho, r1(3, 3), tmp_m((p+1)*(p+1)), tmp_m2((p+1)*(p+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        call fmm_m2m_ztranslate_adj(c(3), src_r, dst_r, p, vscales, vfact, &
            & alpha, src_m, beta, dst_m)
        return
    end if
    call coord_reflect_get_mat(c, rho, r1)
    call fmm_sph_transform(p, r1, alpha, src_m, zero, tmp_m)
    call fmm_m2m_ztranslate_adj(rho, src_r, dst_r, p, vscales, vfact, one, &
        & tmp_m, zero, tmp_m2)
    call fmm_sph_transform(p, r1, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_reflection2_adj

!> Save matrices of M2M operation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[out] transform_mat: Matrix of reflection
!! @param[out] ztranslate_mat: Matrix of OZ translation
subroutine fmm_m2m_reflection_get_mat(c, src_r, dst_r, p, vscales, vfact, &
        & transform_mat, ztranslate_mat)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1)
    integer, intent(in) :: p
    ! Outputs
    real(dp), intent(out) :: transform_mat((p+1)*(2*p+1)*(2*p+3)/6)
    real(dp), intent(out) :: ztranslate_mat((p+1)*(p+2)*(p+3)/6)
    ! Local variables
    real(dp) :: rho, r1(3, 3)
    ! If no need for transformation, just generate translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If even z coordinate is zero do nothing, as this case requires a
        ! simple scaling
        if (c(3) .ne. zero) then
            call fmm_m2m_ztranslate_get_mat(c(3), src_r, dst_r, p, vscales, &
                & vfact, ztranslate_mat)
        end if
        return
    end if
    call coord_reflect_get_mat(c, rho, r1)
    call fmm_sph_transform_get_mat(p, r1, transform_mat)
    call fmm_m2m_ztranslate_get_mat(rho, src_r, dst_r, p, vscales, vfact, &
        & ztranslate_mat)
end subroutine fmm_m2m_reflection_get_mat

!> Apply matrices of M2M translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] transform_mat: Matrix of a reflection
!! @param[in] ztranslate_mat: Matrix of an OZ translation
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion of old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion of new harmonics
subroutine fmm_m2m_reflection_use_mat(c, src_r, dst_r, p, transform_mat, &
        & ztranslate_mat, alpha, src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, &
        & transform_mat((p+1)*(2*p+1)*(2*p+3)/3), &
        & ztranslate_mat((p+1)*(p+2)*(p+3)/6), &
        & alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: tmp_m((p+1)*(p+1)), tmp_m2((p+1)*(p+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If centers are the same just scale
        if (c(3) .eq. 0) then
            call fmm_m2m_scale(src_r, dst_r, p, alpha, src_m, beta, dst_m)
        ! Otherwise apply ztranslate matrix
        else
            call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, alpha, src_m, &
                & beta, dst_m)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_transform_use_mat(p, transform_mat, alpha, src_m, zero, tmp_m)
    call fmm_m2m_ztranslate_use_mat(p, ztranslate_mat, one, tmp_m, zero, &
        & tmp_m2)
    call fmm_sph_transform_use_mat(p, transform_mat, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_reflection_use_mat

!> Adjoint apply matrices of M2M translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] transform_mat: Matrix of a reflection
!! @param[in] ztranslate_mat: Matrix of an OZ translation
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion of old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion of new harmonics
subroutine fmm_m2m_reflection_use_mat_adj(c, src_r, dst_r, p, transform_mat, &
        & ztranslate_mat, alpha, src_m, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, &
        & transform_mat((p+1)*(2*p+1)*(2*p+3)/3), &
        & ztranslate_mat((p+1)*(p+2)*(p+3)/6), &
        & alpha, src_m((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_m((p+1)*(p+1))
    ! Local variables
    real(dp) :: tmp_m((p+1)*(p+1)), tmp_m2((p+1)*(p+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If centers are the same just scale
        if (c(3) .eq. 0) then
            call fmm_m2m_scale_adj(src_r, dst_r, p, alpha, src_m, beta, dst_m)
        ! Otherwise apply ztranslate matrix
        else
            call fmm_m2m_ztranslate_use_mat_adj(p, ztranslate_mat, alpha, &
                & src_m, beta, dst_m)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_transform_use_mat(p, transform_mat, alpha, src_m, zero, tmp_m)
    call fmm_m2m_ztranslate_use_mat_adj(p, ztranslate_mat, one, tmp_m, zero, &
        & tmp_m2)
    call fmm_sph_transform_use_mat(p, transform_mat, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2m_reflection_use_mat_adj

!> Direct L2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of local-to-local
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_l2l_ztranslate(z, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, indj, indn
    ! Init output
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    ! If harmonics have different centers
    if (z .ne. 0) then
        r1 = z / src_r
        r2 = dst_r / z
        pow_r1(1) = 1
        pow_r2(1) = 1
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = alpha * pow_r2(j+1) / vfact(j-k+1) / vfact(j+k+1) * &
                    & vscales(indj)
                do n = j, p
                    indn = n*n + n + 1
                    tmp2 = tmp1 * pow_r1(n+1) / vscales(indn) * &
                        & vfact(n-k+1) * vfact(n+k+1) / vfact(n-j+1) / &
                        & vfact(n-j+1)
                    if (mod(n+j, 2) .eq. 1) then
                        tmp2 = -tmp2
                    end if
                    if (k .eq. 0) then
                        dst_l(indj) = dst_l(indj) + tmp2*src_l(indn)
                    else
                        dst_l(indj+k) = dst_l(indj+k) + tmp2*src_l(indn+k)
                        dst_l(indj-k) = dst_l(indj-k) + tmp2*src_l(indn-k)
                    end if
                end do
            end do
        end do
    ! If harmonics are located at the same point
    else
        r1 = dst_r / src_r
        tmp1 = alpha
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_l(k) = dst_l(k) + src_l(k)*tmp1
            end do
            tmp1 = tmp1 * r1
        end do
    end if
end subroutine fmm_l2l_ztranslate

!> Adjoint L2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L^\top \f$ is an adjoint matrix of
!! local-to-local translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @parma[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_l2l_ztranslate_adj(z, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, indj, indn
    ! Init output
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    ! If harmonics have different centers
    if (z .ne. 0) then
        r1 = -z / dst_r
        r2 = -src_r / z
        pow_r1(1) = one
        pow_r2(1) = one
        do j = 2, p+1
            pow_r1(j) = pow_r1(j-1) * r1
            pow_r2(j) = pow_r2(j-1) * r2
        end do
        do j = 0, p
            indj = j*j + j + 1
            do k = 0, j
                tmp1 = alpha * pow_r2(j+1) / vfact(j-k+1) / vfact(j+k+1) * &
                    & vscales(indj)
                do n = j, p
                    indn = n*n + n + 1
                    tmp2 = tmp1 * pow_r1(n+1) / vscales(indn) * &
                        & vfact(n-k+1) * vfact(n+k+1) / vfact(n-j+1) / &
                        & vfact(n-j+1)
                    if (mod(n+j, 2) .eq. 1) then
                        tmp2 = -tmp2
                    end if
                    if (k .eq. 0) then
                        dst_l(indn) = dst_l(indn) + tmp2*src_l(indj)
                    else
                        dst_l(indn+k) = dst_l(indn+k) + tmp2*src_l(indj+k)
                        dst_l(indn-k) = dst_l(indn-k) + tmp2*src_l(indj-k)
                    end if
                end do
            end do
        end do
    ! If harmonics are located at the same point
    else
        r1 = src_r / dst_r
        tmp1 = alpha
        do j = 0, p
            indj = j*j + j + 1
            do k = indj-j, indj+j
                dst_l(k) = dst_l(k) + src_l(k)*tmp1
            end do
            tmp1 = tmp1 * r1
        end do
    end if
end subroutine fmm_l2l_ztranslate_adj

!> Save matrix of L2L translation along OZ axis
!!
!! In a case input `z` is zero no translation matrix is computed, as @ref
!! fmm_l2l_scale shall be used in this case without any precomputed matrices
!!
!!
!! @param[in] z: the OZ coordinate of the radius-vector from new to old centers
!!      of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$
!! @param[in] vfact: Square roots of factorials
!! @param[in] mat: Translation matrix for spherical harmonics
subroutine fmm_l2l_ztranslate_get_mat(z, src_r, dst_r, p, vscales, vfact, mat)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((p+1)*(p+1)), vfact(2*p+1)
    integer, intent(in) :: p
    ! Output
    real(dp), intent(out) :: mat((p+1)*(p+2)*(p+3)/6)
    ! Local variables
    real(dp) :: r1, r2, tmp1, tmp2, pow_r1(p+1), pow_r2(p+1)
    integer :: j, k, n, indj, indn, indmat
    if (z .eq. zero) then
        ! Fill matrix with NaNs
        tmp1 = zero
        mat = tmp1 / tmp1
        return
    end if
    r1 = z / src_r
    r2 = dst_r / z
    pow_r1(1) = one
    pow_r2(1) = one
    do j = 2, p+1
        pow_r1(j) = pow_r1(j-1) * r1
        pow_r2(j) = pow_r2(j-1) * r2
    end do
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        do k = 0, j
            tmp1 = pow_r2(j+1) / vfact(j-k+1) / vfact(j+k+1) * vscales(indj)
            do n = j, p
                indn = n*n + n + 1
                mat(indmat) = tmp1 * pow_r1(n+1) / vscales(indn) * &
                    & vfact(n-k+1) * vfact(n+k+1) / vfact(n-j+1) / &
                    & vfact(n-j+1)
                if (mod(n+j, 2) .eq. 1) then
                    mat(indmat) = -mat(indmat)
                end if
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_l2l_ztranslate_get_mat

!> Apply matrix of L2L translation along OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of local-to-local
!! translation over OZ axis.
!!
!! This function shall not be used in case if centers of spherical harmonics
!! are identical. The OZ translation matrix is not computed by @ref
!! fmm_l2l_ztranslate_get_mat in this case as @ref fmm_l2l_scale is intended to
!! treat the case.
!!
!!
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] mat: The OZ translation matrix for spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Multipole expansion in new harmonics
subroutine fmm_l2l_ztranslate_use_mat(p, mat, alpha, src_l, beta, dst_l)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: mat((p+1)*(p+2)*(p+3)/6), alpha, &
        & src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    integer :: j, k, n, indj, indn, indmat
    real(dp) :: tmp1, tmp2
    ! Cycle over matrix elements stored in a sparse way
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        ! Init output properly
        if (beta .eq. zero) then
            dst_l(indj-j:indj+j) = zero
        else
            dst_l(indj-j:indj+j) = beta * dst_l(indj-j:indj+j)
        end if
        ! k = 0
        tmp1 = zero
        do n = j, p
            indn = n*n + n + 1
            tmp1 = tmp1 + mat(indmat)*src_l(indn)
            indmat = indmat + 1
        end do
        dst_l(indj) = dst_l(indj) + alpha*tmp1
        ! k > 0
        do k = 1, j
            tmp1 = zero
            tmp2 = zero
            do n = j, p
                indn = n*n + n + 1
                tmp1 = tmp1 + mat(indmat)*src_l(indn+k)
                tmp2 = tmp2 + mat(indmat)*src_l(indn-k)
                indmat = indmat + 1
            end do
            dst_l(indj+k) = dst_l(indj+k) + alpha*tmp1
            dst_l(indj-k) = dst_l(indj-k) + alpha*tmp2
        end do
    end do
end subroutine fmm_l2l_ztranslate_use_mat

!> Apply adjoint matrix of L2L translation along OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L^\top \f$ is an adjoint matrix of a
!! local-to-local translation over OZ axis.
!!
!! This function shall not be used in case if centers of spherical harmonics
!! are identical. The OZ translation matrix is not computed by @ref
!! fmm_l2l_ztranslate_get_mat in this case as @ref fmm_l2l_scale_adj is
!! intended to treat the case.
!!
!!
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] mat: The OZ translation matrix for spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion in new harmonics
!!
!!
!! TODO: improve performance
subroutine fmm_l2l_ztranslate_use_mat_adj(p, mat, alpha, src_l, beta, dst_l)
    ! Inputs
    integer, intent(in) :: p
    real(dp), intent(in) :: mat((p+1)*(p+2)*(p+3)/6), alpha, &
        & src_l((p+1)*(p+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    integer :: j, k, n, indj, indn, indmat
    real(dp) :: tmp1, tmp2
    ! Init output properly
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    ! Cycle over matrix elements stored in a sparse way
    indmat = 1
    do j = 0, p
        indj = j*j + j + 1
        ! k = 0
        do n = j, p
            indn = n*n + n + 1
            dst_l(indn) = dst_l(indn) + alpha*mat(indmat)*src_l(indj)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            do n = j, p
                indn = n*n + n + 1
                dst_l(indn+k) = dst_l(indn+k) + alpha*mat(indmat)*src_l(indj+k)
                dst_l(indn-k) = dst_l(indn-k) + alpha*mat(indmat)*src_l(indj-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_l2l_ztranslate_use_mat_adj

!> Scale L2L, when spherical harmonics are centered in the same point
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of a local-to-local
!! translation when harmonics are centered in the same point.
!!
!!
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_l: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_l: Multipole expansion in new harmonics
subroutine fmm_l2l_scale(src_r, dst_r, p, alpha, src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: src_r, dst_r, alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, tmp1
    integer :: j, k, indj
    ! Init output
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    r1 = dst_r / src_r
    tmp1 = alpha
    do j = 0, p
        indj = j*j + j + 1
        do k = indj-j, indj+j
            dst_l(k) = dst_l(k) + src_l(k)*tmp1
        end do
        tmp1 = tmp1 * r1
    end do
end subroutine fmm_l2l_scale

!> Adjoint scale L2L, when spherical harmonics are centered in the same point
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L^\top \f$ is an adjoint matrix of a
!! local-to-local translation when harmonics are centered in the same point.
!!
!!
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_l: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_l: Multipole expansion in new harmonics
subroutine fmm_l2l_scale_adj(src_r, dst_r, p, alpha, src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: src_r, dst_r, alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: r1, tmp1
    integer :: j, k, indj
    ! Init output
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    r1 = src_r / dst_r
    tmp1 = alpha
    do j = 0, p
        indj = j*j + j + 1
        do k = indj-j, indj+j
            dst_l(k) = dst_l(k) + src_l(k)*tmp1
        end do
        tmp1 = tmp1 * r1
    end do
end subroutine fmm_l2l_scale_adj

!> Direct L2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L \f$ is a matrix of a local-to-local
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_l2l_rotation(c, src_r, dst_r, p, vscales, vfact, alpha, src_l, &
        & beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1), &
        & vmsin(p+1), tmp_l((p+1)*(p+1)), r1xz(2, 2), tmp_l2((p+1)*(p+1))
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. zero) then
        call fmm_l2l_ztranslate(c(3), src_r, dst_r, p, vscales, vfact, alpha, &
            & src_l, beta, dst_l)
        return
    end if
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis
    vmsin = -vsin
    call fmm_sph_rotate_oz(p, vcos, vmsin, alpha, src_l, zero, tmp_l)
    ! Perform rotation in the OXZ plane
    r1xz(1, 1) = ctheta
    r1xz(1, 2) = -stheta
    r1xz(2, 1) = stheta
    r1xz(2, 2) = ctheta
    call fmm_sph_transform_oxz(p, r1xz, one, tmp_l, zero, tmp_l2)
    ! OZ translation
    call fmm_l2l_ztranslate(rho, src_r, dst_r, p, vscales, vfact, one, &
        & tmp_l2, zero, tmp_l)
    ! Backward rotation in the OXZ plane
    r1xz(1, 2) = stheta
    r1xz(2, 1) = -stheta
    call fmm_sph_transform_oxz(p, r1xz, one, tmp_l, zero, tmp_l2)
    ! Backward rotation around OZ axis
    call fmm_sph_rotate_oz(p, vcos, vsin, one, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_rotation

!> Adjoint L2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_L^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_L^\top \f$ is an adjoint matrix of a
!! local-to-local translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_l2l_rotation_adj(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, vcos(p+1), vsin(p+1), &
        & vmsin(p+1), tmp_l((p+1)*(p+1)), r1xz(2, 2), tmp_l2((p+1)*(p+1))
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. zero) then
        call fmm_l2l_ztranslate_adj(c(3), src_r, dst_r, p, vscales, vfact, &
            & alpha, src_l, beta, dst_l)
        return
    end if
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, p, vcos, vsin)
    ! Rotate around OZ axis
    vmsin = -vsin
    call fmm_sph_rotate_oz(p, vcos, vmsin, alpha, src_l, zero, tmp_l)
    ! Perform rotation in the OXZ plane
    r1xz(1, 1) = ctheta
    r1xz(1, 2) = -stheta
    r1xz(2, 1) = stheta
    r1xz(2, 2) = ctheta
    call fmm_sph_transform_oxz(p, r1xz, one, tmp_l, zero, tmp_l2)
    ! OZ translation
    call fmm_l2l_ztranslate_adj(rho, src_r, dst_r, p, vscales, vfact, one, &
        & tmp_l2, zero, tmp_l)
    ! Backward rotation in the OXZ plane
    r1xz(1, 2) = stheta
    r1xz(2, 1) = -stheta
    call fmm_sph_transform_oxz(p, r1xz, one, tmp_l, zero, tmp_l2)
    ! Backward rotation around OZ axis
    call fmm_sph_rotate_oz(p, vcos, vsin, one, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_rotation_adj

!> Direct L2L translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_l2l_reflection(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: max123, ssq123, rho, c1(3), c1_norm, r1(3, 3), &
        & tmp_l((p+1)*(p+1)), tmp_l2((p+1)*(p+1))
    integer :: m, n
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        call fmm_l2l_ztranslate(c(3), src_r, dst_r, p, vscales, vfact, alpha, &
            & src_l, beta, dst_l)
        return
    end if
    ! Compute rho, c(1:2) != zero already
    if (abs(c(1)) .gt. abs(c(2))) then
        max123 = abs(c(1))
        ssq123 = one + (c(2)/c(1))**2
    else
        max123 = abs(c(2))
        ssq123 = one + (c(1)/c(2))**2
    end if
    if (abs(c(3)) .gt. max123) then
        rho = one + ssq123*(max123/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq123 + (c(3)/max123)**2
        rho = max123 * sqrt(rho)
    end if
    ! Reorder (x,y,z) -> (y,z,x) and set proper reflection normal vector c1
    if (c(3) .ge. zero) rho = -rho
    c1(1) = c(2)
    c1(2) = c(3) - rho ! -rho and c(3) have the same sign
    c1(3) = c(1)
    ! Normalize vector c1. We know for sure c1(2) is maximum and is non-zero
    max123 = abs(c1(2))
    ssq123 = one + (c1(1)/c1(2))**2 + (c1(3)/c1(2))**2
    c1_norm = max123 * sqrt(ssq123)
    c1 = c1 / c1_norm
    r1 = 0
    do m = 1, 3
        r1(m, m) = 1
        do n = 1, 3
            r1(n, m) = r1(n, m) - two*c1(n)*c1(m)
        end do
    end do
    call fmm_sph_transform(p, r1, alpha, src_l, zero, tmp_l)
    call fmm_l2l_ztranslate(rho, src_r, dst_r, p, vscales, vfact, one, tmp_l, &
        & zero, tmp_l2)
    call fmm_sph_transform(p, r1, one, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_reflection

!> Adjoint L2L translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_l2l_reflection_adj(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: max123, ssq123, rho, c1(3), c1_norm, r1(3, 3), &
        & tmp_l((p+1)*(p+1)), tmp_l2((p+1)*(p+1))
    integer :: m, n
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        call fmm_l2l_ztranslate_adj(c(3), src_r, dst_r, p, vscales, vfact, &
            & alpha, src_l, beta, dst_l)
        return
    end if
    ! Compute rho, c(1:2) != zero already
    if (abs(c(1)) .gt. abs(c(2))) then
        max123 = abs(c(1))
        ssq123 = one + (c(2)/c(1))**2
    else
        max123 = abs(c(2))
        ssq123 = one + (c(1)/c(2))**2
    end if
    if (abs(c(3)) .gt. max123) then
        rho = one + ssq123*(max123/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq123 + (c(3)/max123)**2
        rho = max123 * sqrt(rho)
    end if
    ! Reorder (x,y,z) -> (y,z,x) and set proper reflection normal vector c1
    if (c(3) .ge. zero) rho = -rho
    c1(1) = c(2)
    c1(2) = c(3) - rho ! -rho and c(3) have the same sign
    c1(3) = c(1)
    ! Normalize vector c1. We know for sure c1(2) is maximum and is non-zero
    max123 = abs(c1(2))
    ssq123 = one + (c1(1)/c1(2))**2 + (c1(3)/c1(2))**2
    c1_norm = max123 * sqrt(ssq123)
    c1 = c1 / c1_norm
    r1 = 0
    do m = 1, 3
        r1(m, m) = 1
        do n = 1, 3
            r1(n, m) = r1(n, m) - two*c1(n)*c1(m)
        end do
    end do
    call fmm_sph_transform(p, r1, alpha, src_l, zero, tmp_l)
    call fmm_l2l_ztranslate_adj(rho, src_r, dst_r, p, vscales, vfact, one, &
        & tmp_l, zero, tmp_l2)
    call fmm_sph_transform(p, r1, one, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_reflection_adj

!> Save matrices of L2L operation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[out] transform_mat: Matrix of reflection
!! @param[out] ztranslate_mat: Matrix of OZ translation
subroutine fmm_l2l_reflection_get_mat(c, src_r, dst_r, p, vscales, vfact, &
        & transform_mat, ztranslate_mat)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1)
    integer, intent(in) :: p
    ! Outputs
    real(dp), intent(out) :: transform_mat((p+1)*(2*p+1)*(2*p+3)/6)
    real(dp), intent(out) :: ztranslate_mat((p+1)*(p+2)*(p+3)/6)
    ! Local variables
    real(dp) :: rho, r1(3, 3)
    ! If no need for transformation, just generate translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If even z coordinate is zero do nothing, as this case requires a
        ! simple scaling
        if (c(3) .ne. zero) then
            call fmm_l2l_ztranslate_get_mat(c(3), src_r, dst_r, p, vscales, &
                & vfact, ztranslate_mat)
        end if
        return
    end if
    call coord_reflect_get_mat(c, rho, r1)
    call fmm_sph_transform_get_mat(p, r1, transform_mat)
    call fmm_l2l_ztranslate_get_mat(rho, src_r, dst_r, p, vscales, vfact, &
        & ztranslate_mat)
end subroutine fmm_l2l_reflection_get_mat

!> Apply matrices of L2L translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] transform_mat: Matrix of a reflection
!! @param[in] ztranslate_mat: Matrix of an OZ translation
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Multipole expansion of old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Multipole expansion of new harmonics
subroutine fmm_l2l_reflection_use_mat(c, src_r, dst_r, p, transform_mat, &
        & ztranslate_mat, alpha, src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, &
        & transform_mat((p+1)*(2*p+1)*(2*p+3)/3), &
        & ztranslate_mat((p+1)*(p+2)*(p+3)/6), &
        & alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: tmp_l((p+1)*(p+1)), tmp_l2((p+1)*(p+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If centers are the same just scale
        if (c(3) .eq. 0) then
            call fmm_l2l_scale(src_r, dst_r, p, alpha, src_l, beta, dst_l)
        ! Otherwise apply ztranslate matrix
        else
            call fmm_l2l_ztranslate_use_mat(p, ztranslate_mat, alpha, src_l, &
                & beta, dst_l)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_transform_use_mat(p, transform_mat, alpha, src_l, zero, tmp_l)
    call fmm_l2l_ztranslate_use_mat(p, ztranslate_mat, one, tmp_l, zero, &
        & tmp_l2)
    call fmm_sph_transform_use_mat(p, transform_mat, one, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_reflection_use_mat

!> Adjoint apply matrices of L2L translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] transform_mat: Matrix of a reflection
!! @param[in] ztranslate_mat: Matrix of an OZ translation
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Multipole expansion of old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Multipole expansion of new harmonics
subroutine fmm_l2l_reflection_use_mat_adj(c, src_r, dst_r, p, transform_mat, &
        & ztranslate_mat, alpha, src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, &
        & transform_mat((p+1)*(2*p+1)*(2*p+3)/3), &
        & ztranslate_mat((p+1)*(p+2)*(p+3)/6), &
        & alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: tmp_l((p+1)*(p+1)), tmp_l2((p+1)*(p+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If centers are the same just scale
        if (c(3) .eq. 0) then
            call fmm_l2l_scale_adj(src_r, dst_r, p, alpha, src_l, beta, dst_l)
        ! Otherwise apply ztranslate matrix
        else
            call fmm_l2l_ztranslate_use_mat_adj(p, ztranslate_mat, alpha, &
                & src_l, beta, dst_l)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_transform_use_mat(p, transform_mat, alpha, src_l, zero, tmp_l)
    call fmm_l2l_ztranslate_use_mat_adj(p, ztranslate_mat, one, tmp_l, zero, &
        & tmp_l2)
    call fmm_sph_transform_use_mat(p, transform_mat, one, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_reflection_use_mat_adj

!> Direct L2L translation by 2 reflections and 1 translation
!!
!! Slightly shorter implementation of @ref fmm_l2l_reflection
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_l2l_reflection2(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: rho, r1(3, 3), tmp_l((p+1)*(p+1)), tmp_l2((p+1)*(p+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        call fmm_l2l_ztranslate(c(3), src_r, dst_r, p, vscales, vfact, alpha, &
            & src_l, beta, dst_l)
        return
    end if
    call coord_reflect_get_mat(c, rho, r1)
    call fmm_sph_transform(p, r1, alpha, src_l, zero, tmp_l)
    call fmm_l2l_ztranslate(rho, src_r, dst_r, p, vscales, vfact, one, tmp_l, &
        & zero, tmp_l2)
    call fmm_sph_transform(p, r1, one, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_reflection2

!> Adjoint L2L translation by 2 reflections and 1 translation
!!
!! Slightly shorter implementation of @ref fmm_l2l_reflection_adj
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of spherical harmonics
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_m: Expansion in new harmonics
subroutine fmm_l2l_reflection2_adj(c, src_r, dst_r, p, vscales, vfact, alpha, &
        & src_l, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((p+1)*(p+1)), &
        & vfact(2*p+1), alpha, src_l((p+1)*(p+1)), beta
    integer, intent(in) :: p
    ! Output
    real(dp), intent(inout) :: dst_l((p+1)*(p+1))
    ! Local variables
    real(dp) :: rho, r1(3, 3), tmp_l((p+1)*(p+1)), tmp_l2((p+1)*(p+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        call fmm_l2l_ztranslate_adj(c(3), src_r, dst_r, p, vscales, vfact, &
            & alpha, src_l, beta, dst_l)
        return
    end if
    call coord_reflect_get_mat(c, rho, r1)
    call fmm_sph_transform(p, r1, alpha, src_l, zero, tmp_l)
    call fmm_l2l_ztranslate_adj(rho, src_r, dst_r, p, vscales, vfact, one, &
        & tmp_l, zero, tmp_l2)
    call fmm_sph_transform(p, r1, one, tmp_l2, beta, dst_l)
end subroutine fmm_l2l_reflection2_adj

!> Direct M2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M \f$ is a matrix of multipole-to-local
!! translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old (multipole) harmonics
!! @param[in] dst_r: Radius of new (local) harmonics
!! @parma[in] pm: Maximal degree of multipole spherical harmonics
!! @parma[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multipler for `src_m`
!! @param[in] src_m: Expansion in old (multipole) harmonics
!! @param[in] beta: Scalar multipler for `dst_l`
!! @param[inout] dst_l: Expansion in new (local) harmonics
subroutine fmm_m2l_ztranslate(z, src_r, dst_r, pm, pl, vscales, vfact, alpha, &
        & src_m, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((pm+pl+1)*(pm+pl+1)), &
        & vfact(2*(pm+pl)+1), alpha, src_m((pm+1)*(pm+1)), beta
    integer, intent(in) :: pm, pl
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Local variables
    real(dp) :: tmp1, tmp2, pow_r1(pm+1), pow_r2(pl+1), r1, r2
    integer :: j, k, n, indj, indn
    ! Init output
    if (beta .eq. zero) then
        dst_l = zero
    else
        dst_l = beta * dst_l
    end if
    ! z cannot be zero, as input sphere (multipole) must not intersect with
    ! output sphere (local)
    if (z .eq. 0) then
        return
    end if
    r1 = src_r / z
    r2 = dst_r / z
    ! This abs(r1) makes it possible to work with negative z to avoid
    ! unnecessary rotation to positive z
    pow_r1(1) = abs(r1)
    pow_r2(1) = one
    do j = 2, pm+1
        pow_r1(j) = pow_r1(j-1) * r1
    end do
    do j = 2, pl+1
        pow_r2(j) = pow_r2(j-1) * r2
    end do
    do j = 0, pl
        indj = j*j + j + 1
        do k = 0, j
            tmp1 = alpha * vscales(indj) * pow_r2(j+1) / vfact(j-k+1) / &
                & vfact(j+k+1)
            do n = k, pm
                indn = n*n + n + 1
                tmp2 = tmp1 * pow_r1(n+1) / vscales(indn) / vfact(n-k+1) / &
                    & vfact(n+k+1) * (vfact(j+n+1)**2)
                if (mod(n+k, 2) .eq. 1) then
                    tmp2 = -tmp2
                end if
                if (k .eq. 0) then
                    dst_l(indj) = dst_l(indj) + tmp2*src_m(indn)
                else
                    dst_l(indj+k) = dst_l(indj+k) + tmp2*src_m(indn+k)
                    dst_l(indj-k) = dst_l(indj-k) + tmp2*src_m(indn-k)
                end if
            end do
        end do
    end do
end subroutine fmm_m2l_ztranslate

!> Adjoint M2L translation over OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M^\top \f$ is an adjoint matrix of
!! multipole-to-local translation over OZ axis.
!!
!!
!! @param[in] z: OZ coordinate from new to old centers of harmonics
!! @param[in] src_r: Radius of old (local) harmonics
!! @param[in] dst_r: Radius of new (multipole) harmonics
!! @parma[in] pl: Maximal degree of local spherical harmonics
!! @parma[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] vscales: Normalization constants for harmonics
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multipler for `src_l`
!! @param[in] src_l: Expansion in old (local) harmonics
!! @param[in] beta: Scalar multipler for `dst_m`
!! @param[inout] dst_m: Expansion in new (multipole) harmonics
subroutine fmm_m2l_ztranslate_adj(z, src_r, dst_r, pl, pm, vscales, vfact, &
        & alpha, src_l, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((pm+pl+1)*(pm+pl+1)), &
        & vfact(2*(pm+pl)+1), alpha, src_l((pl+1)*(pl+1)), beta
    integer, intent(in) :: pl, pm
    ! Output
    real(dp), intent(inout) :: dst_m((pm+1)*(pm+1))
    ! Local variables
    real(dp) :: tmp1, tmp2, pow_r1(pm+1), pow_r2(pl+1), r1, r2
    integer :: j, k, n, indj, indn
    ! Init output
    if (beta .eq. zero) then
        dst_m = zero
    else
        dst_m = beta * dst_m
    end if
    ! z cannot be zero, as input sphere (multipole) must not intersect with
    ! output sphere (local)
    if (z .eq. 0) then
        return
    end if
    r1 = -dst_r / z
    r2 = -src_r / z
    ! This abs(r1) makes it possible to work with negative z to avoid
    ! unnecessary rotation to positive z
    pow_r1(1) = abs(r1)
    pow_r2(1) = one
    do j = 2, pm+1
        pow_r1(j) = pow_r1(j-1) * r1
    end do
    do j = 2, pl+1
        pow_r2(j) = pow_r2(j-1) * r2
    end do
    do j = 0, pl
        indj = j*j + j + 1
        do k = 0, j
            tmp1 = alpha * vscales(indj) * pow_r2(j+1) / vfact(j-k+1) / &
                & vfact(j+k+1)
            do n = k, pm
                indn = n*n + n + 1
                tmp2 = tmp1 * pow_r1(n+1) / vscales(indn) / vfact(n-k+1) / &
                    & vfact(n+k+1) * (vfact(j+n+1)**2)
                if (mod(n+k, 2) .eq. 1) then
                    tmp2 = -tmp2
                end if
                if (k .eq. 0) then
                    dst_m(indn) = dst_m(indn) + tmp2*src_l(indj)
                else
                    dst_m(indn+k) = dst_m(indn+k) + tmp2*src_l(indj+k)
                    dst_m(indn-k) = dst_m(indn-k) + tmp2*src_l(indj-k)
                end if
            end do
        end do
    end do
end subroutine fmm_m2l_ztranslate_adj

!> Save matrix of M2L translation along OZ axis
!!
!! In a case input `z` is zero no translation matrix is computed, as @ref
!! fmm_m2m_scale shall be used in this case without any precomputed matrices
!!
!!
!! @param[in] z: the OZ coordinate of the radius-vector from new to old centers
!!      of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$
!! @param[in] vfact: Square roots of factorials
!! @param[in] mat: Translation matrix for spherical harmonics
subroutine fmm_m2l_ztranslate_get_mat(z, src_r, dst_r, pm, pl, vscales, &
        & vfact, mat)
    ! Inputs
    real(dp), intent(in) :: z, src_r, dst_r, vscales((pm+pl+1)*(pm+pl+1)), &
        & vfact(2*(pm+pl)+1)
    integer, intent(in) :: pm, pl
    ! Output
    real(dp), intent(out) :: mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6)
    ! Local variables
    real(dp) :: r1, r2, tmp1, pow_r1(pm+1), pow_r2(pl+1)
    integer :: j, k, n, indj, indn, indjn, indmat
    if (z .eq. zero) then
        ! M and L cannot be located in the same place
        return
    end if
    r1 = src_r / z
    r2 = dst_r / z
    ! This abs(r1) makes it possible to work with negative z to avoid
    ! unnecessary rotation to positive z
    pow_r1(1) = abs(r1)
    pow_r2(1) = one
    do j = 2, pm+1
        pow_r1(j) = pow_r1(j-1) * r1
    end do
    do j = 2, pl+1
        pow_r2(j) = pow_r2(j-1) * r2
    end do
    indmat = 1
    do j = 0, pl
        indj = j*j + j + 1
        do k = 0, j
            tmp1 = vscales(indj) * pow_r2(j+1) / vfact(j-k+1) / vfact(j+k+1)
            do n = k, pm
                indn = n*n + n + 1
                mat(indmat) = tmp1 * pow_r1(n+1) / vscales(indn) / &
                    & vfact(n-k+1) / vfact(n+k+1) * vfact(j+n+1) * vfact(j+n+1)
                if (mod(n+k, 2) .eq. 1) then
                    mat(indmat) = -mat(indmat)
                end if
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2l_ztranslate_get_mat

!> Apply matrix of M2L translation along OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M \f$ is a matrix of multipole-to-multipole
!! translation over OZ axis.
!!
!! This function shall not be used in case if centers of spherical harmonics
!! are identical. The OZ translation matrix is not computed by @ref
!! fmm_m2l_ztranslate_get_mat in this case as @ref fmm_m2l_scale is intended to
!! treat the case.
!!
!!
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] mat: The OZ translation matrix for spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_m2l_ztranslate_use_mat(pm, pl, mat, alpha, src_m, beta, dst_l)
    ! Inputs
    integer, intent(in) :: pm, pl
    real(dp), intent(in) :: mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6), alpha, src_m((pm+1)*(pm+1)), beta
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Local variables
    integer :: j, k, n, indj, indn, indjn, indmat
    real(dp) :: tmp1, tmp2
    ! Cycle over matrix elements stored in a sparse way
    indmat = 1
    do j = 0, pl
        indj = j*j + j + 1
        ! Init output properly
        if (beta .eq. zero) then
            dst_l(indj-j:indj+j) = zero
        else
            dst_l(indj-j:indj+j) = beta * dst_l(indj-j:indj+j)
        end if
        ! k = 0
        tmp1 = zero
        do n = 0, pm
            indn = n*n + n + 1
            tmp1 = tmp1 + mat(indmat)*src_m(indn)
            indmat = indmat + 1
        end do
        dst_l(indj) = dst_l(indj) + alpha*tmp1
        ! k > 0
        do k = 1, j
            tmp1 = zero
            tmp2 = zero
            do n = k, pm
                indn = n*n + n + 1
                tmp1 = tmp1 + mat(indmat)*src_m(indn+k)
                tmp2 = tmp2 + mat(indmat)*src_m(indn-k)
                indmat = indmat + 1
            end do
            dst_l(indj+k) = dst_l(indj+k) + alpha*tmp1
            dst_l(indj-k) = dst_l(indj-k) + alpha*tmp2
        end do
    end do
end subroutine fmm_m2l_ztranslate_use_mat

!> Apply adjoint matrix of M2L translation along OZ axis
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M^\top \f$ is an adjoint matrix of a
!! multipole-to-local translation over OZ axis.
!!
!! This function shall not be used in case if centers of spherical harmonics
!! are identical. The OZ translation matrix is not computed by @ref
!! fmm_m2l_ztranslate_get_mat in this case as @ref fmm_m2l_scale_adj is
!! intended to treat the case.
!!
!!
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] mat: The OZ translation matrix for spherical harmonics
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_m: Multipole expansion in new harmonics
!!
!!
!! TODO: improve performance
subroutine fmm_m2l_ztranslate_use_mat_adj(pl, pm, mat, alpha, src_l, beta, &
        & dst_m)
    ! Inputs
    integer, intent(in) :: pl, pm
    real(dp), intent(in) :: mat((min(pm,pl)+1) * (min(pm,pl)+2) &
        & * (3*max(pm,pl)+3-min(pm,pl)) / 6), alpha, src_l((pl+1)*(pl+1)), &
        & beta
    ! Output
    real(dp), intent(out) :: dst_m((pm+1)*(pm+1))
    ! Local variables
    integer :: j, k, n, indj, indn, indjn, indmat
    ! Init output properly
    if (beta .eq. zero) then
        dst_m = zero
    else
        dst_m = beta * dst_m
    end if
    ! Cycle over matrix elements stored in a sparse way
    indmat = 1
    do j = 0, pl
        indj = j*j + j + 1
        ! k = 0
        do n = 0, pm
            indn = n*n + n + 1
            dst_m(indn) = dst_m(indn) + alpha*mat(indmat)*src_l(indj)
            indmat = indmat + 1
        end do
        ! k > 0
        do k = 1, j
            do n = k, pm
                indn = n*n + n + 1
                dst_m(indn+k) = dst_m(indn+k) + alpha*mat(indmat)*src_l(indj+k)
                dst_m(indn-k) = dst_m(indn-k) + alpha*mat(indmat)*src_l(indj-k)
                indmat = indmat + 1
            end do
        end do
    end do
end subroutine fmm_m2l_ztranslate_use_mat_adj

!> Direct M2L translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] p: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: expansion in new harmonics
subroutine fmm_m2l_reflection(c, src_r, dst_r, pm, pl, vscales, vfact, alpha, &
        & src_m, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & vfact(2*(pm+pl)+1), alpha, src_m((pm+1)*(pm+1)), beta
    integer, intent(in) :: pm, pl
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Local variables
    real(dp) :: max123, ssq123, rho, c1(3), c1_norm, r1(3, 3), &
        & tmp_m((pm+1)*(pm+1)), tmp_l((pl+1)*(pl+1))
    integer :: m, n
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        call fmm_m2l_ztranslate(c(3), src_r, dst_r, pm, pl, vscales, vfact, &
            & alpha, src_m, beta, dst_l)
        return
    end if
    ! Compute rho, c(1:2) != zero already
    if (abs(c(1)) .gt. abs(c(2))) then
        max123 = abs(c(1))
        ssq123 = one + (c(2)/c(1))**2
    else
        max123 = abs(c(2))
        ssq123 = one + (c(1)/c(2))**2
    end if
    if (abs(c(3)) .gt. max123) then
        rho = one + ssq123*(max123/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq123 + (c(3)/max123)**2
        rho = max123 * sqrt(rho)
    end if
    ! Reorder (x,y,z) -> (y,z,x) and set proper reflection normal vector c1
    if (c(3) .ge. zero) rho = -rho
    c1(1) = c(2)
    c1(2) = c(3) - rho ! -rho and c(3) have the same sign
    c1(3) = c(1)
    ! Normalize vector c1. We know for sure c1(2) is maximum and is non-zero
    max123 = abs(c1(2))
    ssq123 = one + (c1(1)/c1(2))**2 + (c1(3)/c1(2))**2
    c1_norm = max123 * sqrt(ssq123)
    c1 = c1 / c1_norm
    r1 = 0
    do m = 1, 3
        r1(m, m) = 1
        do n = 1, 3
            r1(n, m) = r1(n, m) - two*c1(n)*c1(m)
        end do
    end do
    call fmm_sph_transform(pm, r1, alpha, src_m, zero, tmp_m)
    call fmm_m2l_ztranslate(rho, src_r, dst_r, pm, pl, vscales, vfact, one, &
        & tmp_m, zero, tmp_l)
    call fmm_sph_transform(pl, r1, one, tmp_l, beta, dst_l)
end subroutine fmm_m2l_reflection

!> Adjoint M2L translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: expansion in new harmonics
subroutine fmm_m2l_reflection_adj(c, src_r, dst_r, pl, pm, vscales, vfact, &
        & alpha, src_l, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & vfact(2*(pm+pl)+1), alpha, src_l((pl+1)*(pl+1)), beta
    integer, intent(in) :: pl, pm
    ! Output
    real(dp), intent(inout) :: dst_m((pm+1)*(pm+1))
    ! Local variables
    real(dp) :: max123, ssq123, rho, c1(3), c1_norm, r1(3, 3), &
        & tmp_m((pm+1)*(pm+1)), tmp_l((pl+1)*(pl+1))
    integer :: m, n
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        call fmm_m2l_ztranslate_adj(c(3), src_r, dst_r, pl, pm, vscales, &
            & vfact, alpha, src_l, beta, dst_m)
        return
    end if
    ! Compute rho, c(1:2) != zero already
    if (abs(c(1)) .gt. abs(c(2))) then
        max123 = abs(c(1))
        ssq123 = one + (c(2)/c(1))**2
    else
        max123 = abs(c(2))
        ssq123 = one + (c(1)/c(2))**2
    end if
    if (abs(c(3)) .gt. max123) then
        rho = one + ssq123*(max123/c(3))**2
        rho = abs(c(3)) * sqrt(rho)
    else
        rho = ssq123 + (c(3)/max123)**2
        rho = max123 * sqrt(rho)
    end if
    ! Reorder (x,y,z) -> (y,z,x) and set proper reflection normal vector c1
    if (c(3) .ge. zero) rho = -rho
    c1(1) = c(2)
    c1(2) = c(3) - rho ! -rho and c(3) have the same sign
    c1(3) = c(1)
    ! Normalize vector c1. We know for sure c1(2) is maximum and is non-zero
    max123 = abs(c1(2))
    ssq123 = one + (c1(1)/c1(2))**2 + (c1(3)/c1(2))**2
    c1_norm = max123 * sqrt(ssq123)
    c1 = c1 / c1_norm
    r1 = 0
    do m = 1, 3
        r1(m, m) = 1
        do n = 1, 3
            r1(n, m) = r1(n, m) - two*c1(n)*c1(m)
        end do
    end do
    call fmm_sph_transform(pl, r1, alpha, src_l, zero, tmp_l)
    call fmm_m2l_ztranslate_adj(rho, src_r, dst_r, pl, pm, vscales, vfact, &
        & one, tmp_l, zero, tmp_m)
    call fmm_sph_transform(pm, r1, one, tmp_m, beta, dst_m)
end subroutine fmm_m2l_reflection_adj

!> Direct M2L translation by 2 reflections and 1 translation
!!
!! Slightly shorter implementation of @ref fmm_m2l_reflection
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Local expansion in new harmonics
subroutine fmm_m2l_reflection2(c, src_r, dst_r, pm, pl, vscales, vfact, &
        & alpha, src_m, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & vfact(2*(pm+pl)+1), alpha, src_m((pm+1)*(pm+1)), beta
    integer, intent(in) :: pm, pl
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Local variables
    real(dp) :: rho, r1(3, 3), tmp_m((pm+1)*(pm+1)), tmp_l((pl+1)*(pl+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        call fmm_m2l_ztranslate(c(3), src_r, dst_r, pm, pl, vscales, vfact, &
            & alpha, src_m, beta, dst_l)
        return
    end if
    call coord_reflect_get_mat(c, rho, r1)
    call fmm_sph_transform(pm, r1, alpha, src_m, zero, tmp_m)
    call fmm_m2l_ztranslate(rho, src_r, dst_r, pm, pl, vscales, vfact, one, &
        & tmp_m, zero, tmp_l)
    call fmm_sph_transform(pl, r1, one, tmp_l, beta, dst_l)
end subroutine fmm_m2l_reflection2

!> Adjoint M2L translation by 2 reflections and 1 translation
!!
!! Slightly shorter implementation of @ref fmm_m2l_reflection
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] vscales: Normalization constants for \f$ Y_\ell^m \f$
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Multipole expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion in new harmonics
subroutine fmm_m2l_reflection2_adj(c, src_r, dst_r, pl, pm, vscales, vfact, &
        & alpha, src_l, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & vfact(2*(pm+pl)+1), alpha, src_l((pl+1)*(pl+1)), beta
    integer, intent(in) :: pl, pm
    ! Output
    real(dp), intent(inout) :: dst_m((pm+1)*(pm+1))
    ! Local variables
    real(dp) :: rho, r1(3, 3), tmp_m((pm+1)*(pm+1)), tmp_l((pl+1)*(pl+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        call fmm_m2l_ztranslate_adj(c(3), src_r, dst_r, pl, pm, vscales, &
            & vfact, alpha, src_l, beta, dst_m)
        return
    end if
    call coord_reflect_get_mat(c, rho, r1)
    call fmm_sph_transform(pl, r1, alpha, src_l, zero, tmp_l)
    call fmm_m2l_ztranslate_adj(rho, src_r, dst_r, pl, pm, vscales, vfact, &
        & one, tmp_l, zero, tmp_m)
    call fmm_sph_transform(pm, r1, one, tmp_m, beta, dst_m)
end subroutine fmm_m2l_reflection2_adj

!> Save matrices of M2L operation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[out] transform_mat: Matrix of reflection
!! @param[out] ztranslate_mat: Matrix of OZ translation
subroutine fmm_m2l_reflection_get_mat(c, src_r, dst_r, pm, pl, vscales, &
        & vfact, transform_mat, ztranslate_mat)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & vfact(2*(pm+pl)+1)
    integer, intent(in) :: pm, pl
    ! Outputs
    real(dp), intent(out) :: transform_mat((max(pm,pl)+1)*(2*max(pm,pl)+1)* &
        & (2*max(pm,pl)+3)/6)
    real(dp), intent(out) :: ztranslate_mat((min(pm,pl)+1)*(min(pm,pl)+2)* &
        & (3*max(pm,pl)+3-min(pm,pl))/6)
    ! Local variables
    real(dp) :: rho, r1(3, 3)
    ! If no need for transformation, just generate translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If M and L harmonics are at the same place ignore M2L completely
        if (c(3) .ne. zero) then
            call fmm_m2l_ztranslate_get_mat(c(3), src_r, dst_r, pm, pl, &
                & vscales, vfact, ztranslate_mat)
        end if
        return
    end if
    call coord_reflect_get_mat(c, rho, r1)
    call fmm_sph_transform_get_mat(max(pm,pl), r1, transform_mat)
    call fmm_m2l_ztranslate_get_mat(rho, src_r, dst_r, pm, pl, vscales, &
        & vfact, ztranslate_mat)
end subroutine fmm_m2l_reflection_get_mat

!> Apply matrices of M2L translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] transform_mat: Matrix of a reflection
!! @param[in] ztranslate_mat: Matrix of an OZ translation
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Multipole expansion of old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Multipole expansion of new harmonics
subroutine fmm_m2l_reflection_use_mat(c, src_r, dst_r, pm, pl, transform_mat, &
        & ztranslate_mat, alpha, src_m, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, alpha, src_m((pm+1)*(pm+1)), &
        & beta
    integer, intent(in) :: pm, pl
    real(dp), intent(in) :: &
        & transform_mat((max(pm,pl)+1)*(2*max(pm,pl)+1)*(2*max(pm,pl)+3)/3), &
        & ztranslate_mat((min(pm,pl)+1)*(min(pm,pl)+2)* &
        & (3*max(pm,pl)+3-min(pm,pl))/6)
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Local variables
    real(dp) :: tmp_m((pm+1)*(pm+1)), tmp_l((pl+1)*(pl+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! Only apply ztranslate matrix if c is non-zero
        if (c(3) .ne. 0) then
            call fmm_m2l_ztranslate_use_mat(pm, pl, ztranslate_mat, alpha, &
                & src_m, beta, dst_l)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_transform_use_mat(pm, transform_mat, alpha, src_m, zero, &
        & tmp_m)
    call fmm_m2l_ztranslate_use_mat(pm, pl, ztranslate_mat, one, tmp_m, zero, &
        & tmp_l)
    call fmm_sph_transform_use_mat(pl, transform_mat, one, tmp_l, beta, dst_l)
end subroutine fmm_m2l_reflection_use_mat

!> Adjoint apply matrices of M2L translation by 2 reflections and 1 translation
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] transform_mat: Matrix of a reflection
!! @param[in] ztranslate_mat: Matrix of an OZ translation
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: Local expansion of old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: Multipole expansion of new harmonics
subroutine fmm_m2l_reflection_use_mat_adj(c, src_r, dst_r, pl, pm, &
        & transform_mat, ztranslate_mat, alpha, src_l, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, alpha, src_l((pl+1)*(pl+1)), &
        & beta
    integer, intent(in) :: pl, pm
    real(dp), intent(in) :: &
        & transform_mat((max(pm,pl)+1)*(2*max(pm,pl)+1)*(2*max(pm,pl)+3)/3), &
        & ztranslate_mat((min(pm,pl)+1)*(min(pm,pl)+2)* &
        & (3*max(pm,pl)+3-min(pm,pl))/6)
    ! Output
    real(dp), intent(inout) :: dst_m((pm+1)*(pm+1))
    ! Local variables
    real(dp) :: tmp_m((pm+1)*(pm+1)), tmp_l((pl+1)*(pl+1))
    ! If no need for transformation, just do translation along z
    if ((c(1) .eq. zero) .and. (c(2) .eq. zero)) then
        ! If centers are the same ignore M2L completely
        if (c(3) .ne. 0) then
            call fmm_m2l_ztranslate_use_mat_adj(pl, pm, ztranslate_mat, &
                & alpha, src_l, beta, dst_m)
        end if
        return
    end if
    ! Apply reflection->translation->reflection
    call fmm_sph_transform_use_mat(pl, transform_mat, alpha, src_l, zero, &
        & tmp_l)
    call fmm_m2l_ztranslate_use_mat_adj(pl, pm, ztranslate_mat, one, tmp_l, &
        & zero, tmp_m)
    call fmm_sph_transform_use_mat(pm, transform_mat, one, tmp_m, beta, dst_m)
end subroutine fmm_m2l_reflection_use_mat_adj

!> Direct M2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M \f$ is a matrix of a multipole-to-local
!! translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] vscales: Normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_m`
!! @param[in] src_m: Expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_l`
!! @param[inout] dst_l: Expansion in new harmonics
subroutine fmm_m2l_rotation(c, src_r, dst_r, pm, pl, vscales, vfact, alpha, &
        & src_m, beta, dst_l)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & vfact(2*(pm+pl)+1), alpha, src_m((pm+1)*(pm+1)), beta
    integer, intent(in) :: pm, pl
    ! Output
    real(dp), intent(inout) :: dst_l((pl+1)*(pl+1))
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, vcos(max(pm,pl)+1), &
        & vsin(max(pm,pl)+1), vmsin(max(pm,pl)+1), tmp_m((pm+1)*(pm+1)), &
        & tmp_m2((pm+1)*(pm+1)), r1xz(2, 2), tmp_l((pl+1)*(pl+1)), &
        & tmp_l2((pl+1)*(pl+1))
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. zero) then
        call fmm_m2l_ztranslate(c(3), src_r, dst_r, pm, pl, vscales, vfact, &
            & alpha, src_m, beta, dst_l)
        return
    end if
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, max(pm,pl), vcos, vsin)
    ! Rotate around OZ axis
    vmsin = -vsin
    call fmm_sph_rotate_oz(pm, vcos, vmsin, alpha, src_m, zero, tmp_m)
    ! Perform rotation in the OXZ plane
    r1xz(1, 1) = ctheta
    r1xz(1, 2) = -stheta
    r1xz(2, 1) = stheta
    r1xz(2, 2) = ctheta
    call fmm_sph_transform_oxz(pm, r1xz, one, tmp_m, zero, tmp_m2)
    ! OZ translation
    call fmm_m2l_ztranslate(rho, src_r, dst_r, pm, pl, vscales, vfact, one, &
        & tmp_m2, zero, tmp_l)
    ! Backward rotation in the OXZ plane
    r1xz(1, 2) = stheta
    r1xz(2, 1) = -stheta
    call fmm_sph_transform_oxz(pl, r1xz, one, tmp_l, zero, tmp_l2)
    ! Backward rotation around OZ axis
    call fmm_sph_rotate_oz(pl, vcos, vsin, one, tmp_l2, beta, dst_l)
end subroutine fmm_m2l_rotation

!> Adjoint M2L translation by 4 rotations and 1 translation
!!
!! Compute the following matrix-vector product:
!! \f[
!!      \mathrm{dst} = \beta \mathrm{dst} + \alpha L_M^\top \mathrm{src},
!! \f]
!! where \f$ \mathrm{dst} \f$ is a vector of coefficients of output spherical
!! harmonics, \f$ \mathrm{src} \f$ is a vector of coefficients of input
!! spherical harmonics and \f$ L_M^\top \f$ is an adjoint matrix of a
!! multipole-to-local translation.
!!
!! Rotates around OZ and OY axes, translates over OZ and then rotates back
!! around OY and OZ axes.
!!
!!
!! @param[in] c: Radius-vector from new to old centers of harmonics
!! @param[in] src_r: Radius of old harmonics
!! @param[in] dst_r: Radius of new harmonics
!! @param[in] pl: Maximal degree of local spherical harmonics
!! @param[in] pm: Maximal degree of multipole spherical harmonics
!! @param[in] vscales: normalization constants for Y_lm
!! @param[in] vfact: Square roots of factorials
!! @param[in] alpha: Scalar multiplier for `src_l`
!! @param[in] src_l: expansion in old harmonics
!! @param[in] beta: Scalar multiplier for `dst_m`
!! @param[inout] dst_m: expansion in new harmonics
subroutine fmm_m2l_rotation_adj(c, src_r, dst_r, pl, pm, vscales, vfact, &
        & alpha, src_l, beta, dst_m)
    ! Inputs
    real(dp), intent(in) :: c(3), src_r, dst_r, vscales((pm+pl+1)**2), &
        & vfact(2*(pm+pl)+1), alpha, src_l((pl+1)*(pl+1)), beta
    integer, intent(in) :: pl, pm
    ! Output
    real(dp), intent(inout) :: dst_m((pm+1)*(pm+1))
    ! Local variables
    real(dp) :: rho, ctheta, stheta, cphi, sphi, vcos(max(pm,pl)+1), &
        & vsin(max(pm,pl)+1), vmsin(max(pm,pl)+1), tmp_m((pm+1)*(pm+1)), &
        & tmp_m2((pm+1)*(pm+1)), r1xz(2, 2), tmp_l((pl+1)*(pl+1)), &
        & tmp_l2((pl+1)*(pl+1))
    ! Covert Cartesian coordinates into spherical
    call carttosph(c, rho, ctheta, stheta, cphi, sphi)
    ! If no need for rotations, just do translation along z
    if (stheta .eq. 0) then
        call fmm_m2l_ztranslate_adj(c(3), src_r, dst_r, pl, pm, vscales, &
            & vfact, alpha, src_l, beta, dst_m)
        return
    end if
    ! Compute arrays of cos and sin that are needed for rotations of harmonics
    call trgev(cphi, sphi, max(pm,pl), vcos, vsin)
    ! Rotate around OZ axis
    vmsin = -vsin
    call fmm_sph_rotate_oz(pl, vcos, vmsin, alpha, src_l, zero, tmp_l)
    ! Perform rotation in the OXZ plane
    r1xz(1, 1) = ctheta
    r1xz(1, 2) = -stheta
    r1xz(2, 1) = stheta
    r1xz(2, 2) = ctheta
    call fmm_sph_transform_oxz(pl, r1xz, one, tmp_l, zero, tmp_l2)
    ! OZ translation
    call fmm_m2l_ztranslate_adj(rho, src_r, dst_r, pl, pm, vscales, vfact, &
        & one, tmp_l2, zero, tmp_m)
    ! Backward rotation in the OXZ plane
    r1xz(1, 2) = stheta
    r1xz(2, 1) = -stheta
    call fmm_sph_transform_oxz(pm, r1xz, one, tmp_m, zero, tmp_m2)
    ! Backward rotation around OZ axis
    call fmm_sph_rotate_oz(pm, vcos, vsin, one, tmp_m2, beta, dst_m)
end subroutine fmm_m2l_rotation_adj

!> Build a recursive inertial binary tree
!!
!! Uses inertial bisection in a recursive manner until each leaf node has only
!! one input sphere inside. Number of tree nodes is always 2*nsph-1.
!!
!!
!! @param[in] nsph: Number of input spheres
!! @param[in] csph: Centers of input spheres
!! @param[in] rsph: Radii of input spheres
!! @param[out] order: Ordering of input spheres to make spheres of one cluster
!!      contiguous.
!! @param[out] cluster: Indexes in `order` array of the first and the last
!!      spheres of each cluster
!! @param[out] children: Indexes of the first and the last children nodes. If
!!      both indexes are equal this means there is only 1 child node and if
!!      value of the first child node is 0, there are no children nodes.
!! @param[out] parent: Parent of each cluster. Value 0 means there is no parent
!!      node which corresponds to the root node (with index 1).
!! @param[out] cnode: Center of a bounding sphere of each node
!! @param[out] rnode: Radius of a bounding sphere of each node
!! @param[out] snode: Array of leaf nodes containing input spheres
subroutine tree_rib_build(nsph, csph, rsph, order, cluster, children, parent, &
        & cnode, rnode, snode)
    ! Inputs
    integer, intent(in) :: nsph
    real(dp), intent(in) :: csph(3, nsph), rsph(nsph)
    ! Outputs
    integer, intent(out) :: order(nsph), cluster(2, 2*nsph-1), &
        & children(2, 2*nsph-1), parent(2*nsph-1), snode(nsph)
    real(dp), intent(out) :: cnode(3, 2*nsph-1), rnode(2*nsph-1)
    ! Local variables
    integer :: nclusters, i, j, n, s, e, div
    real(dp) :: r, r1, r2, c(3), c1(3), c2(3), d, maxc, ssqc
    !! At first construct the tree
    nclusters = 2*nsph - 1
    ! Init the root node
    cluster(1, 1) = 1
    cluster(2, 1) = nsph
    parent(1) = 0
    ! Index of the first unassigned node
    j = 2
    ! Init spheres ordering
    do i = 1, nsph
        order(i) = i
    end do
    ! Divide nodes until a single spheres
    do i = 1, nclusters
        ! The first sphere in i-th node
        s = cluster(1, i)
        ! The last sphere in i-th node
        e = cluster(2, i)
        ! Number of spheres in i-th node
        n = e - s + 1
        ! Divide only if there are 2 or more spheres
        if (n .gt. 1) then
            ! Use inertial bisection to reorder spheres and cut into 2 halfs
            call tree_rib_node_bisect(nsph, csph, n, order(s:e), div)
            ! Assign the first half to the j-th node
            cluster(1, j) = s
            cluster(2, j) = s + div - 1
            ! Assign the second half to the (j+1)-th node
            cluster(1, j+1) = s + div
            cluster(2, j+1) = e
            ! Update list of children of i-th node
            children(1, i) = j
            children(2, i) = j + 1
            ! Set parents of new j-th and (j+1)=th node
            parent(j) = i
            parent(j+1) = i
            ! Shift index of the first unassigned node
            j = j + 2
        ! Set information for a leaf node
        else
            ! No children nodes
            children(:, i) = 0
            ! i-th node contains sphere ind(s)
            snode(order(s)) = i
        end if
    end do
    !! Compute bounding spheres of each node of the tree
    ! Bottom-to-top pass over all nodes of the tree
    do i = nclusters, 1, -1
        ! In case of a leaf node just use corresponding input sphere as a
        ! bounding sphere of the node
        if (children(1, i) .eq. 0) then
            ! Get correct index of the corresponding input sphere
            j = order(cluster(1, i))
            ! Get corresponding center and radius
            cnode(:, i) = csph(:, j)
            rnode(i) = rsph(j)
        ! In case of a non-leaf node get minimal sphere that contains bounding
        ! spheres of children nodes
        else
            ! The first child
            j = children(1, i)
            c1 = cnode(:, j)
            r1 = rnode(j)
            ! The second child
            j = children(2, i)
            c2 = cnode(:, j)
            r2 = rnode(j)
            ! Distance between centers of bounding spheres of children nodes
            c = c1 - c2
            ! Compute distance by scale and sum of scaled squares
            maxc = max(abs(c(1)), abs(c(2)), abs(c(3)))
            if (maxc .eq. zero) then
                d = zero
            else
                d = (c(1)/maxc)**2 + (c(2)/maxc)**2 + (c(3)/maxc)**2
                d = maxc * sqrt(d)
            end if
            ! If sphere #2 is completely inside sphere #1 use the first sphere
            if (r1 .ge. (r2+d)) then
                c = c1
                r = r1
            ! If sphere #1 is completely inside sphere #2 use the second sphere
            else if(r2 .ge. (r1+d)) then
                c = c2
                r = r2
            ! Otherwise use special formula to find a minimal sphere
            else
                r = (r1+r2+d) / 2
                c = c2 + c/d*(r-r2)
            end if
            ! Assign computed bounding sphere
            cnode(:, i) = c
            rnode(i) = r
        end if
    end do
end subroutine tree_rib_build

!> Divide given cluster of spheres into two subclusters by inertial bisection
!!
!! @param[in] nsph: Number of all input spheres
!! @param[in] csph: Centers of all input spheres
!! @param[in] n: Number of spheres in a given cluster
!! @param[inout] order: Indexes of spheres in a given cluster. On exit, indexes
!!      `order(1:div)` correspond to the first subcluster and indexes
!!      `order(div+1:n)` correspond to the second subcluster.
!! @param[out] div: Break point of `order` array between two clusters.
subroutine tree_rib_node_bisect(nsph, csph, n, order, div)
    ! Inputs
    integer, intent(in) :: nsph, n
    real(dp), intent(in) :: csph(3, nsph)
    ! Outputs
    integer, intent(inout) :: order(n)
    integer, intent(out) :: div
    ! Local variables
    real(dp) :: c(3), tmp_csph(3, n), s(3)
    real(dp), allocatable :: work(:)
    external :: dgemm, dsyev, dgemv
    integer :: i, l, r, lwork, info, tmp_order(n), istat
    ! Get average coordinate
    c = zero
    do i = 1, n
        c = c + csph(:, order(i))
    end do
    c = c / n
    ! Get coordinates minus average in a contiguous array
    do i = 1, n
        tmp_csph(:, i) = csph(:, order(i)) - c
    end do
    !! Find right singular vectors
    ! Get proper size of temporary workspace
    lwork = -1
    call dgesvd('N', 'O', 3, n, tmp_csph, 3, s, tmp_csph, 3, tmp_csph, 3, &
        & s, lwork, info)
    lwork = s(1)
    allocate(work(lwork), stat=istat)
    if (istat .ne. 0) stop "allocation failed"
    ! Get right singular vectors
    call dgesvd('N', 'O', 3, n, tmp_csph, 3, s, tmp_csph, 3, tmp_csph, 3, &
        & work, lwork, info)
    if (info .ne. 0) stop "dgesvd did not converge"
    deallocate(work, stat=istat)
    if (istat .ne. 0) stop "deallocation failed"
    !! Sort spheres by sign of the above scalar product, which is equal to
    !! the leading right singular vector scaled by the leading singular value.
    !! However, we only care about sign, so we take into account only the
    !! leading right singular vector.
    ! First empty index from the left
    l = 1
    ! First empty index from the right
    r = n
    ! Cycle over values of the singular vector
    do i = 1, n
        ! Positive scalar products are moved to the beginning of `order`
        if (tmp_csph(1, i) .ge. 0) then
            tmp_order(l) = order(i)
            l = l + 1
        ! Negative scalar products are moved to the end of `order`
        else
            tmp_order(r) = order(i)
            r = r - 1
        end if
    end do
    ! Set divider and update order
    div = r
    order = tmp_order
end subroutine tree_rib_node_bisect

end module dd_core

