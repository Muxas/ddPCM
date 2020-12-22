!> @copyright (c) 2020-2020 RWTH Aachen. All rights reserved.
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
real(kind=rp), parameter :: zero = 0d0, one = 1d0, two = 2d0, four = 4d0
real(kind=rp), parameter :: sqrt2 = sqrt(two)
real(kind=rp), parameter :: pi4 = atan(one)
real(kind=rp), parameter :: pi = four * pi4
real(kind=rp), parameter :: twopi = two * pi
real(kind=rp), parameter :: sqrt4pi = four * sqrt(atan(one))
!> Number of supported Lebedev-Laikov grids
integer, parameter :: nllg = 32
!> Number of grid points of each Lebedev-Laikov grid
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
    real(kind=rp), allocatable :: csph(:, :)
    !> Array of radii of atoms of dimension (nsph)
    real(kind=rp), allocatable :: rsph(:)
    !> Dielectric permittivity inside
    real(kind=rp) :: epsin
    !> Dielectric permittivity outside
    real(kind=rp) :: epsout
    !> Parameter kappa (meaning?)
    real(kind=rp) :: kappa
    !> Regularization parameter
    real(kind=rp) :: eta
    !> Relative error threshold for iterative solvers
    real(kind=rp) :: tol
    !> Maximal degree of modeling spherical harmonics
    integer :: lmax
    !> Number of Lebedev grid points on each sphere
    integer :: ngrid
    !> Whether to compute (1) or not (0) forces
    integer :: force
    !> Enable (1) or disable (0) use of FMM techniques
    integer :: fmm
    !!!!!!!!!!!!!! Constants, initialized by ddinit
    !> Maximal degree of used real normalized spherical harmonics
    !!
    !! `dmax=lmax` in case `fmm=0` and `dmax=pm+pl` in case `fmm=1`. If FMM is
    !! used, its M2L operation requires computing spherical harmonics of all
    !! degrees up to `pm+pl`. If `force=1` then this parameter might be
    !! increased by 1 because certain implementations of analytical gradients
    !! might rely on spherical harmonics of +1 degree.
    integer :: dmax
    !> Scales of real normalized spherical harmonics of degree up to dmax
    !!
    !! This array has `(dmax+1)**2` entries
    real(kind=rp), allocatable :: vscales(:)
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
    !> Values of spherical harmonics at Lebedev grid points
    !!
    !! Dimensions of this array are (vgrid_dmax, ngrid)
    real(kind=rp), allocatable :: vgrid(:, :)
    !> Maximum sqrt of factorial precomputed
    !!
    !! Just like with `dmax` parameter, number of used factorials is either
    !! `2*lmax+1` or `2*(pm+pl)+1` depending on whether FMM is used or not
    integer :: nfact
    !> Array of square roots of factorials of dimension (nfact)
    real(kind=rp), allocatable :: fact(:)
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
    real*8, allocatable :: cnode(:, :)
    ! Radius of bounding sphere of each cluster
    real*8, allocatable :: rnode(:)
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

!> Initialize ddX input and parameters
!!
!!
!! @param[in] n: Number of atoms
!! @param[in] x: \f$ x \f$ coordinates of atoms. Dimension is `(n)`
!! @param[in] y: \f$ y \f$ coordinates of atoms. Dimension is `(n)`
!! @param[in] z: \f$ z \f$ coordinates of atoms. Dimension is `(n)`
!! @param[in] rvdw: Van-der-Waals radii of atoms. Dimension is `(n)`
!! @param[in] model: Choose model: 1 for COSMO, 2 for PCM and 3 for LPB
!! @param[inout] ngrid: Approximate number of Lebedev grid points on input and
!!      actual number of grid points on exit. `ngrid` > 0
!! @param[in] lmax: Maximal degree of modeling spherical harmonics. `lmax` >= 0
!! @param[in] force: 1 if forces are required and 0 otherwise
!! @param[in] fmm: 1 to use FMM acceleration and 0 otherwise
!! @param[in] pm: Maximal degree of multipole spherical harmonics. `pm` >= 0
!! @param[in] pl: Maximal degree of local spherical harmonics. `pl` >= 0
!! @param[out] dd_data: Object containing all inputs
subroutine ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pm, pl, &
        & dd_data)
    ! Inputs
    integer, intent(in) :: n, model, lmax, force, fmm, pm, pl
    real(kind=rp) :: x(n), y(n), z(n), rvdw(n)
    ! Output
    type(dd_data_type) :: dd_data
    ! Inouts
    integer, intent(inout) :: ngrid
    ! Local variables
    integer :: istat, i, inear, jnear, igrid, jgrid
    !!! Set input parameters
    dd_data % nsph = n
    allocate(dd_data % csph(3, n), dd_data % rsph(n), stat=istat)
    if (istat .ne. 0) then
        write(*, *) "ddinit: [1] allocation failed!"
        stop
    end if
    dd_data % csph(1, :) = x
    dd_data % csph(2, :) = y
    dd_data % csph(3, :) = z
    dd_data % rsph = rvdw
    if ((model .lt. 1) .and. (model .gt. 3)) then
        write(*, *) "ddinit: wrong value of parameter `model`"
        stop
    end if
    dd_data % model = model
    if (lmax .lt. 0) then
        write(*, *) "ddinit: wrong value of parameter `lmax`"
        stop
    end if
    dd_data % lmax = lmax
    if ((force .lt. 0) .and. (force .gt. 1)) then
        write(*, *) "ddinit: wrong value of parameter `force`"
        stop
    end if
    dd_data % force = force
    if ((fmm .lt. 0) .and. (fmm .gt. 1)) then
        write(*, *) "ddinit: wrong value of parameter `fmm`"
        stop
    end if
    dd_data % fmm = fmm
    ! Set FMM parameters
    if(fmm .eq. 1) then
        if(pm .lt. 0) then
            write(*, *) "ddinit: wrong value of parameter `pm`"
            stop
        end if
        if(pl .lt. 0) then
            write(*, *) "ddinit: wrong value of parameter `pl`"
            stop
        end if
    end if
    !!! Generate constants
    ! At first compute sizes of auxiliary arrays for `fmm=0`
    dd_data % dmax = lmax
    dd_data % nfact = max(2*lmax+1, 2)
    ! Update sizes of arrays if `fmm=1`
    if(fmm .eq. 1) then
        dd_data % dmax = max(pm+pl, dd_data % lmax)
        dd_data % nfact = max(2*(pm+pl)+1, dd_data % nfact)
    end if
    ! Compute scaling factors of spherical harmonics
    call ylmscale(dd_data % dmax, dd_data % vscales)
    ! Precompute square roots of factorials
    dd_data % fact(1) = 1
    do i = 2, dd_data % nfact
        dd_data % fact(i) = dd_data % fact(i-1) * sqrt(dble(i-1))
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
    call llgrid(ngrid, dd_data % wgrid, dd_data % cgrid)
end subroutine ddinit

!> Deallocate object with corresponding data
!!
!! @param[inout] dd_data
subroutine ddfree(dd_data)
    ! Input/output
    type(dd_data_type) :: dd_data
    ! Local variables
    integer :: istat
    deallocate(dd_data % csph, dd_data % rsph, stat=istat)
    if (istat .ne. 0) then
        write(*, *) "ddfree: [1] deallocation failed!"
        stop
    end if
end subroutine ddfree

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
        rho = x(3) * sqrt(rho)
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
!! varies in \f$ [1-eta/2, 1+\eta/2] \f$ which is an centered shift.
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
    x = t - (s + 1.d0)*eta / 2.d0
    ! Lower bound of switch region
    flow = one - eta
    ! Define switch function \chi
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

