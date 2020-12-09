module ddpcm_lib
use ddcosmo, only: nbasis, nsph, ngrid, ncav, lmax, iconv, iprint, &
    & wghpot, intrhs, prtsph, zero, pt5, one, two, four, pi, basis, &
    & eps, csph, rsph, grid, w, ui, ndiis, sprod, ylmbas, facl, ddmkxi, &
    & ptcart, fdoka, fdokb, fdoga, nsph
use pcm_fmm

implicit none

real*8, allocatable :: rx_prc(:,:,:)
real*8, allocatable :: rhs(:,:), phieps(:,:), xs(:,:)
real*8, allocatable :: s(:,:), y(:,:), q(:,:)
logical :: dodiag
! Measure time
real*8 :: start_time, finish_time
! FMM-related global variables. I myself avoid global variables, but this
! program is written in this non-portable style. IT will take long time to
! change it in entire proglem, so I just follow it for now.
! Maximum degree of spherical harmonics for M (multipole) expansion
integer :: pm
! Maximum degree of spherical harmonics for L (local) expansion
integer :: pl
! Scalar factors for M2M, M2L and L2L operations
real*8, allocatable :: vscales(:)
! Values of L (local) spherical harmonics in grid points
real*8, allocatable :: vgrid(:, :)
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

contains

  subroutine ddpcm_init()
  ! initialize ddpcm module by allocating the preconditioner and
  ! various arrays, then build the preconditioner
  implicit none
  integer :: istatus
  allocate(rx_prc(nbasis,nbasis,nsph),s(nbasis,nsph),y(nbasis,nsph), &
    & xs(nbasis,nsph),phieps(nbasis,nsph),q(nbasis,nsph),rhs(nbasis,nsph), &
    & stat=istatus)
  if (istatus.ne.0) then
    write(*, "(A)") "ERROR: ddpcm_init allocation failed"
    stop
  end if
  call mkprec
  end subroutine ddpcm_init

  subroutine ddpcm_finalize()
  ! deallocate various arrays
  implicit none
  integer :: istatus
  deallocate(rx_prc,s,y,xs,phieps,q,rhs,stat=istatus)
  if (istatus.ne.0) then
    write(*, "(A)") "ERROR: ddpcm_finalize deallocation failed"
    stop
  end if
  end subroutine ddpcm_finalize

  subroutine ddpcm(do_adjoint, phi, psi, esolv)
  ! main ddpcm driver, given the potential at the exposed cavity
  ! points and the psi vector, computes the solvation energy
  implicit none
  real*8, intent(in) :: phi(ncav), psi(nbasis,nsph)
  logical, intent(in) :: do_adjoint
  real*8, intent(out) :: esolv
  integer :: istatus
  real*8 :: tol, fac
  integer :: isph, n_iter
  real*8, allocatable :: g(:,:), phiinf(:,:)
  logical :: ok
  external :: lx, ldm1x, lstarx, hnorm
  
  ! Allocated arrays are initialized later in corresponding routines
  allocate(phiinf(nbasis,nsph), g(ngrid,nsph), stat=istatus)
  if (istatus.ne.0) then
    write(*, "(A)") "ERROR: ddpcm allocation failed"
    stop
  end if
  tol = 10.0d0**(-iconv)

  ! build rhs
  ! `g` is initialized here
  call wghpot(phi,g)
  do isph = 1, nsph
  ! `rhs` is initialized here
    call intrhs(isph,g(:,isph),rhs(:,isph))
  end do

  ! rinf rhs
  ! `phiinf` is initialized here
  dodiag = .true.
  call rinfx(nbasis*nsph,rhs,phiinf)

  ! solve the ddpcm linear system
  n_iter = 200
  dodiag = .false.
  ! TODO: initial guess `phieps` must be updated properly
  phieps = rhs
  call cpu_time(start_time)
  call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,phiinf,phieps,n_iter, &
      & ok,rx,apply_rx_prec,hnorm)
  call cpu_time(finish_time)
  if (iprint.ge.1) then
      write(*, "(A,ES11.4E2,A)") " ddpcm step time:", &
        & finish_time-start_time, " seconds"
      write(*, "(A,I0)") " ddpcm step iterations: ", n_iter
  endif

  ! solve the ddcosmo linear system
  n_iter = 200
  dodiag = .false.
  call cpu_time(start_time)
  call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,phieps,xs,n_iter, &
      & ok,lx,ldm1x,hnorm)
  call cpu_time(finish_time)
  if (iprint.ge.1) then
      write(*, "(A,ES11.4E2,A)") " ddcosmo step time:", &
        & finish_time-start_time, " seconds"
      write(*, "(A,I0)") " ddcosmo step iterations: ", n_iter
  endif
  if (iprint.ge.2) call prtsph('x',nsph,0,xs)

  ! compute the energy
  esolv = pt5*sprod(nsph*nbasis,xs,psi)

  ! Solve adjoint system if needed
  if (do_adjoint) then

    ! solve ddcosmo adjoint system
    n_iter = 200
    dodiag = .false.
    call cpu_time(start_time)
    call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,psi,s,n_iter, &
      & ok,lstarx,ldm1x,hnorm)
    call cpu_time(finish_time)
    if (iprint.ge.1) then
        write(*, "(A,ES11.4E2,A)") " adjoint ddcosmo step time:", &
            & finish_time-start_time, " seconds"
        write(*, "(A,I0)") " adjoint ddcosmo step iterations: ", n_iter
    endif
    if (iprint.ge.2) call prtsph('s',nsph,0,s)

    ! solve ddpcm adjoint system
    n_iter = 200
    dodiag = .false.
    call cpu_time(start_time)
    call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,s,y,n_iter, &
      & ok,rstarx,apply_rstarx_prec,hnorm)
    call cpu_time(finish_time)
    if (iprint.ge.1) then
        write(*,"(A,ES11.4E2,A)") " adjoint ddpcm step time:", &
            & finish_time-start_time, " seconds"
        write(*,"(A,I0)") " adjoint ddpcm step iterations: ", &
            & n_iter
    endif
    if (iprint.ge.2) call prtsph('y',nsph,0,y)

    ! recover effect of rinf^*
    fac = two*pi*(one - (eps + one)/(eps - one))
    q = s + fac*y

    if (iprint.ge.2) call prtsph('q',nsph,0,q)
  end if
  deallocate(phiinf,g)
  end subroutine ddpcm

  subroutine ddpcm_fmm(do_adjoint, phi, psi, esolv, fmm_pm, fmm_pl, nadm)
  ! FMM-accelerated ddpcm driver, given the potential at the exposed cavity
  ! points and the psi vector, computes the solvation energy
  !
  ! Parameters:
  !     do_adjoint: flag whether to solve adjoint systems
  !     phi:
  !     psi:
  !     esolv: output solvation emergy
  !     pm: degree of multipole harmonics of the FMM
  !     pl: degree of local harmonics of the FMM
  !     nadm: upper bound on number of admissible far-field and near-field
  !         pairs per tree cluster. If this value is too small the procedure
  !         will raise and error and stop execution.
  implicit none
  real*8, intent(in) :: phi(ncav), psi(nbasis,nsph)
  real*8, intent(inout) :: esolv
  logical, intent(in) :: do_adjoint
  integer, intent(in) :: fmm_pm, fmm_pl ! ignored
  integer, intent(in) :: nadm
  integer :: istatus
  real*8 :: tol, resid, res0, rel_tol, fac
  integer :: isph, n_iter, l, ll, m
  real*8, allocatable :: g(:,:), phiinf(:,:)
  logical :: ok
  external :: lx, ldm1x, lstarx
  real*8, external :: hnorm, dnrm2
  ! Temporary variables to get list of all admissible pairs
  integer :: lwork, iwork, jwork
  integer, allocatable :: work(:, :)
  ! Sizes of FMM matrices
  integer :: m2m_reflect_size, m2m_ztrans_size
  integer :: l2l_reflect_size, l2l_ztrans_size
  integer :: m2l_reflect_size, m2l_ztrans_size
  ! Set globally degrees of FMM harmonics
  pm = pmax
  pl = pmax
  ! Prepare FMM tree and other things. This can be moved out from this
  ! function to rerun ddpcm_fmm with the same molecule without recomputing
  ! FMM-related variables.
  nclusters = 2*nsph - 1
  ! Magic constant that defines size of temporary buffer. Sometimes this size
  ! is not enough and therefore there will be an error message displayed,
  ! asking to increase this magic constant
  lwork = nclusters * nadm
  ! Allocate space for FMM constants
  allocate(vscales((pm+pl+1)*(pm+pl+1)), vgrid((pl+1)*(pl+1), ngrid), &
      & stat=istatus)
  if (istatus.ne.0) then
    write(*, "(A)") "ERROR: ddpcm_fmm allocation #1 failed"
    stop
  end if
  ! Init constants
  call init_globals(pm, pl, vscales, ngrid, w, grid, vgrid)
  ! Allocate cluster tree
  allocate(ind(nsph), cluster(2, nclusters), children(2, nclusters), &
      & parent(nclusters), cnode(3, nclusters), rnode(nclusters), &
      & snode(nsph), stat=istatus)
  if (istatus.ne.0) then
    write(*, "(A)") "ERROR: ddpcm_fmm allocation #2 failed"
    stop
  end if
  ! Init order of spheres
  do isph = 1, nsph
    ind(isph) = isph
  end do
  ! Build binary cluster tree
  call btree_init(nsph, csph, rsph, ind, cluster, children, parent, &
        cnode, rnode, snode)
  ! Allocate space for admissible far-field and near-field pairs
  ! More allocations will be done later regarding this step
  allocate(nfar(nclusters), nnear(nclusters), work(3, lwork), &
      & sfar(nclusters+1), snear(nclusters+1), stat=istatus)
  if (istatus.ne.0) then
    write(*, "(A)") "ERROR: ddpcm_fmm allocation #3 failed"
    stop
  end if
  ! Try to find all admissibly far and near pairs of tree nodes
  iwork = 0 ! init with zero for first call to tree_get_farnear_work
  call tree_get_farnear_work(nclusters, children, cnode, rnode, lwork, &
        iwork, jwork, work, nnfar, nfar, nnnear, nnear)
  ! Increase size of work array if needed and run again. Function
  ! tree_get_farnear_work uses previously computed work array, so it will not
  ! do the same work several times.
  if (iwork .ne. jwork+1) then
    write(*, "(A)") "Size of a temporary buffer for tree construction is &
        &too small, please increase parameter `nadm` of the ddpcm_fmm &
        &routine"
    stop
  end if
  ! Allocate arrays for admissible far and near pairs
  allocate(far(nnfar), near(nnnear), stat=istatus)
  if (istatus.ne.0) then
    write(*, "(A)") "ERROR: ddpcm_fmm allocation #4 failed"
    stop
  end if
  ! Get list of admissible pairs from temporary work array
  call tree_get_farnear(jwork, lwork, work, nclusters, nnfar, nfar, sfar, &
      far, nnnear, nnear, snear, near)
  ! Allocate compact external grid storage (the first portion)
  allocate(ngrid_ext_sph(nsph), stat=istatus)
  if (istatus.ne.0) then
    write(*, "(A)") "ERROR: ddpcm_fmm allocation #5 failed"
    stop
  end if
  ! Get number of external grid points
  call get_ngrid_ext(nsph, ngrid, ui, ngrid_ext, ngrid_ext_sph)
  ! Allocate compact external grid storage (the last portion)
  allocate(grid_ext_ia(nsph+1), grid_ext_ja(ngrid_ext), stat=istatus)
  if (istatus.ne.0) then
    write(*, "(A)") "ERROR: ddpcm_fmm allocation #6 failed"
    stop
  end if
  ! Get external grid points in a compact way
  call get_grid_ext_ind(nsph, ngrid, ui, ngrid_ext, ngrid_ext_sph, &
      & grid_ext_ia, grid_ext_ja)
  ! Get number of near-field M2P interactions
  call get_ngrid_ext_near(nsph, ngrid_ext_sph, nclusters, nnear, snode, &
      & ngrid_ext_near)
  ! Compute all transfer matrices only if needed
  if(ifmm_onfly .eq. 0) then
      ! Allocate FMM near-field M2P and far-field L2P matrices
      allocate(l2p_mat((pl+1)*(pl+1), ngrid_ext), &
          & m2p_mat((pm+1)*(pm+1), ngrid_ext_near), stat=istatus)
      if (istatus.ne.0) then
          write(*, "(A)") "ERROR: ddpcm_fmm allocation #7 failed"
          stop
      end if
      ! Allocate FMM far-field M2M matrices
      m2m_reflect_size = (pm+1) * (2*pm+1) * (2*pm+3) / 3
      m2m_ztrans_size = (pm+1) * (pm+2) * (pm+3) / 6
      allocate(m2m_reflect_mat(m2m_reflect_size, nclusters-1), &
          & m2m_ztrans_mat(m2m_ztrans_size, nclusters-1), stat=istatus)
      if (istatus.ne.0) then
          write(*, "(A)") "ERROR: ddpcm_fmm allocation #8 failed"
          stop
      end if
      ! Allocate FMM far-field L2L matrices
      l2l_reflect_size = (pl+1) * (2*pl+1) * (2*pl+3) / 3
      l2l_ztrans_size = (pl+1) * (pl+2) * (pl+3) / 6
      allocate(l2l_reflect_mat(l2l_reflect_size, nclusters-1), &
          & l2l_ztrans_mat(l2l_ztrans_size, nclusters-1), stat=istatus)
      if (istatus.ne.0) then
          write(*, "(A)") "ERROR: ddpcm_fmm allocation #9 failed"
          stop
      end if
      ! Allocate FMM far-field M2L matrices
      m2l_reflect_size = max(m2m_reflect_size, l2l_reflect_size)
      m2l_ztrans_size = (min(pm,pl)+1) * (min(pm,pl)+2) &
          & * (3*max(pm,pl)+3-min(pm,pl)) / 6
      allocate(m2l_reflect_mat(m2l_reflect_size, nnfar), &
          & m2l_ztrans_mat(m2l_ztrans_size, nnfar), stat=istatus)
      if (istatus.ne.0) then
          write(*, "(A)") "ERROR: ddpcm_fmm allocation #10 failed"
          stop
      end if
      ! Precompute all M2M, M2L and L2L reflection and OZ translation matrices
      call pcm_matvec_grid_fmm_get_mat(nsph, csph, rsph, ngrid, grid, w, vgrid, &
          & ui, pm, pl, vscales, ind, nclusters, cluster, snode, children, cnode, &
          & rnode, nnfar, sfar, far, nnnear, snear, near, m2m_reflect_mat, &
          & m2m_ztrans_mat, m2l_reflect_mat, m2l_ztrans_mat, l2l_reflect_mat, &
          & l2l_ztrans_mat, ngrid_ext, ngrid_ext_sph, grid_ext_ia, grid_ext_ja, &
          & l2p_mat, ngrid_ext_near, m2p_mat)
      write(*,*) "Allocated and computed FMM matrices"
  endif

  ! Continue with ddpcm
  allocate(phiinf(nbasis,nsph), g(ngrid,nsph), stat=istatus)
  if (istatus.ne.0) then
    write(*, "(A)") "ERROR: ddpcm_fmm allocation #11 failed"
    stop
  end if
  tol = 10.0d0**(-iconv)

  ! build rhs
  ! `g` is initialized here
  call wghpot(phi,g)
  do isph = 1, nsph
  ! `rhs` is initialized here
    call intrhs(isph,g(:,isph),rhs(:,isph))
  end do

  ! rinf rhs
  ! `phiinf` is initialized here
  dodiag = .true.
  call rinfx(nbasis*nsph,rhs,phiinf)

  ! solve the ddpcm linear system
  n_iter = 100000
  dodiag = .false.
  ! TODO: initial guess `phieps` must be updated properly
  phieps = rhs
  call cpu_time(start_time)
  call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,phiinf,phieps,n_iter, &
      & ok,rx_fmm,apply_rx_prec,hnorm)
  call cpu_time(finish_time)
  if (iprint.ge.1) then
      write(*, "(A,ES11.4E2,A)") " ddpcm step time:", &
        & finish_time-start_time, " seconds"
      write(*, "(A,I0)") " ddpcm step iterations: ", n_iter
  endif
  ! Print matvec statistics
  if (iprint.ge.2) call pcm_matvec_print_stats

  ! solve the ddcosmo linear system
  n_iter = 100000
  dodiag = .false.
  call cpu_time(start_time)
  call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,phieps,xs,n_iter, &
      & ok,lx,ldm1x,hnorm)
  call cpu_time(finish_time)
  if (iprint.ge.1) then
      write(*, "(A,ES11.4E2,A)") " ddcosmo step time:", &
        & finish_time-start_time, " seconds"
      write(*, "(A,I0)") " ddcosmo step iterations: ", n_iter
  endif
  if (iprint.ge.2) call prtsph('x',nsph,0,xs)

  ! compute the energy
  esolv = pt5*sprod(nsph*nbasis,xs,psi)

  ! Solve adjoint systems
  if (do_adjoint) then

    ! solve ddcosmo adjoint system
    n_iter = 200
    dodiag = .false.
    call cpu_time(start_time)
    call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,psi,s,n_iter, &
      & ok,lstarx,ldm1x,hnorm)
    call cpu_time(finish_time)
    if (iprint.ge.1) then
        write(*, "(A,ES11.4E2,A)") " adjoint ddcosmo step time:", &
            & finish_time-start_time, " seconds"
        write(*, "(A,I0)") " adjoint ddcosmo step iterations: ", n_iter
    endif
    if (iprint.ge.2) call prtsph('s',nsph,0,s)

    ! solve ddpcm adjoint system
    n_iter = 200
    dodiag = .false.
    call cpu_time(start_time)
    call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,s,y,n_iter, &
      & ok,rstarx_fmm,apply_rstarx_prec,hnorm)
    call cpu_time(finish_time)
    if (iprint.ge.1) then
        write(*,"(A,ES11.4E2,A)") " adjoint ddpcm step time:", &
            & finish_time-start_time, " seconds"
        write(*,"(A,I0)") " adjoint ddpcm step iterations: ", &
            & n_iter
    endif
    if (iprint.ge.2) call prtsph('y',nsph,0,y)

    ! recover effect of Rinf^*
    fac = two*pi*(one - (eps + one)/(eps - one))
    q = s + fac*y

    if (iprint.ge.2) call prtsph('q',nsph,0,q)
  end if

  ! DO NOT DEALLOCATE HERE, AS IT IS USED IN MAIN_FMM.f90
  ! Deallocate ddpcm-related arrays
  !deallocate(phiinf,g)
  ! Deallocate FMM-related arrays
  !deallocate(ind, snode)
  !deallocate(cluster, children, parent, cnode, rnode)
  !deallocate(vscales, vgrid)
  !deallocate(nfar, nnear)
  !deallocate(work)
  !deallocate(far, sfar, near, snear)
  !deallocate(ngrid_ext_sph)
  !deallocate(grid_ext_ia, grid_ext_ja)
  ! Deallocate transfer matrices if they were allocated/computed
  !if(ifmm_onfly .eq. 2) then
  !    deallocate(l2p_mat)
  !    deallocate(m2p_mat)
  !    deallocate(m2m_reflect_mat)
  !    deallocate(m2m_ztrans_mat)
  !    deallocate(l2l_reflect_mat)
  !    deallocate(l2l_ztrans_mat)
  !    deallocate(m2l_reflect_mat)
  !    deallocate(m2l_ztrans_mat)
  !endif
  end subroutine ddpcm_fmm

  subroutine ddpcm_fmm_adj(phi, psi, esolv)
  ! FMM-accelerated ddpcm driver, given the potential at the exposed cavity
  ! points and the psi vector, computes the solvation energy
  implicit none
  real*8, intent(in) :: phi(ncav), psi(nbasis,nsph)
  real*8, intent(inout) :: esolv
  real*8 :: tol, resid, res0, rel_tol, fac
  integer :: isph, n_iter, l, ll, m
  logical :: ok
  external :: lx, ldm1x
  real*8, external :: hnorm, dnrm2
  real(kind=8), allocatable :: fmm_mat(:, :, :, :)
  real(kind=8), allocatable :: coef_sph(:, :), coef_out(:, :)
  real(kind=8) :: total_norm, total_err
  integer :: i, j
  ! Temporary variables to get list of all admissible pairs
  integer :: lwork, iwork, jwork
  integer, allocatable :: work(:, :)
  ! Prepare FMM tree and other things. This can be moved out from this
  ! function to rerun ddpcm_fmm with the same molecule without recomputing
  ! FMM-related variables.
  nclusters = 2*nsph-1
  allocate(ind(nsph), cluster(2, nclusters), children(2, nclusters), &
      parent(nclusters), cnode(3, nclusters), rnode(nclusters), snode(nsph))
  pm = lmax ! This must be updated to achieve required accuracy of FMM
  pl = lmax ! This must be updated to achieve required accuracy of FMM
  allocate(vscales((pm+pl+1)*(pm+pl+1)), vgrid((pl+1)*(pl+1), ngrid))
  ! Init order of spheres
  do isph = 1, nsph
    ind(isph) = isph
  end do
  ! Init constants vscales
  call init_globals(pm, pl, vscales, ngrid, w, grid, vgrid)
  ! Build binary tree
  call btree_init(nsph, csph, rsph, ind, cluster, children, parent, &
        cnode, rnode, snode)
  nclusters = 2*nsph-1
  ! Allocate space to find all admissible pairs
  allocate(nfar(nclusters), nnear(nclusters))
  ! Try to find all admissibly far and near pairs of tree nodes
  lwork = nclusters*1000 ! Some magic constant which might be changed
  allocate(work(3, lwork))
  iwork = 0 ! init with zero for first call to tree_get_farnear_work
  call tree_get_farnear_work(nclusters, children, cnode, rnode, lwork, &
        iwork, jwork, work, nnfar, nfar, nnnear, nnear)
  ! Increase size of work array if needed and run again. Function
  ! tree_get_farnear_work uses previously computed work array, so it will not
  ! do the same work several times.
  if (iwork .ne. jwork+1) then
    write(*,*) 'Please increase lwork on line 162 of ddpcm_lib.f90 in the code'
    stop
  end if
  ! Allocate arrays for admissible far and near pairs
  allocate(far(nnfar), near(nnnear), sfar(nclusters+1), snear(nclusters+1))
  ! Get list of admissible pairs from temporary work array
  call tree_get_farnear(jwork, lwork, work, nclusters, nnfar, nfar, sfar, &
      far, nnnear, nnear, snear, near)
  ! Get external grid points
  allocate(ngrid_ext_sph(nsph))
  call get_ngrid_ext(nsph, ngrid, ui, ngrid_ext, ngrid_ext_sph)
  allocate(grid_ext_ia(nsph+1), grid_ext_ja(ngrid_ext))
  call get_grid_ext_ind(nsph, ngrid, ui, ngrid_ext, ngrid_ext_sph, &
      & grid_ext_ia, grid_ext_ja)
  ! Allocate far-field L2P matrices
  allocate(l2p_mat((pl+1)*(pl+1), ngrid_ext))
  ! Get near-field M2P data
  call get_ngrid_ext_near(nsph, ngrid_ext_sph, nclusters, nnear, snode, &
      & ngrid_ext_near)
  ! Allocate near-field M2P matrices
  allocate(m2p_mat((pm+1)*(pm+1), ngrid_ext_near))
  ! Allocate reflection and OZ translation matrices for M2M, M2L and L2L
  allocate(m2m_reflect_mat((pm+1)*(2*pm+1)*(2*pm+3)/3, nclusters-1))
  allocate(m2m_ztrans_mat((pm+1)*(pm+2)*(pm+3)/6, nclusters-1))
  allocate(l2l_reflect_mat((pl+1)*(2*pl+1)*(2*pl+3)/3, nclusters-1))
  allocate(l2l_ztrans_mat((pl+1)*(pl+2)*(pl+3)/6, nclusters-1))
  allocate(m2l_reflect_mat((max(pm,pl)+1)*(2*max(pm,pl)+1) &
      & *(2*max(pm,pl)+3)/3, nnfar))
  allocate(m2l_ztrans_mat((min(pm,pl)+1)*(min(pm,pl)+2) &
      & *(3*max(pm,pl)+3-min(pm,pl))/6, nnfar))
  ! Precompute all M2M, M2L and L2L reflection and OZ translation matrices
  call pcm_matvec_grid_fmm_get_mat(nsph, csph, rsph, ngrid, grid, w, vgrid, &
      & ui, pm, pl, vscales, ind, nclusters, cluster, snode, children, cnode, &
      & rnode, nnfar, sfar, far, nnnear, snear, near, m2m_reflect_mat, &
      & m2m_ztrans_mat, m2l_reflect_mat, m2l_ztrans_mat, l2l_reflect_mat, &
      & l2l_ztrans_mat, ngrid_ext, ngrid_ext_sph, grid_ext_ia, grid_ext_ja, &
      & l2p_mat, ngrid_ext_near, m2p_mat)


  ! Allocate memory for the matrix
  allocate(fmm_mat((pl+1)*(pl+1), nsph, (pm+1)*(pm+1), nsph))
  fmm_mat = 0
  call tree_matrix_fmm(nsph, csph, rsph, ngrid, grid, w, &
        & vgrid, ui, pm, pl, vscales, ind, nclusters, cluster, snode, &
        & children, &
        & cnode, rnode, nnfar, sfar, far, nnnear, snear, near, &
        & m2m_reflect_mat, m2m_ztrans_mat, m2l_reflect_mat, m2l_ztrans_mat, &
        & l2l_reflect_mat, l2l_ztrans_mat, ngrid_ext, ngrid_ext_sph, &
        & grid_ext_ia, grid_ext_ja, l2p_mat, ngrid_ext_near, m2p_mat, &
        & fmm_mat)
  write(*,*) "FMM matrix generated"
  total_norm = dnrm2((pm+1)*(pm+1)*nsph*(pl+1)*(pl+1)*nsph, fmm_mat, 1)

  ! Perform adjoint matvec
  allocate(coef_sph((pl+1)*(pl+1), nsph), coef_out((pm+1)*(pm+1), nsph))
  do i = 1, (pl+1)*(pl+1)
    do j = 1, nsph
      coef_sph = 0
      coef_sph(i, j) = 1
      coef_out = 0
      call pcm_matvec_adj_grid_fmm_use_mat(nsph, csph, rsph, ngrid, grid, w, &
            & vgrid, ui, pm, pl, vscales, ind, nclusters, cluster, snode, &
            & children, &
            & cnode, rnode, nnfar, sfar, far, nnnear, snear, near, &
            & m2m_reflect_mat, m2m_ztrans_mat, m2l_reflect_mat, m2l_ztrans_mat, &
            & l2l_reflect_mat, l2l_ztrans_mat, ngrid_ext, ngrid_ext_sph, &
            & grid_ext_ia, grid_ext_ja, l2p_mat, ngrid_ext_near, m2p_mat, &
            & coef_sph, coef_out)
      fmm_mat(i, j, :, :) = fmm_mat(i, j, :, :) - coef_out
    end do
  end do
  total_err = dnrm2((pm+1)*(pm+1)*nsph*(pl+1)*(pl+1)*nsph, fmm_mat, 1)
  write(*,*) "Error in adjoint matrix: ", total_err, total_norm
  ! Deallocate matrix
  deallocate(fmm_mat)
  deallocate(coef_sph, coef_out)
  ! Deallocate FMM-related arrays
  deallocate(ind, snode)
  deallocate(cluster, children, parent, cnode, rnode)
  deallocate(vscales, vgrid)
  deallocate(nfar, nnear)
  deallocate(work)
  deallocate(far, sfar, near, snear)
  deallocate(ngrid_ext_sph)
  deallocate(grid_ext_ia, grid_ext_ja)
  deallocate(l2p_mat)
  deallocate(m2p_mat)
  deallocate(m2m_reflect_mat)
  deallocate(m2m_ztrans_mat)
  deallocate(l2l_reflect_mat)
  deallocate(l2l_ztrans_mat)
  deallocate(m2l_reflect_mat)
  deallocate(m2l_ztrans_mat)
  return
  end subroutine ddpcm_fmm_adj

  subroutine mkprec
  ! assemble the diagonal blocks of the reps matrix
  ! then invert them to build the preconditioner
  implicit none
  integer :: isph, lm, ind, l1, m1, ind1, its, istatus
  real*8  :: f, f1
  integer, allocatable :: ipiv(:)
  real*8,  allocatable :: work(:)

  allocate(ipiv(nbasis),work(nbasis*nbasis),stat=istatus)
  if (istatus.ne.0) then
    write(*,*) 'mkprec : allocation failed !'
    stop
  endif

  rx_prc(:,:,:) = zero

  ! off diagonal part
  do isph = 1, nsph
    do its = 1, ngrid
      f = two*pi*ui(its,isph)*w(its)
      do l1 = 0, lmax
        ind1 = l1*l1 + l1 + 1
        do m1 = -l1, l1
          f1 = f*basis(ind1 + m1,its)/(two*dble(l1) + one)
          do lm = 1, nbasis
            rx_prc(lm,ind1 + m1,isph) = rx_prc(lm,ind1 + m1,isph) + &
                & f1*basis(lm,its)
          end do
        end do
      end do
    end do
  end do

  ! add diagonal
  f = two*pi*(eps + one)/(eps - one)
  do isph = 1, nsph
    do lm = 1, nbasis
      rx_prc(lm,lm,isph) = rx_prc(lm,lm,isph) + f
    end do
  end do

  ! invert the blocks
  do isph = 1, nsph
    call dgetrf(nbasis,nbasis,rx_prc(:,:,isph),nbasis,ipiv,istatus)
    if (istatus.ne.0) then 
      write(6,*) 'lu failed in mkprc'
      stop
    end if
    call dgetri(nbasis,rx_prc(:,:,isph),nbasis,ipiv,work, &
        & nbasis*nbasis,istatus)
    if (istatus.ne.0) then 
      write(6,*) 'inversion failed in mkprc'
      stop
    end if
  end do

  deallocate (work,ipiv,stat=istatus)
  if (istatus.ne.0) then
    write(*,*) 'mkprec : deallocation failed !'
    stop
  end if
  return
  endsubroutine mkprec


  subroutine rx(n,x,y)
  ! computes y = reps x =
  ! = (2*pi*(eps + 1)/(eps - 1) - d) x 
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(out) :: y(nbasis,nsph)
  real*8 :: fac

  ! output `y` is cleaned here
  call dx(n,x,y)
  y = -y

  if (dodiag) then
    fac = two*pi*(eps + one)/(eps - one)
    y = y + fac*x
  end if
  end subroutine rx

  subroutine rx_fmm(n, x, y)
  ! Use FMM to compute Y = Reps X =
  ! = (2*pi*(eps + 1)/(eps - 1) - D) X 
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(out) :: y(nbasis,nsph)
  real*8 :: fac
  integer :: isph

  call dx_fmm(n,x,y)
  y = -y

  if (dodiag) then
    fac = two*pi*(eps + one)/(eps - one)
    y = y + fac*x
  end if
  do isph = 1, nsph
    !y(:, isph) = y(:, isph) * rsph(isph)**2
  end do
  end subroutine rx_fmm

  subroutine rinfx(n,x,y)
  ! computes y = rinf x = 
  ! = (2*pi -  d) x
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(out) :: y(nbasis,nsph)
  real*8 :: fac

  ! output `y` is cleaned here
  call dx(n,x,y)
  y = -y

  if (dodiag) then
    fac = two*pi
    y = y + fac*x
  end if
  end subroutine rinfx

  subroutine rinfx_fmm(n,x,y)
  ! Computes Y = Rinf X = 
  ! = (2*pi -  D) X
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(out) :: y(nbasis,nsph)
  real*8 :: fac
  integer :: isph

  call dx_fmm(n,x,y)
  y = -y

  if (dodiag) then
    fac = two*pi
    y = y + fac*x
  end if
  end subroutine rinfx_fmm

  subroutine dx(n,x,y)
  ! computes y = d x
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(out) :: y(nbasis,nsph)
  real*8, allocatable :: vts(:), vplm(:), basloc(:), vcos(:), vsin(:)
  real*8 :: c(3), vij(3), sij(3)
  real*8 :: vvij, tij, fourpi, tt, f, f1
  integer :: its, isph, jsph, l, m, ind, lm, istatus
  
  allocate(vts(ngrid),vplm(nbasis),basloc(nbasis),vcos(lmax + 1), &
      & vsin(lmax + 1),stat=istatus)
  if (istatus.ne.0) then
    write(6,*) 'dx: allocation failed !'
    stop
  end if
  y = zero
  fourpi = four*pi
  ! this loop is easily parallelizable
  ! !!$omp parallel do default(none) schedule(dynamic) &
  ! !!$omp private(isph,its,jsph,basloc,vplm,vcos,vsin,vij, &
  ! !!$omp vvij,tij,sij,tt,l,ind,f,m,vts,c) &
  ! !!$omp shared(nsph,ngrid,ui,csph,rsph,grid, &
  ! !!$omp lmax,fourpi,dodiag,x,y,basis)
  do isph = 1, nsph
    ! compute the "potential" from the other spheres
    ! at the exposed lebedv points of the i-th sphere 
    vts = zero
    do its = 1, ngrid
      if (ui(its,isph).gt.zero) then
        c = csph(:,isph) + rsph(isph)*grid(:,its)
        do jsph = 1, nsph
          if (jsph.ne.isph) then
            ! build the geometrical variables
            vij = c - csph(:,jsph)
            vvij = sqrt(dot_product(vij,vij))
            tij = vvij/rsph(jsph)
            sij = vij/vvij 
            ! build the local basis
            call ylmbas(sij,basloc,vplm,vcos,vsin)
            ! with all the required stuff, finally compute
            ! the "potential" at the point 
            tt = one/tij 
            do l = 0, lmax
              ind = l*l + l + 1
              f = fourpi*dble(l)/(two*dble(l) + one)*tt
              do m = -l, l
                vts(its) = vts(its) + f*x(ind + m,jsph) * &
                    & basloc(ind + m)
              end do
              tt = tt/tij
            end do
          else if (dodiag) then
            ! add the diagonal contribution
            do l = 0, lmax
              ind = l*l + l + 1
              f = (two*dble(l) + one)/fourpi
              do m = -l, l
                vts(its) = vts(its) - pt5*x(ind + m,isph) * &
                    & basis(ind + m,its)/f
              end do
            end do
          end if 
        end do
        vts(its) = ui(its,isph)*vts(its) 
      end if
    end do
    ! now integrate the potential to get its modal representation
    call intrhs(isph,vts,y(:,isph))
  end do

  deallocate(vts,vplm,basloc,vcos,vsin,stat=istatus)
  if (istatus.ne.0) then
    write(6,*) 'dx: deallocation failed !'
    stop
  end if
  end subroutine dx

  
  subroutine dx_fmm(n,x,y)
  ! Computes Y = D X with help of FMM
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(out) :: y(nbasis,nsph)
  real*8 :: fmm_x((pm+1)*(pm+1), nsph) ! Input for FMM
  real*8 :: fmm_y((pl+1)*(pl+1), nsph) ! Output of FMM
  integer :: isph, igrid, l, lind, m
  real*8 :: f, vts(ngrid), fourpi
  if(nbasis .le. (pm+1)*(pm+1)) then
    fmm_x(1:nbasis, :) = x
    fmm_x(nbasis+1:, :) = zero
  else
    fmm_x(:, :) = x(1:(pm+1)*(pm+1), :)
  endif
  ! Apply diagonal contribution if needed
  if(dodiag) then
    fourpi = four * pi
    do isph = 1, nsph
      do igrid = 1, ngrid
        vts(igrid) = zero
        do l = 0, lmax
          lind = l*l + l + 1
          f = (two*dble(l)+one) / fourpi
          do m = -l, l
            vts(igrid) = vts(igrid) - pt5*x(lind+m,isph)*basis(lind+m,igrid)/f
          end do
        end do
        vts(igrid) = ui(igrid, isph) * vts(igrid)
      end do
      call intrhs(isph, vts, y(:, isph))
    end do
  else
    y = 0
  end if
  ! Do actual FMM matvec
  if(ifmm_onfly .eq. 1) then
      call pcm_matvec_grid_fmm_fast(nsph, csph, rsph, ngrid, grid, w, vgrid, ui, &
          pm, pl, vscales, ind, nclusters, cluster, children, cnode, rnode, &
          nnfar, sfar, far, nnnear, snear, near, fmm_x, fmm_y)
  else
      call pcm_matvec_grid_fmm_use_mat(nsph, csph, rsph, ngrid, grid, w, vgrid, &
          & ui, pm, pl, vscales, ind, nclusters, cluster, snode, children, cnode, &
          & rnode, nnfar, sfar, far, nnnear, snear, near, m2m_reflect_mat, &
          & m2m_ztrans_mat, m2l_reflect_mat, m2l_ztrans_mat, l2l_reflect_mat, &
          & l2l_ztrans_mat, ngrid_ext, ngrid_ext_sph, grid_ext_ia, grid_ext_ja, &
          & l2p_mat, ngrid_ext_near, m2p_mat, fmm_x, fmm_y)
  endif
  !write(*,*) 'I did actual FMM matvec'
  ! Cut fmm_y to fit into output y
  if(nbasis .le. (pl+1)*(pl+1)) then
      y = y + fmm_y(1:nbasis, :)
  else
      y(1:(pl+1)*(pl+1), :) = y(1:(pl+1)*(pl+1), :) + fmm_y
  endif
  end subroutine dx_fmm

  subroutine apply_rx_prec(n,x,y)
  ! apply the block diagonal preconditioner
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  integer :: isph
  ! simply do a matrix-vector product with the stored preconditioner 
  ! !!$omp parallel do default(shared) schedule(dynamic) &
  ! !!$omp private(isph)
  do isph = 1, nsph
    call dgemm('n','n',nbasis,1,nbasis,one,rx_prc(:,:,isph),nbasis, &
      & x(:,isph),nbasis,zero,y(:,isph),nbasis)
  end do
  end subroutine apply_rx_prec

  subroutine rstarx(n,x,y)
  ! computes y = reps x =
  ! = (2*pi*(eps + 1)/(eps - 1) - d^*) x 
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  real*8 :: fac

  call dstarx(n,x,y)
  y = -y

  if (dodiag) then
    fac = two*pi*(eps + one)/(eps - one)
    y = y + fac*x
  end if
  end subroutine rstarx

  subroutine rstarx_fmm(n,x,y)
  ! Computes Y = Reps X =
  ! = (2*pi*(eps + 1)/(eps - 1) - D^*) X 
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  real*8 :: fac

  call dstarx_fmm(n,x,y)
  y = -y

  if (dodiag) then
    fac = two*pi*(eps + one)/(eps - one)
    y = y + fac*x
  end if
  end subroutine rstarx_fmm

  subroutine dstarx(n,x,y)
  ! computes y = d^* x
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  real*8, allocatable :: vts(:), vplm(:), basloc(:), vcos(:), vsin(:)
  real*8 :: c(3), vji(3), sji(3)
  real*8 :: vvji, tji, fourpi, tt, f, f1
  integer :: its, isph, jsph, l, m, ind, lm, istatus
  
  allocate(vts(ngrid),vplm(nbasis),basloc(nbasis),vcos(lmax + 1), &
      & vsin(lmax + 1),stat=istatus)
  if (istatus.ne.0) then
    write(6,*) 'dx: allocation failed !'
    stop
  end if
  y = zero
  fourpi = four*pi

  !$omp parallel do default(none) schedule(dynamic) &
  !$omp private(isph,its,jsph,basloc,vplm,vcos,vsin,vji, &
  !$omp vvji,tji,sji,tt,l,ind,f,m,vts,c) &
  !$omp shared(nsph,ngrid,ui,csph,rsph,grid,facl, &
  !$omp lmax,fourpi,dodiag,x,y,basis,w,nbasis)
  do isph = 1, nsph
    do jsph = 1, nsph
      if (jsph.ne.isph) then
        do its = 1, ngrid
          if (ui(its,jsph).gt.zero) then
            ! build the geometrical variables
            vji = csph(:,jsph) + rsph(jsph)*grid(:,its) - csph(:,isph)
            vvji = sqrt(dot_product(vji,vji))
            tji = vvji/rsph(isph)
            sji = vji/vvji
            ! build the local basis
            call ylmbas(sji,basloc,vplm,vcos,vsin)
            tt = w(its)*ui(its,jsph)*dot_product(basis(:,its),x(:,jsph))/tji
            do l = 0, lmax
              ind = l*l + l + 1
              f = dble(l)*tt/facl(ind)
              do m = -l, l
                y(ind+m,isph) = y(ind+m,isph) + f*basloc(ind+m)
              end do
              tt = tt/tji
            end do
          end if
        end do
      else if (dodiag) then
        do its = 1, ngrid
          f = pt5*w(its)*ui(its,jsph)*dot_product(basis(:,its),x(:,jsph))
          do ind = 1, nbasis
            y(ind,isph) = y(ind,isph) - f*basis(ind,its)/facl(ind)
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

  subroutine dstarx_fmm(n,x,y)
  ! Computes Y = D^* X with help of FMM
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  real*8 :: fmm_x((pl+1)*(pl+1), nsph) ! Input for adjoint FMM
  real*8 :: fmm_y((pm+1)*(pm+1), nsph) ! Output of adjoint FMM
  integer :: isph, igrid, lind
  real*8 :: f, vts(ngrid), fourpi
  fmm_x(1:nbasis, :) = x
  fmm_x(nbasis+1:, :) = zero
  ! Do actual adjoint FMM matvec
  call pcm_matvec_adj_grid_fmm_use_mat(nsph, csph, rsph, ngrid, grid, w, vgrid, &
      & ui, pm, pl, vscales, ind, nclusters, cluster, snode, children, cnode, &
      & rnode, nnfar, sfar, far, nnnear, snear, near, m2m_reflect_mat, &
      & m2m_ztrans_mat, m2l_reflect_mat, m2l_ztrans_mat, l2l_reflect_mat, &
      & l2l_ztrans_mat, ngrid_ext, ngrid_ext_sph, grid_ext_ia, grid_ext_ja, &
      & l2p_mat, ngrid_ext_near, m2p_mat, fmm_x, fmm_y)
  ! Apply diagonal contribution if needed
  if(dodiag) then
    fourpi = four * pi
    do isph = 1, nsph
      do igrid = 1, ngrid
        f = pt5*w(igrid)*ui(igrid,isph)*dot_product(basis(:,igrid),x(:,isph))
        do lind = 1, nbasis
          y(lind,isph) = y(lind,isph) - f*basis(lind,igrid)/facl(lind)
        end do
      end do
    end do
    ! Cut fmm_y to fit into output y
    y = y + fmm_y(1:nbasis, :)
  else
    ! Cut fmm_y to fit into output y
    y = fmm_y(1:nbasis, :)
  end if
  !write(*,*) 'I did actual FMM matvec'
  end subroutine dstarx_fmm

  subroutine apply_rstarx_prec(n,x,y)
  ! apply the block diagonal preconditioner
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  integer :: isph
  ! simply do a matrix-vector product with the transposed preconditioner 
  !$omp parallel do default(shared) schedule(dynamic) &
  !$omp private(isph)
  do isph = 1, nsph
    call dgemm('t','n',nbasis,1,nbasis,one,rx_prc(:,:,isph),nbasis, &
      & x(:,isph),nbasis,zero,y(:,isph),nbasis)
  end do
  end subroutine apply_rstarx_prec

  subroutine ddpcm_forces(phi,fx)
  ! compute the geometrical contribution to the forces
  implicit none
  real*8, intent(in) :: phi(ncav)
  real*8, intent(inout) :: fx(3,nsph)
  real*8, allocatable :: vsin(:), vcos(:), vplm(:), basloc(:), &
    & dbsloc(:,:), phiexp(:,:), phieexp(:,:), f(:,:)
  ! cartesian representation of various adjoint solutions
  real*8, allocatable :: scr(:,:), ycr(:,:), qcr(:,:)
  integer :: istatus, ii, isph, ig
  ! debug
  real*8 :: fx_test(3,nsph)
  real*8, allocatable :: fscr(:,:)
  real*8 :: fac
  real*8, external :: dnrm2
  allocate(fscr(3,nsph))
  allocate(vsin(lmax+1),vcos(lmax+1),vplm(nbasis),basloc(nbasis), &
    & dbsloc(3,nbasis),scr(ngrid,nsph), &
    & phiexp(ngrid,nsph),ycr(ngrid,nsph),qcr(ngrid,nsph),stat=istatus)
  if (istatus.ne.0) write(6,*) 'ddpcm forces allocation failed'

  fx = zero

  ! expand the adjoints
  !$omp parallel do default(shared) private(isph,ig)
  do isph = 1, nsph
    do ig = 1, ngrid
      scr(ig,isph) = dot_product(s(:,isph),basis(:,ig))
      ycr(ig,isph) = dot_product(y(:,isph),basis(:,ig))
    end do
  end do
  !$omp end parallel do
  fac = two*pi*(one - (eps + one)/(eps - one))
  qcr = scr + fac*ycr

  ! expand the potential on a sphere-by-sphere basis (needed for parallelism):
  ii = 0
  phiexp = zero
  do isph = 1, nsph
    do ig = 1, ngrid
      if (ui(ig,isph).gt.zero) then
        ii = ii + 1
        phiexp(ig,isph) = phi(ii)
      end if
    end do
  end do

  ! compute the geometrical contributions from the ddcosmo matrix
  ! (fdoka, fdokb), from the ddpcm matrix (gradr) and from the
  ! geometrical part of the rhs (fdoga) 
  fx = zero
  fscr = zero
  write(6,*) 'fdoka'
  do isph = 1, nsph
    call fdoka(isph,xs,scr(:,isph),basloc,dbsloc,vplm,vcos,vsin,fscr(:,isph)) 
    write(6,'(1x,i5,3f16.8)') isph, -pt5*fscr(:,isph)
  end do
  fx = fx + fscr
  fscr = zero
  write(6,*) 'fdokb'
  do isph = 1, nsph
    call fdokb(isph,xs,scr,basloc,dbsloc,vplm,vcos,vsin,fscr(:,isph)) 
    write(6,'(1x,i5,3f16.8)') isph, -pt5*fscr(:,isph)
  end do
  fx = fx + fscr
  fscr = zero
  write(6,*) 'gradr'
  fx_test = 0
  call pcm_force_grid_fmm_use_mat(nsph, csph, rsph, ngrid, grid, w, &
        & vgrid, ui, pm, pl, lmax, vscales, ind, nclusters, cluster, snode, &
        & children, &
        & cnode, rnode, nnfar, sfar, far, nnnear, snear, near, &
        & m2m_reflect_mat, m2m_ztrans_mat, m2l_reflect_mat, m2l_ztrans_mat, &
        & l2l_reflect_mat, l2l_ztrans_mat, ngrid_ext, ngrid_ext_sph, &
        & grid_ext_ia, grid_ext_ja, l2p_mat, ngrid_ext_near, m2p_mat, &
        & rhs-phieps, ycr, fx_test)
  do isph = 1, nsph
    call gradr(isph,vplm,vcos,vsin,basloc,dbsloc,rhs-phieps,ycr,fscr(:,isph))
    write(6,'(1x,i5,3e16.8,3e16.8)') isph, -pt5*fscr(:,isph), -pt5*fx_test(:,isph)
  end do
  write(*,*) "Rel.error of gradr: ", dnrm2(3*nsph, fx_test-fscr, 1) / &
      & dnrm2(3*nsph, fscr, 1)
  fx = fx + fscr
  !fx = fx + fx_test
  fscr = zero
  write(6,*) 'fdoga'
  do isph = 1, nsph
    call fdoga(isph,qcr,phiexp,fscr(:,isph)) 
    write(6,'(1x,i5,3f16.8)') isph, -pt5*fscr(:,isph)
  end do
  fx = fx + fscr

  fx = -pt5*fx

  write(6,*) 'geometrical forces'
  do isph = 1, nsph
    write(6,'(1x,i5,3e16.8)') isph, fx(:,isph)
  end do

  deallocate(vsin,vcos,vplm,basloc,dbsloc,scr,phiexp,ycr,qcr,stat=istatus)
  if (istatus.ne.0) write(6,*) 'ddpcm forces deallocation failed'
  ! debug
  deallocate(fscr)
  end subroutine ddpcm_forces

  subroutine ddpcm_zeta(zeta)
  ! returns adjoint solution expansion at the grid points
  ! needed for analytical derivatives
  implicit none
  real*8, intent(out) :: zeta(ncav)
  integer :: ii, isph, its
  ii = 0
  do isph = 1, nsph
    do its = 1, ngrid
      if (ui(its,isph) .gt. zero) then
        ii = ii + 1
        zeta(ii) = w(its)*ui(its,isph)*dot_product(basis(:,its),q(:,isph))
      end if
    end do
  end do
  zeta = -pt5*zeta 
  end subroutine ddpcm_zeta

  subroutine gradr(isph,vplm,vcos,vsin,basloc,dbsloc,g,y,fx)
  use ddcosmo
  implicit none
  ! compute the gradient of ddpcm r and contract it
  ! < y, grad r (phie - phi) >

  ! physical quantities
  real*8 :: g(nbasis,*), y(ngrid,*), fx(*)
  ! various scratch arrays
  real*8 vik(3), sik(3), vki(3), ski(3), vkj(3), skj(3), vji(3), &
    & sji(3), va(3), vb(3), a(3), d(3)
  ! jacobian matrix
  real*8 sjac(3,3)
  ! other scratch arrays
  real*8 vplm(*), vcos(*), vsin(*), basloc(*), dbsloc(3,*)
  ! indexes
  integer isph, its, ik, ksph, l, m, ind, jsph, icomp, jcomp
  ! various scalar quantities
  real*8 cx, cy, cz, vvki, tki, gg, fl, fac, vvkj, tkj
  real*8 tt, fcl, dij, fjj, gi, fii, vvji, tji, qji
  real*8 b, vvik, tik, qik, tlow, thigh, duj

  tlow  = one - pt5*(one - se)*eta
  thigh = one + pt5*(one + se)*eta

  ! first set of contributions:
  ! diagonal block, kc and part of kb

  do its = 1, ngrid
    ! sum over ksph in neighbors of isph
    do ik = inl(isph), inl(isph+1) - 1
      ksph = nl(ik)
      ! build geometrical quantities
      cx = csph(1,ksph) + rsph(ksph)*grid(1,its)
      cy = csph(2,ksph) + rsph(ksph)*grid(2,its)
      cz = csph(3,ksph) + rsph(ksph)*grid(3,its)
      vki(1) = cx - csph(1,isph)
      vki(2) = cy - csph(2,isph)
      vki(3) = cz - csph(3,isph)
      vvki = sqrt(vki(1)*vki(1) + vki(2)*vki(2) + &
        & vki(3)*vki(3))
      tki  = vvki/rsph(isph)

      ! contributions involving grad i of uk come from the switching
      ! region.
      ! note: ui avoids contributions from points that are in the
      ! switching between isph and ksph but are buried in a third
      ! sphere.
      if ((tki.gt.tlow).and.(tki.lt.thigh) .and. &
        & ui(its,ksph).gt.zero) then
        ! other geometrical quantities
        ski = vki/vvki

        ! diagonal block kk contribution, with k in n(i)
        gg = zero
        do l = 0, lmax
          ind = l*l + l + 1
          fl = dble(l)
          fac = two*pi/(two*fl + one)
          do m = -l, l 
            !! DEBUG comment
            gg = gg + fac*basis(ind+m,its)*g(ind+m,ksph)
          end do
        end do

        ! kc contribution
        do jsph = 1, nsph
          if (jsph.ne.ksph .and. jsph.ne.isph) then 
            vkj(1) = cx - csph(1,jsph)
            vkj(2) = cy - csph(2,jsph)
            vkj(3) = cz - csph(3,jsph)
            vvkj = sqrt(vkj(1)*vkj(1) + vkj(2)*vkj(2) + &
              & vkj(3)*vkj(3))
            tkj  = vvkj/rsph(jsph)
            skj  = vkj/vvkj
            call ylmbas(skj,basloc,vplm,vcos,vsin)
            tt = one/tkj
            do l = 0, lmax
              ind = l*l + l + 1
              fcl = - four*pi*dble(l)/(two*dble(l)+one)*tt
              do m = -l, l
                !! DEBUG comment
                gg = gg + fcl*g(ind+m,jsph)*basloc(ind+m)
              end do
              tt = tt/tkj
            end do
          end if
        end do

        ! part of kb contribution
        call ylmbas(ski,basloc,vplm,vcos,vsin)
        tt = one/tki
        do l = 0, lmax
          ind = l*l + l + 1
          fcl = - four*pi*dble(l)/(two*dble(l)+one)*tt
          do m = -l, l
            !! DEBUG comment
            gg = gg + fcl*g(ind+m,isph)*basloc(ind+m)
          end do
          tt = tt/tki
        end do

        ! common step, product with grad i uj
        duj = dfsw(tki,se,eta)/rsph(isph)
        fjj = duj*w(its)*gg*y(its,ksph)
        fx(1) = fx(1) - fjj*ski(1)
        fx(2) = fx(2) - fjj*ski(2)
        fx(3) = fx(3) - fjj*ski(3)
      end if
    end do

    ! diagonal block ii contribution
    if (ui(its,isph).gt.zero.and.ui(its,isph).lt.one) then
      gi = zero
      do l = 0, lmax
        ind = l*l + l + 1
        fl = dble(l)
        fac = two*pi/(two*fl + one)
        do m = -l, l 
          !! DEBUG comment
          gi = gi + fac*basis(ind+m,its)*g(ind+m,isph)
        end do
      end do
      fii = w(its)*gi*y(its,isph)
      fx(1) = fx(1) + fii*zi(1,its,isph)
      fx(2) = fx(2) + fii*zi(2,its,isph)
      fx(3) = fx(3) + fii*zi(3,its,isph)
    end if
  end do

  ! second set of contributions:
  ! part of kb and ka
  do its = 1, ngrid

    ! run over all the spheres except isph 
    do jsph = 1, nsph
      if (ui(its,jsph).gt.zero .and. jsph.ne.isph) then
        ! build geometrical quantities
        cx = csph(1,jsph) + rsph(jsph)*grid(1,its)
        cy = csph(2,jsph) + rsph(jsph)*grid(2,its)
        cz = csph(3,jsph) + rsph(jsph)*grid(3,its)
        vji(1) = cx - csph(1,isph)
        vji(2) = cy - csph(2,isph)
        vji(3) = cz - csph(3,isph)
        vvji = sqrt(vji(1)*vji(1) + vji(2)*vji(2) + &
          &  vji(3)*vji(3))
        tji = vvji/rsph(isph)
        qji = one/vvji
        sji = vji/vvji

        ! build the jacobian of sji
        sjac = zero
        sjac(1,1) = - one
        sjac(2,2) = - one
        sjac(3,3) = - one
        do icomp = 1, 3
          do jcomp = 1, 3
            sjac(icomp,jcomp) = qji*(sjac(icomp,jcomp) &
              & + sji(icomp)*sji(jcomp))
          end do
        end do

        ! assemble the local basis and its gradient
        call dbasis(sji,basloc,dbsloc,vplm,vcos,vsin)

        ! assemble the contribution
        a = zero
        tt = one/(tji)
        do l = 0, lmax
          ind = l*l + l + 1
          fl = dble(l)
          fcl = - tt*four*pi*fl/(two*fl + one)
          do m = -l, l
            fac = fcl*g(ind+m,isph)
            b = (fl + one)*basloc(ind+m)/(rsph(isph)*tji)

            ! apply the jacobian to grad y
            va(1) = sjac(1,1)*dbsloc(1,ind+m) + &
              & sjac(1,2)*dbsloc(2,ind+m) + sjac(1,3)*dbsloc(3,ind+m)
            va(2) = sjac(2,1)*dbsloc(1,ind+m) + &
              & sjac(2,2)*dbsloc(2,ind+m) + sjac(2,3)*dbsloc(3,ind+m)
            va(3) = sjac(3,1)*dbsloc(1,ind+m) + &
              & sjac(3,2)*dbsloc(2,ind+m) + sjac(3,3)*dbsloc(3,ind+m)
            a(1) = a(1) + fac*(sji(1)*b + va(1))
            a(2) = a(2) + fac*(sji(2)*b + va(2))
            a(3) = a(3) + fac*(sji(3)*b + va(3))
          end do
          tt = tt/tji
        end do
        fac = ui(its,jsph)*w(its)*y(its,jsph)
        fx(1) = fx(1) - fac*a(1)
        fx(2) = fx(2) - fac*a(2)
        fx(3) = fx(3) - fac*a(3)
      end if
    end do
  end do

  ! ka contribution
  do its = 1, ngrid
    cx = csph(1,isph) + rsph(isph)*grid(1,its)
    cy = csph(2,isph) + rsph(isph)*grid(2,its)
    cz = csph(3,isph) + rsph(isph)*grid(3,its)
    a = zero

    ! iterate on all the spheres except isph
    do ksph = 1, nsph
      if (ui(its,isph).gt.zero .and. ksph.ne.isph) then
        ! geometrical stuff
        vik(1) = cx - csph(1,ksph)
        vik(2) = cy - csph(2,ksph)
        vik(3) = cz - csph(3,ksph)
        vvik = sqrt(vik(1)*vik(1) + vik(2)*vik(2) + & 
          & vik(3)*vik(3))
        tik = vvik/rsph(ksph)
        qik = one/vvik
        sik = vik/vvik

        ! build the jacobian of sik
        sjac = zero
        sjac(1,1) = one
        sjac(2,2) = one
        sjac(3,3) = one
        do icomp = 1, 3
          do jcomp = 1, 3
            sjac(icomp,jcomp) = qik*(sjac(icomp,jcomp) &
              & - sik(icomp)*sik(jcomp))
          end do
        end do

        ! if we are in the switching region, recover grad_i u_i
        vb = zero
        if (ui(its,isph).lt.one) then
          vb(1) = zi(1,its,isph)
          vb(2) = zi(2,its,isph)
          vb(3) = zi(3,its,isph)
        end if

        ! assemble the local basis and its gradient
        call dbasis(sik,basloc,dbsloc,vplm,vcos,vsin)

        ! assemble the contribution
        tt = one/(tik)
        do l = 0, lmax
          ind = l*l + l + 1
          fl = dble(l)
          fcl = - tt*four*pi*fl/(two*fl + one)
          do m = -l, l
            fac = fcl*g(ind+m,ksph)
            fac = - fac*basloc(ind+m)
            a(1) = a(1) + fac*vb(1)
            a(2) = a(2) + fac*vb(2) 
            a(3) = a(3) + fac*vb(3)

            fac = ui(its,isph)*fcl*g(ind+m,ksph)
            b = - (fl + one)*basloc(ind+m)/(rsph(ksph)*tik)

            ! apply the jacobian to grad y
            va(1) = sjac(1,1)*dbsloc(1,ind+m) + &
              & sjac(1,2)*dbsloc(2,ind+m) + sjac(1,3)*dbsloc(3,ind+m)
            va(2) = sjac(2,1)*dbsloc(1,ind+m) + &
              & sjac(2,2)*dbsloc(2,ind+m) + sjac(2,3)*dbsloc(3,ind+m)
            va(3) = sjac(3,1)*dbsloc(1,ind+m) + &
              & sjac(3,2)*dbsloc(2,ind+m) + sjac(3,3)*dbsloc(3,ind+m)
            a(1) = a(1) + fac*(sik(1)*b + va(1))
            a(2) = a(2) + fac*(sik(2)*b + va(2))
            a(3) = a(3) + fac*(sik(3)*b + va(3))
          end do
          tt = tt/tik
        end do
      end if
    end do
    fac = w(its)*y(its,isph)
    fx(1) = fx(1) - fac*a(1)
    fx(2) = fx(2) - fac*a(2)
    fx(3) = fx(3) - fac*a(3)
  end do
  end subroutine gradr

end module ddpcm_lib
