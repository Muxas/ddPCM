module ddpcm_lib
use ddcosmo, only: nbasis, nsph, ngrid, ncav, lmax, iconv, iprint, &
    & wghpot, intrhs, prtsph, zero, pt5, one, two, four, pi, basis, &
    & eps, csph, rsph, grid, w, ui, ndiis, sprod, ylmbas
use pcm_fmm
!use ddcosmo
implicit none

real*8, allocatable :: rx_prc(:,:,:)
real*8, allocatable :: rhs(:,:), phieps(:,:), xs(:,:)
real*8, allocatable :: g(:,:)
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
! Temporary variables to get list of all admissible pairs
integer :: lwork, iwork, jwork
integer, allocatable :: work(:, :)
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

  subroutine ddpcm(phi, psi, esolv)
  ! main ddpcm driver, given the potential at the exposed cavity
  ! points and the psi vector, computes the solvation energy
  implicit none
  real*8, intent(in) :: phi(ncav), psi(nbasis,nsph)
  real*8, intent(inout) :: esolv
  real*8 :: tol
  integer :: isph, n_iter
  logical :: ok
  external :: lx, ldm1x, hnorm
  
  allocate(rx_prc(nbasis,nbasis,nsph))
  allocate(rhs(nbasis,nsph),phieps(nbasis,nsph),xs(nbasis,nsph))
  allocate(g(ngrid,nsph))
  tol = 10.0d0**(-iconv)

  ! build the preconditioner
  call mkprec

  ! build the RHS
  !write(6,*) 'pot', ncav
  !do isph = 1, ncav
  !  write(6,*) phi(isph)
  !end do
  g = zero
  xs = zero
  call wghpot(phi,g)
  do isph = 1, nsph
    call intrhs(isph,g(:,isph),xs(:,isph))
  end do

  !call prtsph('phi',nsph,0,xs)
  !call prtsph('psi',nsph,0,psi)

  ! rinf rhs
  dodiag = .true.
  call rinfx(nbasis*nsph,xs,rhs)
  !call prtsph('rhs',nsph,0,rhs)

  ! solve the ddpcm linear system
  n_iter = 200
  dodiag = .false.
  call cpu_time(start_time)
  call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,rhs,phieps,n_iter, &
      & ok,rx,apply_rx_prec,hnorm)
  call cpu_time(finish_time)
  write(6,*) 'ddpcm solve time: ', finish_time-start_time
  !call prtsph('phie',nsph,0,phieps)

  ! solve the ddcosmo linear system
  n_iter = 200
  dodiag = .false.
  call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,phieps,xs,n_iter, &
      & ok,lx,ldm1x,hnorm)
  !call prtsph('x',nsph,0,xs)

  ! compute the energy
  esolv = pt5*sprod(nsph*nbasis,xs,psi)

  return
  end subroutine ddpcm

  subroutine ddpcm_fmm(phi, psi, esolv)
  ! FMM-accelerated ddpcm driver, given the potential at the exposed cavity
  ! points and the psi vector, computes the solvation energy
  implicit none
  real*8, intent(in) :: phi(ncav), psi(nbasis,nsph)
  real*8, intent(inout) :: esolv
  real*8 :: tol
  integer :: isph, n_iter
  logical :: ok
  external :: lx, ldm1x, hnorm
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
  lwork = nclusters*400 ! Some magic constant which might be changed
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

  ! Continue with ddpcm
  allocate(rx_prc(nbasis,nbasis,nsph))
  allocate(rhs(nbasis,nsph),phieps(nbasis,nsph),xs(nbasis,nsph))
  allocate(g(ngrid,nsph))
  tol = 10.0d0**(-iconv)

  ! build the preconditioner
  call mkprec

  ! build the RHS
  !write(6,*) 'pot', ncav
  !do isph = 1, ncav
  !  write(6,*) phi(isph)
  !end do
  g = zero
  xs = zero
  call wghpot(phi,g)
  do isph = 1, nsph
    call intrhs(isph,g(:,isph),xs(:,isph))
  end do

  !call prtsph('phi',nsph,0,xs)
  !call prtsph('psi',nsph,0,psi)

  ! rinf rhs
  dodiag = .true.
  call rinfx_fmm(nbasis*nsph,xs,rhs)
  !call prtsph('rhs',nsph,0,rhs)

  ! solve the ddpcm linear system
  n_iter = 200
  dodiag = .false.
  call cpu_time(start_time)
  call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,rhs,phieps,n_iter, &
      & ok,rx_fmm,apply_rx_prec,hnorm)
  call cpu_time(finish_time)
  write(6,*) 'ddpcm_fmm solve time: ', finish_time-start_time

  ! Print matvec statistics
  !call pcm_matvec_print_stats

  !call prtsph('phie',nsph,0,phieps)

  ! solve the ddcosmo linear system
  n_iter = 200
  dodiag = .false.
  call jacobi_diis(nsph*nbasis,iprint,ndiis,4,tol,phieps,xs,n_iter, &
      & ok,lx,ldm1x,hnorm)
  !call prtsph('x',nsph,0,xs)

  ! compute the energy
  esolv = pt5*sprod(nsph*nbasis,xs,psi)

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
  end subroutine ddpcm_fmm


  subroutine mkprec
  ! Assemble the diagonal blocks of the Reps matrix
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
    call DGETRF(nbasis,nbasis,rx_prc(:,:,isph),nbasis,ipiv,istatus)
    if (istatus.ne.0) then 
      write(6,*) 'LU failed in mkprc'
      stop
    end if
    call DGETRI(nbasis,rx_prc(:,:,isph),nbasis,ipiv,work, &
        & nbasis*nbasis,istatus)
    if (istatus.ne.0) then 
      write(6,*) 'Inversion failed in mkprc'
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
  ! Computes Y = Reps X =
  ! = (2*pi*(eps + 1)/(eps - 1) - D) X 
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  real*8 :: fac

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
  real*8, intent(inout) :: y(nbasis,nsph)
  real*8 :: fac

  call dx_fmm(n,x,y)
  y = -y

  if (dodiag) then
    fac = two*pi*(eps + one)/(eps - one)
    y = y + fac*x
  end if
  end subroutine rx_fmm


  subroutine rinfx(n,x,y)
  ! Computes Y = Rinf X = 
  ! = (2*pi -  D) X
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
  real*8 :: fac

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
  real*8, intent(inout) :: y(nbasis,nsph)
  real*8 :: fac
  !real*8 :: z(nbasis, nsph), norm, diff
  !integer :: i, j

  call dx_fmm(n,x,y)
  !call dx(n, x, z)
  !norm = 0
  !diff = 0
  !do i = 1, nbasis
  !  do j = 1, nsph
  !    norm = norm + z(i, j)**2
  !    diff = diff + (z(i, j)-y(i, j))**2
  !  end do
  !end do
  !write(*,*) 'error/norm of fmm matvec=', sqrt(diff), '/', sqrt(norm)
  y = -y
  !z = -z

  if (dodiag) then
    fac = two*pi
    y = y + fac*x
    !z = z + fac*x
    !norm = 0
    !diff = 0
    !do i = 1, nbasis
    !  do j = 1, nsph
    !    norm = norm + z(i, j)**2
    !    diff = diff + (z(i, j)-y(i, j))**2
    !  end do
    !end do
    !write(*,*) 'error/norm of fmm matvec (diag=T)=', sqrt(diff), '/', sqrt(norm)
  end if
  end subroutine rinfx_fmm

  subroutine dx(n,x,y)
  ! Computes Y = D X
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: x(nbasis,nsph)
  real*8, intent(inout) :: y(nbasis,nsph)
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
  do isph = 1, nsph
    ! compute the "potential" from the other spheres
    ! at the exposed Lebedv points of the i-th sphere 
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
  real*8, intent(inout) :: y(nbasis,nsph)
  real*8 :: fmm_x((pm+1)*(pm+1), nsph) ! Input for FMM
  real*8 :: fmm_y((pl+1)*(pl+1), nsph) ! Output of FMM
  integer :: isph, igrid, l, lind, m
  real*8 :: f, vts(ngrid), fourpi
  fmm_x(1:nbasis, :) = x
  fmm_x(nbasis+1:, :) = zero
  ! Do actual FMM matvec
  !call pcm_matvec_grid_fmm_fast(nsph, csph, rsph, ngrid, grid, w, vgrid, ui, &
  !    pm, pl, vscales, ind, nclusters, cluster, children, cnode, rnode, &
  !    nnfar, sfar, far, nnnear, snear, near, fmm_x, fmm_y)
  call pcm_matvec_grid_fmm_use_mat(nsph, csph, rsph, ngrid, grid, w, vgrid, &
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
    ! Cut fmm_y to fit into output y
    y = y + fmm_y(1:nbasis, :)
  else
    ! Cut fmm_y to fit into output y
    y = fmm_y(1:nbasis, :)
  end if
  !write(*,*) 'I did actual FMM matvec'
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

end module ddpcm_lib
