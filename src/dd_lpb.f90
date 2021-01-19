!> @copyright (c) 2020-2021 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file src/dd_core.f90
!! Core routines and parameters of entire ddX software
!!
!! @version 1.0.0
!! @author Abhinav Jha and Michele Nottoli
!! @date 2021-01-07

module dd_lpb
use dd_cosmo
use dd_core
use dd_operators
use dd_solvers
implicit none
!!
!! Logical variables for iterations, cosmo solver, and HSP solver
!!
logical :: first_out_iter
logical :: do_diag
integer :: matAB
!!
!! Hardcoded values
!!
integer, parameter :: lmax0 = 6
integer, parameter :: nbasis0 = 49
real(dp),  parameter :: epsp = 1.0d0
!!
!! Taken from Chaoyu's MATLAB code
!!
!! wij          : (?) Maybe Eq. (54) 
real(dp), allocatable :: wij(:,:)
real(dp), allocatable :: coefvec(:,:,:), Pchi(:,:,:), &
                              & Qmat(:,:,:), coefY(:,:,:)
!! SI_ri        : Bessel' function of first kind
!! DI_ri        : Derivative of Bessel' function of first kind
!! SK_ri        : Bessel' function of second kind
!! DK_ri        : Derivative of Bessel' function of second kind
!! NOTE: Tolerance and n_iter_gmres can be given by the user
!! tol_gmres    : Tolerance of GMRES iteration
!! n_iter_gmres : Maximum number of GMRES itertation

real(dp), allocatable :: SI_ri(:,:), DI_ri(:,:), SK_ri(:,:), &
                              & DK_ri(:,:), termimat(:,:)
real(dp)              :: tol_gmres, n_iter_gmres

contains
  !!
  !! ddLPB calculation happens here
  !! @param[in] dd_data : dd Data 
  !! @param[in] phi     : Boundary conditions
  !! @param[in] psi     : Electrostatic potential vector.
  !!                      Use of psi unknown
  !! @param[in] gradphi : Gradient of phi
  !!
  !! @param[out] sigma  : Solution of ddLPB
  !! @param[out] esolv  : Electrostatic solvation energy
  !!
  subroutine ddlpb(dd_data, phi, psi, gradphi, sigma, esolv, charge, ndiis, niter, iconv)
  ! main ddLPB
  implicit none
  type(dd_data_type), intent(in)  :: dd_data
  logical                         :: converged = .false.
  integer                         :: iteration = 1
  integer, intent(in)             :: ndiis
  integer, intent(in)             :: niter
  integer, intent(in)             :: iconv
  real(dp), intent(inout)    :: esolv
  real(dp)                   :: inc, old_esolv
  real(dp), intent(inout)    :: sigma(dd_data % nbasis,dd_data % nsph)
  real(dp), intent(in)       :: phi(dd_data % ncav), &
                                        & gradphi(3,dd_data % ncav)
  real(dp), intent(in)       :: psi(dd_data % nbasis, dd_data % nsph)
  real(dp), intent(in)       :: charge(dd_data % nsph)
  !!
  !! Xr         : Reaction potential solution (Laplace equation)
  !! Xe         : Extended potential solution (HSP equation)
  !! rhs_r      : Right hand side corresponding to Laplace equation
  !! rhs_e      : Right hand side corresponding to HSP equation
  !! rhs_r_init : Initial right hand side corresponding to Laplace equation
  !! rhs_e_init : Initial right hand side corresponding to HSP equation
  !!
  real(dp), allocatable ::   Xr(:,:), Xe(:,:), rhs_r(:,:), rhs_e(:,:), &
                                  & rhs_r_init(:,:), rhs_e_init(:,:)
  !!
  !! g      : Intermediate matrix for computation of g0
  !! f      : Intermediate matrix for computation of f0
  !! g0     : Vector associated to psi_0 Eq.(77) QSM19.SISC
  !! f0     : Vector associated to partial_n_psi_0 Eq.(99) QSM19.SISC
  !! lx     : External routine from matvec.f90. Used for Jacobi solver
  !! ldm1x  : External routine from matvec.f90. Used for Jacobi solver
  !! hnorm  : External routine from matvec.f90. Used for Jacobi solver
  !!          h^-1/2 norm of the increment on each sphere
  !! ok     : Boolean to check convergence of solver
  !! tol    : Tolerance for Jacobi solver
  !! n_iter : Number of iterative steps
  real(dp), allocatable :: g(:,:), f(:,:), g0(:), f0(:)
  integer                    :: isph
  integer                    :: i
  logical                    :: ok = .false.
  real(dp)              :: tol
  integer                    :: n_iter
  integer                    :: its
  !
  ! Allocate Bessel's functions of the first kind and the second kind
  ! and their derivatives
  !
  call ddlpb_init(dd_data)

  allocate(g(dd_data % ngrid,dd_data % nsph),&
           & f(dd_data % ngrid, dd_data % nsph))
  allocate(g0(dd_data % nbasis),f0(dd_data % nbasis))
  allocate(rhs_r(dd_data % nbasis, dd_data % nsph),&
           & rhs_e(dd_data % nbasis, dd_data % nsph))
  allocate(rhs_r_init(dd_data % nbasis, dd_data % nsph),&
           & rhs_e_init(dd_data % nbasis, dd_data % nsph))
  allocate(Xr(dd_data % nbasis, dd_data % nsph),&
           & Xe(dd_data % nbasis, dd_data % nsph))

  !do i = 1, ncav
  !  write(6,'(3F15.8)') phi(i), gradphi(:,i)
  !end do
  !stop

  ! Build the right hand side
  ! do i = 1, ncav
  !   write(6,'(4F20.10)') phi(i), gradphi(:,i)
  ! end do
  !
  !! wghpot: Weigh potential at cavity points. Comes from ddCOSMO
  !1         Intermediate computation of G_0 Eq.(77) QSM19.SISC
  !!
  !! @param[in]  phi : Boundary conditions (This is psi_0 Eq.(20) QSM19.SISC)
  !! @param[out] g   : Boundary conditions on solute-solvent boundary gamma_j_e
  !!
  call wghpot(dd_data, phi,g)
  !!
  !! wghpot_f : Intermediate computation of F_0 Eq.(75) from QSM19.SISC
  !!
  !! @param[in]  gradphi : Gradient of psi_0
  !! @param[out] f       : Boundary conditions scaled by characteristic function
  !!
  call wghpot_f(dd_data, gradphi,f)

  ! do isph = 1, nsph
  !   do i = 1, ngrid
  !     write(6,'(2F20.10)') g(i,isph), f(i,isph)
  !   end do
  ! end do
  
  !!
  !! Integrate Right hand side
  !! rhs_r_init: g0+f0
  !! rhs_e_init: f0
  !!
  do isph = 1, dd_data % nsph
    !! intrhs is a subroutine in dd_operators
    !! @param[in]  isph : Sphere number, used for output
    !! @param[in]  g    : Intermediate right side g
    !! @param[out] g0   : Integrated right side Eq.(77) in QSM19.SISC
    call intrhs(dd_data % iprint, dd_data % ngrid, &
                dd_data % lmax, dd_data % vwgrid, isph, g(:,isph), g0)
    call intrhs(dd_data % iprint, dd_data % ngrid, &
                dd_data % lmax, dd_data % vwgrid, isph,f(:,isph),f0)
    !! rhs 
    rhs_r_init(:,isph) = g0 + f0
    rhs_e_init(:,isph) = f0
  end do

  rhs_r = rhs_r_init
  rhs_e = rhs_e_init

  tol = 10.0d0**(-iconv)

  first_out_iter = .true.

  do while (.not.converged)

    !! Solve the ddCOSMO step
    !! A X_r = RHS_r (= G_X+G_0) 
    !! NOTE: Number of iterative steps can be user provided
    n_iter = 200
    !! Call Jacobi solver
    !! @param[in]      nsph*nylm : Size of matrix
    !! @param[in]      iprint    : Flag for printing
    !! @param[in]      ndiis     : Number of points to be used for 
    !!                            DIIS extrapolation. Set to 25 in ddCOSMO
    !! @param[in]      4         : Norm to be used to evaluate convergence
    !!                             4 refers to user defined norm. Here hnorm
    !! @param[in]      tol       : Convergence tolerance
    !! @param[in]      rhs_r     : Right-hand side
    !! @param[in, out] xr        : Initial guess to solution and final solution
    !! @param[in, out] n_iter    : Number of iterative steps
    !! @param[in, out] ok        : Boolean to check whether the solver converged
    !! @param[in]      lx        : External subroutine to compute matrix 
    !!                             multiplication, i.e., Lx_r, comes from matvec.f90
    !! @param[in]      ldm1x     : External subroutine to apply invert diagonal
    !!                             matrix to vector, i.e., L^{-1}x_r, comes from matvec.f90
    !! @param[in]      hnorm     : User defined norm, comes from matvec.f90
    call jacobi_diis(dd_data, dd_data % n, dd_data % iprint, ndiis, 4, tol, &
                     & rhs_r, Xr, n_iter, ok, lx, ldm1x, hnorm)
    call convert_ddcosmo(dd_data, 1, Xr)
    ! call print_ddvector('xr',xr)
  
    !! Solve ddLPB step
    !! B X_e = RHS_e (= F_0)
    call lpb_hsp(dd_data, rhs_e, Xe)
    ! call print_ddvector('xe',xe)
  
    !! Update the RHS
    !! / RHS_r \ = / g + f \ - / c1 c2 \ / X_r \
    !! \ RHS_e /   \ f     /   \ c1 c2 / \ X_e /
    call update_rhs(dd_data, rhs_r_init, rhs_e_init, rhs_r, rhs_e, Xr, Xe)
    ! call print_ddvector('rhs_r',rhs_r)
    ! call print_ddvector('rhs_e',rhs_e)

    !! Compute energy
    !! esolv = pt5*sprod(nsph*nylm,xr,psi)
    esolv = zero
    do isph = 1, dd_data % nsph
      esolv = esolv + pt5*charge(isph)*Xr(1,isph)*(one/(two*sqrt(pi)))
    end do

    !! Check for convergence
    inc = abs(esolv - old_esolv)/abs(esolv)
    old_esolv = esolv
    if ((iteration.gt.1) .and. (inc.lt.tol)) then
      write(6,*) 'Reach tolerance.'
      converged = .true.
    end if
    write(6,*) iteration, esolv, inc
    iteration = iteration + 1

    ! to be removed
    first_out_iter = .false.
  end do

  return
  end subroutine ddlpb
  !!
  !! Allocate Bessel's functions of the first kind and the second kind
  !! Uses the file bessel.f90
  !! @param[out] SI_ri : Bessel's function of the first kind
  !! @param[out] DI_ri : Derivative of Bessel's function of the first kind
  !! @param[out] SK_ri : Bessel's function of the second kind
  !! @param[out] DK_ri : Derivative of Bessel's function of the second kind
  !! @param[out] NM    : Highest order computed
  !!
  subroutine ddlpb_init(dd_data)
  use bessel
  implicit none
  type(dd_data_type), intent(in)  :: dd_data
  integer                         :: istatus, isph, NM
  allocate(SI_ri(0:dd_data % lmax, dd_data % nsph),&
           & DI_ri(0:dd_data % lmax, dd_data % nsph),&
           & SK_ri(0:dd_data % lmax, dd_data % nsph), &
           & DK_ri(0:dd_data % lmax, dd_data % nsph), &
           & termimat(0:dd_data % lmax, dd_data % nsph), stat=istatus)
  if (istatus.ne.0) then
    write(*,*)'ddinit : [1] allocation failed !'
    stop
  end if
  do isph = 1, dd_data % nsph
    call SPHI_bessel(dd_data % lmax,dd_data % rsph(isph)*dd_data % kappa,&
                     & NM,SI_ri(:,isph), &
                     & DI_ri(:,isph))
    call SPHK_bessel(dd_data % lmax,dd_data % rsph(isph)*dd_data % kappa,&
                     & NM,SK_ri(:,isph), &
                     & DK_ri(:,isph))
  end do
  return
  end subroutine ddlpb_init

  !!
  !! Find intermediate F0 in the RHS of the ddLPB model given in Eq.(82)
  !! @param[in]  gradphi : Gradient of psi_0
  !! @param[out] f       : Intermediate calculation of F0
  !!
  subroutine wghpot_f(dd_data, gradphi, f )
  use bessel
  implicit none
  type(dd_data_type), intent(in)  :: dd_data
  real(dp), dimension(3, dd_data % ncav),       intent(in)  :: gradphi
  real(dp), dimension(dd_data % ngrid, dd_data % nsph),    intent(out) :: f

  integer :: isph, ig, ic, ind, ind0, jg, l, m, jsph
  real(dp) :: nderphi, sumSijn, rijn, coef_Ylm, sumSijn_pre, termi, &
      & termk, term
  real(dp), dimension(3) :: sijn, vij
  real(dp) :: rho, ctheta, stheta, cphi, sphi
  real(dp), allocatable :: SK_rijn(:), DK_rijn(:)

  integer :: l0, m0, NM, kep, istatus
  real(dp), dimension(dd_data % nbasis, dd_data % nsph) :: c0
  real(dp), dimension(0:dd_data % lmax, dd_data % nsph) :: coef_bessel
  real(dp), allocatable :: vplm(:), basloc(:), vcos(:), vsin(:)

  ! initialize
  allocate(vplm(dd_data % nbasis),basloc(dd_data % nbasis),vcos(dd_data % lmax+1),vsin(dd_data % lmax+1))
  allocate(SK_rijn(0:lmax0),DK_rijn(0:lmax0))
  ic = 0 ; f(:,:)=0.d0
  !
  ! Compute c0 Eq.(98) QSM19.SISC
  !
  do isph = 1, dd_data % nsph
    do ig = 1, dd_data % ngrid
      if ( dd_data % ui(ig,isph).ne.zero ) then
        ic = ic + 1
        nderphi = dot_product( gradphi(:,ic),dd_data % cgrid(:,ig) )
        c0(:, isph) = c0(:,isph) + dd_data % wgrid(ig)*dd_data % ui(ig,isph)*nderphi*dd_data % vgrid(:,ig)
      end if
    end do
  end do

  allocate (coefY(dd_data % ncav, nbasis0, dd_data % nsph), stat = istatus)
  ! memuse = memuse + ncav*nbasis0*nsph
  ! memmax = max(memmax,memuse)
  if ( istatus .ne. 0 ) then
      write(*,*)'wghpot_f : [1] allocation failed!'
      stop
  end if
  !
  ! Compute coef_bessel
  ! Here (der_i_l0/i_l0 - der_k_l0/k_l0)^(-1) is computed in Eq.(97)
  ! Here we consider three cases. One is when kappa is too small,
  ! kappa is too big, and kappa is moderate. To find these values,
  ! we use the formulas given in the book, Mathematical Methods for 
  ! Physicsts, Arfken and Weber.
  !
  do jsph = 1, dd_data % nsph
    do l0 = 0, lmax0
      ! Case: kappa > tol_inf
      if (max(DI_ri(l0,jsph), SI_ri(l0,jsph)).gt.tol_inf) then
        termi = dd_data % kappa
      ! Case: kappa < tol_zero
      !       One can ignore the other terms (including the second
      !       term given below) as kappa is so small the
      !       contributions are negligible
      else if (min(DI_ri(l0,jsph), SI_ri(l0,jsph)).lt.tol_zero) then
        termi = l0/dd_data % rsph(jsph) + &
            & (l0 + one)*(dd_data % kappa**2*dd_data % rsph(jsph))/((two*l0 + one) * &
            & (two*l0 + three))
      ! Case: kappa is of normal size.
      ! NOTE: We notice a factor of kappa. The reason being while computing
      !       DI_ri the argument is kappa*r. Hence a factor of kappa.
      else
        termi = DI_ri(l0,jsph)/SI_ri(l0,jsph)*dd_data % kappa
      end if
      !write(*,*) SI_ri(l0,jsph), termi

      ! Similar calculation for SK_ri to SI_ri
      if (SK_ri(l0,jsph).gt.tol_inf) then
        termk = - (l0 + one)/dd_data % rsph(jsph) - &
            & l0*(dd_data % kappa**2*dd_data % rsph(jsph))/((two*l0 - one)*(two*l0 + one))
      else if (SK_ri(l0,jsph).lt.tol_zero) then
        termk = -dd_data % kappa
      else
        termk = DK_ri(l0,jsph)/SK_ri(l0,jsph)*dd_data % kappa
      end if

      !write(*,*) SK_ri(l0,jsph), termk
      coef_bessel(l0,jsph) = one/(termi - termk)
      !write(*,*) DI_ri(l0,jsph), SI_ri(l0,jsph), coef_bessel(l0,jsph)
      !write(*,*) (min(-DK_ri(l0,jsph), SK_ri(l0,jsph)).lt.tol_zero), &
      !    & DK_ri(l0,jsph), termk
    end do
  end do

  ! Computation of F0 using above terms
  ! (?) Use of kep not known (?)
  kep = 0
  do isph = 1, dd_data % nsph
    do ig = 1, dd_data % ngrid
      if (dd_data % ui(ig,isph).gt.zero) then
        kep = kep + 1
        sumSijn = zero
        ! Loop to compute Sijn
        do jsph = 1, dd_data % nsph
          ! (?) Use of sumSijn_pre (?)
          sumSijn_pre = sumSijn
          vij  = dd_data % csph(:,isph) + dd_data % rsph(isph)*dd_data % cgrid(:,ig) - dd_data % csph(:,jsph)
          rijn = sqrt(dot_product(vij,vij))
          sijn = vij/rijn
          
          ! Compute Bessel function of 2nd kind for the coordinates
          ! (s_ijn, r_ijn) and compute the basis function for s_ijn
          call SPHK_bessel(lmax0,rijn*dd_data % kappa,NM,SK_rijn,DK_rijn)
          call ylmbas(sijn , rho, ctheta, stheta, cphi, &
                      & sphi, dd_data % lmax, dd_data % vscales, &
                      & basloc, vplm, vcos, vsin)

          do l0 = 0,lmax0
            ! term: k_l0(r_ijn)/k_l0(r_i)
            ! Uses approximation of SK_ri in double factorial terms when argument
            ! is less than one
            if (SK_ri(l0,jsph).gt.tol_inf) then
            ! Uses approximation of SK_ri in double factorial terms when argument
            ! is greater than l0
              term = (dd_data % rsph(jsph)/rijn)**(l0+1)
            else if (SK_ri(l0,jsph).lt.tol_zero) then
              term = (dd_data % rsph(jsph)/rijn)*exp(-dd_data % kappa*(rijn-dd_data % rsph(jsph)))
            else
              term = SK_rijn(l0)/SK_ri(l0,jsph)
            end if
            ! coef_Ylm : (der_i_l0/i_l0 - der_k_l0/k_l0)^(-1)*k_l0(r_ijn)/k_l0(r_i)
            coef_Ylm =  coef_bessel(l0,jsph)*term
            do m0 = -l0, l0
              ind0 = l0**2 + l0 + m0 + 1
              sumSijn = sumSijn + c0(ind0,jsph)*coef_Ylm*basloc(ind0)
              ! NOTE: (?) Don't know the use of coefY till now (?)
              coefY(kep,ind0,jsph) = coef_Ylm*basloc(ind0)
            end do
          end do
        end do
        !
        ! Here Intermediate value of F_0 is computed Eq. (99)
        ! Mutilplication with Y_lm and weights will happen afterwards
        !write(6,*) sumSijn, epsp, eps, dd_data % ui(ig,isph)
        f(ig,isph) = -(epsp/dd_data % eps)*dd_data % ui(ig,isph) * sumSijn
      end if
    end do
  end do 

  deallocate( vplm, basloc, vcos, vsin, SK_rijn, DK_rijn  )
  return
  end subroutine wghpot_f

  !
  ! Subroutine used for the GMRES solver
  ! NOTE: It is refered as matABx in the GMRES solver.
  !       Fortran is not case sensitive
  ! @param[in]      n : Size of the matrix
  ! @param[in]      x : Input vector
  ! @param[in, out] y : y=A*x
  !
  subroutine matABx(dd_data, n, x, y )
  implicit none 
  type(dd_data_type), intent(in)  :: dd_data
  integer, intent(in) :: n
  real(dp), dimension(dd_data % nbasis, dd_data % nsph), intent(in) :: x
  real(dp), dimension(dd_data % nbasis, dd_data % nsph), intent(inout) :: y
  integer :: isph, istatus
  real(dp), allocatable :: pot(:), vplm(:), basloc(:), vcos(:), vsin(:)
  integer :: i
  ! allocate workspaces
  allocate( pot(dd_data % ngrid), vplm(dd_data % nbasis), basloc(dd_data % nbasis), &
            & vcos(dd_data % lmax+1), vsin(dd_data % lmax+1) , stat=istatus )
  if ( istatus.ne.0 ) then
    write(*,*) 'Bx: allocation failed !'
    stop
  endif

  y = zero
  do isph = 1, dd_data % nsph
    call calcv2_lpb(dd_data, isph, pot, x, basloc, vplm, vcos, vsin )
    ! intrhs comes from ddCOSMO
    call intrhs(dd_data % iprint, dd_data % ngrid, &
                dd_data % lmax, dd_data % vwgrid, isph, pot, y(:,isph) )
    ! Action of off-diagonal blocks
    y(:,isph) = - y(:,isph)
    ! Add action of diagonal block
    y(:,isph) = y(:,isph) + x(:,isph)
  end do

  deallocate( pot, basloc, vplm, vcos, vsin , stat=istatus )
  if ( istatus.ne.0 ) then
    write(*,*) 'matABx: allocation failed !'
    stop
  endif
  end subroutine matABx

  !!
  !! Scale the ddCOSMO solution vector
  !! @param[in]      direction : Direction of the scaling
  !! (?) NOTE: We send direction = 1, then why the case of -1 (?)
  !! @param[in, out] vector    : ddCOSMO solution vector
  !!
  subroutine convert_ddcosmo(dd_data, direction, vector)
  implicit none
  type(dd_data_type), intent(in)  :: dd_data
  integer :: isph, l, m, ind
  integer, intent(in) :: direction
  real(dp), intent(inout) :: vector(dd_data % nbasis, dd_data % nsph)
  real(dp) :: fac
  
  do isph = 1, dd_data % nsph
    do l = 0, dd_data % lmax
      ind = l*l + l + 1
      fac = four*pi/(two*dble(l) + one) 
      if (direction.eq.-1) fac = one/fac
      do m = -l, l
        vector(ind + m,isph) = fac*vector(ind + m,isph)
      end do
    end do
  end do
  return
  end subroutine convert_ddcosmo


  !!
  !! Solve the HSP problem
  !! @param[in]      rhs : Right hand side for HSP
  !! @param[in, out] Xe  : Solution vector
  !!
  subroutine lpb_hsp(dd_data, rhs, Xe)
  implicit none
  type(dd_data_type), intent(in)  :: dd_data
  real(dp), dimension(dd_data % nbasis, dd_data % nsph), intent(in) :: rhs
  real(dp), dimension(dd_data % nbasis, dd_data % nsph), intent(inout) :: Xe
  integer :: isph, istatus, n_iter, info, c1, c2, cr
  real(dp) :: tol, r_norm
  real(dp), allocatable :: work(:,:)
  integer, parameter  :: gmm = 20, gmj = 25
  
  allocate(work(dd_data % nsph*dd_data % nbasis, 0:2*gmj + gmm + 2 - 1),stat=istatus)
  if (istatus.ne.0) then
    write(*,*) ' LPB-HSP: failed allocation for GMRES'
    stop
  endif
   
  work = zero
  Xe = rhs
  
  matAB = 1
  tol_gmres = 1.0d-8
  !!
  !! Call GMRES solver
  !! @param[in]      Residue_print : Set to false by default. Prints the
  !!                                 intermediate residue.
  !! @param[in]      nsph*nylm     : Size of matrix
  !! @param[in]      gmj           : Integer truncation parameter, Default: 25
  !! @param[in]      gmm           : Integer dimension of the GMRES,
  !!                                 Default: 20
  !! @param[in]      rhs           : Right hand side
  !! @param[in, out] Xe            : Initial guess of the problem
  !! @param[in]      work          : Work space of size
  !!                               : nsph*nylm X (2*gmj + gmm + 2)
  !! @param[in]      tol_gmres     : GMRES tolerance
  !! @param[in]      Stopping      : Stopping criteria, Default set to
  !!                                 'rel' for relative. Other option
  !!                                 'abs' for absolute
  !! @param[in]      n_iter_gmres  : Number of GMRES iteration
  !! @param[in]      r_norm        : Residual measure
  !! @param[in]      matABx        : Subroutine A*x. Named matabx in file
  !! @param[in, out] info          : Flag after solve. 0 means within tolerance
  !!                                 1 means max number of iteration
  call gmresr(dd_data, .false., dd_data % nsph*dd_data % nbasis, gmj, gmm, rhs, Xe, work, tol_gmres,'rel', &
      & n_iter_gmres, r_norm, matABx, info)

  deallocate(work)
  endsubroutine lpb_hsp

  !
  ! Intermediate computation of BX_e
  ! @param[in]      isph   : Number of the sphere
  ! @param[in, out] pot    : Array of size ngrid
  ! @param[in]      x      : Input vector (Usually X_e)
  ! @param[in, out] basloc : Used to compute spherical harmonic
  ! @param[in, out] vplm   : Used to compute spherical harmonic
  ! @param[in, out] vcos   : Used to compute spherical harmonic
  ! @param[in, out] vsin   : Used to compute spherical harmonic
  !
  subroutine calcv2_lpb (dd_data, isph, pot, x, basloc, vplm, vcos, vsin )
  type(dd_data_type), intent(in) :: dd_data
  integer, intent(in) :: isph
  real(dp), dimension(dd_data % nbasis, dd_data % nsph), intent(in) :: x
  real(dp), dimension(dd_data % ngrid), intent(inout) :: pot
  real(dp), dimension(dd_data % nbasis), intent(inout) :: basloc
  real(dp), dimension(dd_data % nbasis), intent(inout) :: vplm
  real(dp), dimension(dd_data % lmax+1), intent(inout) :: vcos
  real(dp), dimension(dd_data % lmax+1), intent(inout) :: vsin
  real(dp), dimension(dd_data % nbasis) :: fac_cosmo, fac_hsp
  integer :: its, ij, jsph
  real(dp) :: rho, ctheta, stheta, cphi, sphi
  real(dp) :: vij(3), sij(3)
  real(dp) :: vvij, tij, xij, oij

  pot = zero
  do its = 1, dd_data % ngrid
    if (dd_data % ui(its,isph).lt.one) then
      do ij = dd_data % inl(isph), dd_data % inl(isph+1)-1
        jsph = dd_data % nl(ij)

        ! compute geometrical variables
        vij  = dd_data % csph(:,isph) + dd_data % rsph(isph)*dd_data % cgrid(:,its) - dd_data % csph(:,jsph)
        vvij = sqrt(dot_product(vij,vij))
        tij  = vvij/dd_data % rsph(jsph) 

        if ( tij.lt.one ) then
          sij = vij/vvij
          call ylmbas(sij, rho, ctheta, stheta, cphi, &
                      & sphi, dd_data % lmax, &
                      & dd_data % vscales, basloc, &
                      & vplm, vcos, vsin)
          call inthsp(dd_data, vvij, dd_data % rsph(jsph), jsph, basloc, fac_hsp)
          xij = fsw(tij, dd_data % se, dd_data % eta)
          if (dd_data % fi(its,isph).gt.one) then
            oij = xij/dd_data % fi(its, isph)
          else
            oij = xij
          end if
          pot(its) = pot(its) + oij*dot_product(fac_hsp,x(:,jsph))
        end if
      end do
    end if
  end do
  endsubroutine calcv2_lpb


  !
  ! Intermediate calculation in calcv2_lpb subroutine
  ! @param[in]  rijn    : Radius of sphers x_ijn
  ! @param[in]  ri      : Radius of sphers x_i
  ! @param[in]  isph    : Index of sphere
  ! @param[in]  basloc  : Spherical Harmonic
  ! @param[out] fac_hsp : Return bessel function ratio multiplied by 
  !                       the spherical harmonic Y_l'm'. Array of size nylm
  !
  subroutine inthsp(dd_data, rijn, ri, isph, basloc, fac_hsp)
  use bessel
  implicit none
  type(dd_data_type), intent(in)  :: dd_data
  integer, intent(in) :: isph
  real(dp), intent(in) :: rijn, ri
  real(dp), dimension(dd_data % nbasis), intent(in) :: basloc
  real(dp), dimension(dd_data % nbasis), intent(inout) :: fac_hsp
  real(dp), dimension(0:dd_data % lmax) :: SI_rijn, DI_rijn
  integer :: l, m, ind, NM

  SI_rijn = 0
  DI_rijn = 0
  fac_hsp = 0

  ! Computation of modified spherical Bessel function values      
  call SPHI_bessel(dd_data % lmax,rijn*dd_data % kappa,NM,SI_rijn,DI_rijn)
  
  do l = 0, dd_data % lmax
    do  m = -l, l
      ind = l*l + l + 1 + m
      if ((SI_rijn(l).lt.zero) .or. (SI_ri(l,isph).lt.tol_zero) &
          & .or. (SI_rijn(l)/SI_ri(l,isph).gt.(rijn/ri)**l)) then
        fac_hsp(ind) = (rijn/ri)**l*basloc(ind)
      else if ( SI_ri(l,isph).gt.tol_inf) then
        fac_hsp(ind) = zero
      else
        fac_hsp(ind) = SI_rijn(l)/SI_ri(l,isph)*basloc(ind)
      end if
    end do
  end do
  endsubroutine inthsp


  subroutine intcosmo(dd_data, tij, basloc, fac_cosmo)
  implicit none
  type(dd_data_type), intent(in)  :: dd_data
  real(dp),  intent(in) :: tij
  real(dp), dimension(dd_data % nbasis), intent(in) :: basloc
  real(dp), dimension(dd_data % nbasis), intent(inout) :: fac_cosmo
  integer :: l, m, ind
  do l = 0, dd_data % lmax
    do  m = -l, l
        ind = l*l + l + 1 + m
        fac_cosmo(ind) = tij**l*basloc(ind)
    end do
  end do
  end subroutine intcosmo

  !
  ! Update the RHS in outer iteration
  ! @param[in] rhs_cosmo_init : G_0
  ! @param[in] rhs_hsp_init   : F_0
  ! @param[in, out] rhs_cosmo : -C_1*X_r^(k-1) - -C_2*X_e^(k-1) + G_0 + F_0
  ! @param[in, out] rhs_hsp   : -C_1*X_r^(k-1) - -C_2*X_e^(k-1) + F_0
  ! @param[in] Xr             : X_r^(k-1)
  ! @param[in] Xe             : X_e^(k-1)
  !
  subroutine update_rhs(dd_data, rhs_cosmo_init, rhs_hsp_init, rhs_cosmo, & 
      & rhs_hsp, Xr, Xe)
  use dd_cosmo
  use bessel
  implicit none
  type(dd_data_type), intent(in)  :: dd_data
  real(dp), dimension(dd_data % nbasis,dd_data % nsph), intent(in) :: rhs_cosmo_init, &
      & rhs_hsp_init
  real(dp), dimension(dd_data % nbasis,dd_data % nsph), intent(inout) :: rhs_cosmo, rhs_hsp
  real(dp), dimension(dd_data % nbasis,dd_data % nsph) :: rhs_plus
  real(dp), dimension(dd_data % nbasis,dd_data % nsph), intent(in) :: Xr, Xe
  integer :: isph, jsph, ig, kep, ind, l1,m1, ind1, ind0, count, istatus
  real(dp), dimension(3) :: vij
  real(dp), dimension(dd_data % nbasis,dd_data % nsph) :: diff_re
  real(dp), dimension(nbasis0,dd_data % nsph) :: diff0
  real(dp), dimension(dd_data % nbasis,dd_data % nbasis,dd_data % nsph) :: smat
  real(dp), dimension(dd_data % ncav) :: diff_ep
  real(dp) :: Qval, rijn, val
  integer :: c0, cr, c_qmat, c_init, c_ep0, c_ep1 !, nbasis_appro
      
  if (first_out_iter) then
    allocate(coefvec(dd_data % ngrid, dd_data % nbasis, dd_data % nsph), &
            Pchi(dd_data % nbasis, nbasis0, dd_data % nsph), stat = istatus)
    if (istatus.ne.0) then
      write(*,*)'update_rhs : [1] allocation failed!'
      stop
    end if
  end if

  ! Compute P_chi matrix, Eq.(87)
  ! TODO: probably has to be declared somewhere
  ! and i have to recover mkpmat 
  if (first_out_iter) then      
    do jsph = 1, dd_data % nsph
      call mkpmat(dd_data, jsph, Pchi(:,:,jsph))
    end do    
  end if 
      
  ! Precompute coefvec of Qmat, Cost: linear scaling
  ! TODO: remove all the precomputations
  ! or, if they are really necessary, do them in a separate subroutine
  if (first_out_iter) then 
    do isph = 1,dd_data % nsph
      do ig = 1,dd_data % ngrid
        if (dd_data % ui(ig, isph) .gt. 0) then
          do ind  = 1, dd_data % nbasis 
            coefvec(ig,ind,isph) = dd_data % wgrid(ig)*dd_data % ui(ig,isph)*dd_data % vgrid(ind,ig)
          end do
        end if
      end do
    end do
  end if
      
  ! precompute termimat
  ! TODO: same as before
  if (first_out_iter) then
    do jsph = 1, dd_data % nsph
      do l1 = 0, dd_data % lmax
        if (max(DI_ri(l1,jsph),SI_ri(l1,jsph)).gt.tol_inf) then
          termimat(l1,jsph) = dd_data % kappa
        else if (min(DI_ri(l1,jsph),SI_ri(l1,jsph)).lt.tol_zero) then
          termimat(l1,jsph) = l1/dd_data % rsph(jsph) + &
              & (l1 + 1)*(dd_data % kappa**2*dd_data % rsph(jsph))/((two*l1 + one) * &
              & (two*l1 + three))
        else
          termimat(l1,jsph) = DI_ri(l1,jsph)/SI_ri(l1,jsph)*dd_data % kappa
        end if
      end do
    end do
  end if

  if (first_out_iter) then
    if (dd_data % iprint .gt. 0) then  
      write(*,999) dble(c_init-c0)/dble(cr)
 999  format('Time of initializing Pchi, coefvec, termi: ',f8.3, &
          & ' secs.')
    end if
  end if

  ! diff_re = epsp/eps*l1/ri*Xr - i'(ri)/i(ri)*Xe,
  diff_re = zero 
  do jsph = 1, dd_data % nsph
    do l1 = 0, dd_data % lmax
      do m1 = -l1,l1
        ind1 = l1**2 + l1 + m1 + 1
        diff_re(ind1,jsph) = ((epsp/dd_data % eps)*(l1/dd_data % rsph(jsph)) * &
            & Xr(ind1,jsph) - termimat(l1,jsph)*Xe(ind1,jsph))
      end do
    end do
  end do

  ! diff0 = Pchi * diff_er, linear scaling
  ! TODO: probably doing PX on the fly is better 
  diff0 = zero 
  do jsph = 1, dd_data % nsph
    do ind0 = 1, nbasis0
      diff0(ind0, jsph) = dot_product(diff_re(:,jsph), &
          & Pchi(:,ind0, jsph))
    end do
  end do

  ! diff_ep = diff0 * coefY,    COST: M^2*nbasis*Nleb
  diff_ep = zero
  do kep = 1, dd_data % ncav
    val = zero
    do jsph = 1, dd_data % nsph 
      do ind0 = 1, nbasis0
        val = val + diff0(ind0,jsph)*coefY(kep,ind0,jsph)
      end do
    end do
    diff_ep(kep) = val 
  end do

  rhs_plus = zero
  kep = 0
  do isph = 1, dd_data % nsph
    do ig = 1, dd_data % ngrid
      if (dd_data % ui(ig,isph).gt.zero) then 
        kep = kep + 1
        do ind = 1, dd_data % nbasis
          rhs_plus(ind,isph) = rhs_plus(ind,isph) + &
              & coefvec(ig,ind,isph)*diff_ep(kep)
        end do
      end if
    end do
  end do

  rhs_cosmo = rhs_cosmo_init - rhs_plus
  rhs_hsp = rhs_hsp_init - rhs_plus 

  return
  end subroutine update_rhs  

  !
  ! Computation of P_chi
  ! @param[in]  isph : Sphere number
  ! @param[out] pmat : Matrix of size nylm X (lmax0+1)^2, Fixed lmax0
  !
  subroutine mkpmat(dd_data, isph, pmat )
  implicit none
  type(dd_data_type), intent(in)  :: dd_data
  integer,  intent(in) :: isph
  real(dp), dimension(dd_data % nbasis, (lmax0+1)**2), intent(inout) :: pmat
  integer :: l, m, ind, l0, m0, ind0, its, nbasis0
  real(dp)  :: f, f0

  pmat(:,:) = zero
  do its = 1, dd_data % ngrid
    if (dd_data % ui(its,isph).ne.0) then
      do l = 0, dd_data % lmax
        ind = l*l + l + 1
        do m = -l,l
          f = dd_data % wgrid(its) * dd_data % vgrid(ind+m,its) * dd_data % ui(its,isph)
          do l0 = 0, lmax0
            ind0 = l0*l0 + l0 + 1
            do m0 = -l0, l0
              f0 = dd_data % vgrid(ind0+m0,its)
              pmat(ind+m,ind0+m0) = pmat(ind+m,ind0+m0) + f * f0
            end do
          end do
        end do
      end do
    end if
  end do
  endsubroutine mkpmat

end module
