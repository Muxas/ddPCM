!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! wrapper for the linear solvers for IEFPCM. The IEFPCM equation we want to solve is
!
!     R_\eps \Phi_\eps = R_\infty \Phi.
! 
! the right-hand side is therefore g = R_\infty \Phi.
!
! input:
! 
!   star     logical, true:  solve the adjoint PCM equations,
!                     false: solve the PCM equatiosn
!
!   cart     logical, true:  the right-hand side for the PCM has to be assembled 
!                            inside this routine and the unscaled potential at the 
!                            external points of the cavity is provided in phi. 
!                     false: the right-hand side for the PCM equations is provided
!                            in glm.
!                     cart is not referenced if star is true. 
!
!   doprec   logical, true:  assemble the preconditioner 
!                     false: the preconditioner is already assembled and available
!
!   phi      real,    contains the potential at the external cavity points if star is
!                     false and cart is true.
!                     phi is not referenced in any other case.
!
!   glm      real,    contains the right-hand side for the PCM (adjoint) equations. 
!                     if star is false and cart is true, glm is not referenced. 
!
! output:
!
!   phi_eps: real,    the solution to the COSMO (adjoint) equations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! this routine performs the following operations:
!
!   - allocates memory for the linear solvers, and fixes dodiag.
!     This parameters controls whether the diagonal part of the matrix is considered 
!     in matvec, which depends on the solver used. It is false for jacobi_diis and
!     true for GMRES. 
!
!   - if star is false and cart is true, assembles the right-hand side for the PCM
!     equations. Note that for GMRES, a preconditioner is applied.
!
!   - computes a guess for the solution (using the preconditioner)
!
!   - calls the required iterative solver
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
subroutine pcm( star, cart, doprec, phi, glm, phi_eps )
!
      use ddcosmo , only : ncav, nsph, nbasis, iconv, isolver, zero, prec,       &
                           precm1, ngrid, lmax, wghpot, intrhs, ndiis, do_diag,  &
                           iout, iprint 
!
      implicit none
      logical,                         intent(in)    :: star, cart, doprec
      real*8,  dimension(ncav),        intent(in)    :: phi
      real*8,  dimension(nbasis,nsph), intent(in)    :: glm
      real*8,  dimension(nbasis,nsph), intent(inout) :: phi_eps
!
      integer              :: isph, istatus, n_iter, info, c1, c2, cr
      real*8               :: tol, r_norm
      logical              :: ok
!
      real*8, allocatable  :: g(:,:), rhs(:,:), work(:,:), x(:,:), u(:), ulm(:), basloc(:), &
                              vplm(:), vcos(:), vsin(:)
!
      integer, parameter   :: gmm = 20, gmj = 25
!
      external             :: rx, prx, precx, hnorm, rstarx, prstarx
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     set a few parameters for the solver and matvec routine:
!
      tol    = 10.0d0**(-iconv)
      n_iter = 300
!
!     initialize the timer:
!
      call system_clock(count_rate=cr)
      call system_clock(count=c1)
!
!     allocate additional memory for GMRES:

      if ( isolver.eq.1 ) then
!
!       allocate workspace
        allocate( work(nsph*nbasis,0:2*gmj+gmm+2 -1) , stat=istatus )
        if ( istatus.ne.0 ) then
          write(*,*) ' pcm: [1] failed allocation for GMRES'
          stop
        endif
!
!       initialize
        work = zero
!
      endif
!
!     
!     DIRECT PCM EQUATION  R_\eps \Phi_\eps = R_\infty \Phi
!     -----------------------------------------------------
      if ( .not.star ) then
!
!       allocate workspace for rhs
        allocate( rhs(nbasis,nsph) , stat=istatus )
        if ( istatus.ne.0 ) then
          write(*,*) ' pcm: [2] failed allocation'
        endif
!
!       initialize
        rhs = zero
!
!       if required, set up the preconditioner:
        if ( doprec ) then
!
!         allocate workspaces for preconditioner
          if ( allocated(prec) )    deallocate(prec)
          if ( allocated(precm1) )  deallocate(precm1)
!          
          allocate( prec(nbasis,nbasis,nsph), precm1(nbasis,nbasis,nsph) , stat=istatus )
          if ( istatus.ne.0 ) then
            write(*,*) ' pcm: [3] failed allocation of the preconditioner'
            stop
          endif
!
!         now, build the preconditioner
          do isph = 1,nsph
            call mkprec( isph, .true., prec(:,:,isph), precm1(:,:,isph) )
          enddo
!          
        endif 
!
!       assemble rhs by weighting the potential and integrating
        if ( cart ) then
!                
!         allocate workspaces
          allocate( g(ngrid,nsph), x(nbasis,nsph), u(ngrid), ulm(nbasis), basloc(nbasis), &
                    vplm(nbasis), vcos(lmax+1), vsin(lmax+1) , stat=istatus )
          if ( istatus.ne.0 ) then
            write(*,*) ' pcm: [4] failed allocation'
          endif
!
!         Start weighting the potential...
! 
          call wghpot( phi, g )
!
!         ... and compute its multipolar expansion
!
          do isph = 1,nsph
            call intrhs( isph, g(:,isph), x(:,isph) )
          enddo
!
!         now, apply R_\infty [ do_diag should be set to .true. !!! ] :
!
          do isph = 1,nsph
            call mkrvec( isph, zero, x, rhs(:,isph), ulm, u, basloc, vplm, vcos, vsin )
          enddo
!
!         GMRES
          if ( isolver.eq.1 ) then
!                  
!           apply preconditioner to rhs
            x = rhs
            call precx( nbasis*nsph, x, rhs )
!            
          endif
!
!         dellocate workspaces
          deallocate( g, x, u, ulm, basloc, vplm, vcos, vsin , stat=istatus )
          if ( istatus.ne.0 ) then
            write(*,*) ' pcm: [1] failed deallocation'
          endif
!
!       no need to manipulate rhs
        else
!
!         jacobi/diis
          if ( isolver.eq.0 ) then
!
            rhs = glm
!
!         GMRES 
          elseif ( isolver.eq.1 ) then
!
!           apply preconditioner
            call precx( nbasis*nsph, glm, rhs )
!            
          endif
        endif
!
!       assemble a guess:
!
        call precx( nsph*nbasis, rhs, phi_eps )
!
!       call the solver:
!
!       jacobi/diis
        if ( isolver.eq.0 ) then
!
!         exclude diagonal blocks from action of R_\eps
          do_diag = .false.
!
!         action of  diag^-1 :  precx
!         action of  offdiag :  rx
!          
          call jacobi_diis( nsph*nbasis, iprint, ndiis, 3, tol, rhs, phi_eps, n_iter, ok, rx, precx )
!
!         restore action of diagonal blocks
          do_diag = .true.
!
!       GMRES      
        elseif ( isolver.eq.1 ) then
!
!         the gmres solver can not handle preconditioners, so we will solve 
!  
!           P R_\eps \Phi_\eps = P g,
!
!         where P is a jacobi preconditioner. note thus the plx matrix-vector multiplication routine.
!
          call gmresr( (iprint.gt.0), nsph*nbasis, gmj, gmm, rhs, phi_eps, work, tol, 'abs', n_iter, r_norm, prx, info )
!
!         solver success flag
          ok = ( info.eq.0 )
!
        endif
!
!!!! esolv = pt5 * ((eps - one)/eps) * sprod(nsph*nbasis,sigma,psi)
!
!       deallocate workspace
        deallocate( rhs , stat=istatus )
        if ( istatus.ne.0 ) then
           write(*,*) 'pcm: [2] failed deallocation'
        endif
!
!       deallocate preconditioner
        if ( doprec ) then
          deallocate( prec, precm1 , stat=istatus )
          if ( istatus.ne.0 ) then
             write(*,*) 'pcm: [3] failed deallocation'
          endif
        endif
!
!
!     ADJOINT PCM EQUATION  R_\epsi^* \Phi_\eps = g
!     ---------------------------------------------
      else
!
!       allocate workspace for rhs
        allocate( rhs(nbasis,nsph) , stat=istatus )
        if ( istatus.ne.0 ) then
          write(6,*) 'pcm: [5] allocation failed '
          stop
        endif
!
!       if required, assemble the preconditioner:
!
        if (doprec) then
!
!         allocate workspaces for the preconditionner
          if ( allocated(prec) )    deallocate( prec )
          if ( allocated(precm1) )  deallocate( precm1 )
          allocate( prec(nbasis,nbasis,nsph), precm1(nbasis,nbasis,nsph) , stat=istatus )
          if ( istatus.ne.0 ) then
            write(*,*) ' pcm: [3] failed allocation of the preconditioner'
            stop
          endif
!
!         build the preconditioner
!
          do isph = 1,nsph
            call adjprec( isph, .true., prec(:,:,isph), precm1(:,:,isph) )
          enddo
!          
        endif
!
!       jacobi/diis
        if ( isolver.eq.0 ) then
!                
          rhs = glm
!          
!       GMRES          
        elseif (isolver .eq. 1) then
!
!         apply preconditioner to rhs
          call precx( nbasis*nsph, glm, rhs )
!          
        endif
!
!       assemble a guess:
!
        call precx( nbasis*nsph, rhs, phi_eps )
!
!       call the solver:
!
!       jacobi/diis
        if ( isolver.eq.0 ) then
!
!         exclude diagonal blocks from action of R_\eps^*
          do_diag = .false.
!          
!         action of  diag^-1 :  precx
!         action of  offdiag :  rstarx
!
          call jacobi_diis( nsph*nbasis, iprint, ndiis, 3, tol, rhs, phi_eps, n_iter, ok, rstarx, precx )
!
!         restore action of diagonal blocks
          do_diag = .true.
!
!       GMRES
        elseif ( isolver.eq.1 ) then
!
!       the gmres solver can not handle preconditioners, so we will solve 
!  
!         P R_\eps^* Phi_eps = P g,
!
!       where P is a jacobi preconditioner. note thus the plx matrix-vector multiplication routine.
!
          call gmresr( (iprint.gt.0), nsph*nbasis, gmj, gmm, rhs, phi_eps, work, tol, 'abs', n_iter, r_norm, prstarx, info )
!
!         solver success flag
          ok = ( info.eq.0 )
!
        endif
!
!       deallocate workspace
        deallocate( rhs , stat=istatus )
        if ( istatus.ne.0 ) then
           write(*,*) 'pcm: [4] failed deallocation'
        endif
!
!       deallocate preconditioner
        if ( doprec ) then
          deallocate( prec, precm1 , stat=istatus )
          if ( istatus.ne.0 ) then
             write(*,*) 'pcm: [5] failed deallocation'
          endif
        endif
!
      endif
!
!     deallocate GMRES workspace
      if ( isolver.eq.1 ) then
        deallocate( work , stat=istatus)
        if ( istatus.ne.0 ) then
          write(*,*) 'pcm: [6] failed deallocation'
        endif
      endif
!
!     check solution
      if ( .not.ok ) then
!              
        if ( star ) then
          write(iout,1020)
 1020     format('adjoint ddPCM did not converge! Aborting...')
        else
          write(iout,1021)
 1021     format('ddPCM did not converge! Aborting...')
        endif
        stop
!        
      endif
!
      call system_clock(count=c2)
!
!     printing
      if ( iprint.gt.0 ) then
!
        write(iout,*)
        if (star) then
          write(iout,1010) 'adjoint ', dble(c2-c1)/dble(cr)
 1010     format(' the solution to the ddPCM ',a,'R_\eps equations took ',f8.3,' seconds.')
        else
          write(iout,1010) '', dble(c2-c1)/dble(cr)
        endif
        write(iout,*)
!
      endif
!
!
endsubroutine pcm
