!-------------------------------------------------------------------------------
! Purpose : solve PCM equation R_eps S \sigma = R_oo \Phi by solving the 
!           system :
!
!             { R_eps \Phi_eps = R_oo \Phi
!             {
!             {       S \sigma = \Phi_eps    ( <== COSMO )
!
!           After numerical discretization, we obtain linear systems 
!
!             { A_eps W = A_oo g
!             {
!             { L sigma = W
!
! Arguments :
!
!   sigma_g - IN : g ; OUT : sigma
!-------------------------------------------------------------------------------
subroutine iefpcm(expot, phi, psi, sigma_g, phi_eps , esolv)
!
      use  ddcosmo
!      
      implicit none
!
      real*8, dimension(ncav),        intent(in)    :: expot
      real*8, dimension( ngrid,nsph), intent(in)    :: phi
      real*8, dimension(nbasis,nsph), intent(in)    :: psi
      real*8, dimension(nbasis,nsph), intent(inout) :: sigma_g
      real*8, dimension(nbasis,nsph), intent(out)   :: phi_eps
      real*8,                         intent(inout) :: esolv
!      
!     P. Gatto, Nov 2016      
!     real*8, dimension(ngrid, nsph), intent(inout) :: sigma_g
!
!     local arrays:
      real*8, allocatable :: philm(:,:), wlm(:,:), glm(:,:), vold(:,:), &
                             wlm_new(:,:)
!
!     scratch arrays:
      real*8, allocatable :: basloc(:), vplm(:), vcos(:), vsin(:), xlm(:), x(:)
!
!     diis arrays:
      real*8, allocatable :: xdiis(:,:,:), ediis(:,:,:), bmat(:)
!
      integer :: it, isph, nmat, lenb, istatus, j, jsph,l,m,ind,nt,ns,its
      real*8  :: ene, vrms, vmax, tol, s1, s2, s3, tt, fep
      real*8, allocatable :: err(:), ddiag(:)
      logical :: dodiis
!
      real*8, allocatable :: dphi(:,:,:,:), f(:,:)
!
      integer, parameter :: nitmax=300
      real*8,  parameter :: ten=10.0d0, tredis=1.d-2
!
      real*8 :: pot(ngrid), Awlm(nbasis)

      real*8 :: philm_fmm(nbasis,nsph)
      real*8 :: xs(nsph),ys(nsph),zs(nsph)
      real*8 :: xt(ngrid*nsph),yt(ngrid*nsph),zt(ngrid*nsph),fmm_vec(ngrid*nsph)
      real*8 :: phi_j(nbasis,nsph),ggrid(ngrid,nsph),xx(ngrid,nsph)
      real*8 :: voldgrid(ngrid,nsph),wlm_fmm(nbasis,nsph)
      real*8 :: basloc_fmm(nbasis,ngrid),vrmsold

!!!      logical,parameter :: use_fmm = .true.
!
!-------------------------------------------------------------------------------
!
!     set solver tolerance
      tol = ten**(-iconv)
!      
!     set up DIIS iterative solver
      dodiis = .false.
      nmat   = 1
      lenb   = ndiis + 1
!
!     allocate memory
      allocate( err(   nbasis)          , &
                xlm(   nbasis)          , &
                x(     ngrid)           , &
                vold(  nbasis,nsph)     , &
                philm( nbasis,nsph)     , &
                wlm(   nbasis,nsph)     , &
                wlm_new(   nbasis,nsph) , &
                glm(   nbasis,nsph)     , &
                basloc(nbasis)          , &  ! spherical harmonics (SH)
                vplm(  nbasis)          , &  ! service array for SH
                vcos(  lmax+1)          , &  ! service array for SH
                vsin(  lmax+1)          , &  ! service array for SH
                xdiis(nbasis,nsph,ndiis), &
                ediis(nbasis,nsph,ndiis), &
                bmat(lenb*lenb)         , &
                stat=istatus )
!               prec(nbasis,nbasis,nsph), precm1(nbasis,nbasis,nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'iefpcm : allocation failed ! [1]'
        stop
      endif
!      
!
!===================================================================================
! PCM                                                                              |
!===================================================================================
!
!     STEP 1 : compute rhs phi = A_oo g
!     ---------------------------------
!
!     loop over atoms
      do isph = 1,nsph
!      
!       compute SH expansion glm of phi
        call intrhs( isph, phi(:,isph), glm(:,isph) )
!        
      end do
!
!     P. Gatto, Dec 2016 : why is it not initialized with philm ???
!
!     initial residual : R^0  = g - A_eps * x^0 = g
      vold(:,:) = glm(:,:)
!      
!     $omp parallel do default(shared) private(isph,xlm,x,basloc,vplm,vcos,vsin)
!      
!     initialize
      philm(:,:) = zero ; x(:) = zero ; xlm(:) = zero ; basloc(:) = zero
      vplm(:) = zero ; vcos(:) = zero ; vsin(:) = zero
!      
!
!===================================================================================
!     FMM                                                                          |
!===================================================================================
!
      if (use_fmm) then
!
!       sources for FMM
!       ---------------
        xs(:) = csph(1,:)
        ys(:) = csph(2,:)
        zs(:) = csph(3,:)
        ns = nsph
!
!       targets for FMM
!       ---------------
        nt = 0
        do isph = 1,nsph
          do its = 1,ngrid
!
!           non-zero contribution from target point
            if ( ui(its,isph).gt.zero ) then
!
!             increment number of target points
              nt = nt + 1
!
!             store target point
              xt(nt) = csph(1,isph) + rsph(isph)*grid(1,its)
              yt(nt) = csph(2,isph) + rsph(isph)*grid(2,its)
              zt(nt) = csph(3,isph) + rsph(isph)*grid(3,its)
!
            endif
!            
          enddo
        enddo
!        
!                      m'  4 pi l'     l'+1        m'
!       build [ Phi_j ]  = -------  r_j     [ g_j ]
!                      l'  2l' + 1                 l'
!       ---------------------------------------------
        do jsph = 1,nsph         
!
!         initialize r_j^(l+1) 
          tt = one
!
          do l = 0,lmax
!
            ind = l*l + l + 1
!
!           update r_j^(l+1)
            tt = tt*rsph(jsph)

            do m = -l,l

              phi_j(ind+m,jsph) = tt * glm(ind+m,jsph) / facl(ind+m) * dble(l)

            enddo
          enddo
        enddo
!
!       call to FMM
!       -----------
        call fmm( lmax, ns, phi_j, xs, ys, zs, nt, xt(1:nt), yt(1:nt), zt(1:nt), fmm_vec(1:nt) )
!
!       expand glm at integration points
!       --------------------------------
        ggrid(:,:) = zero ; basloc_fmm(:,:) = zero
        do its=1,ngrid
!
          call ylmbas( grid(:,its), basloc_fmm(:,its), vplm, vcos, vsin )
!
          do isph = 1,nsph
            do l = 0,lmax
!
              ind = l*l + l + 1
!
              do m = -l,l
!
                ggrid(its,isph) = ggrid(its,isph) + basloc_fmm(ind+m,its)*glm(ind+m,isph)
!
              enddo
            enddo
          enddo
        enddo
!
!       explode result of fmm and multiply by U_i^n
!       -------------------------------------------
        xx(:,:) = zero
        nt = 0 
        do isph = 1,nsph
          do its = 1,ngrid
!
!           non-zero contribution from target point
            if ( ui(its,isph).gt.zero ) then
!
!             advance target point index
              nt = nt + 1
!              
!             retrive and multiply
              xx(its,isph) = ui(its,isph)*( fmm_vec(nt) - 2.d0*pi*ggrid(its,isph) )
!              
            endif
          enddo
        enddo
!
!       integrate against SH, add action of identity term
!       -------------------------------------------------
        do isph = 1,nsph

          call intrhs( isph, xx(:,isph), philm_fmm(:,isph) )
          philm_fmm(:,isph) = 2*pi * glm(:,isph) - philm_fmm(:,isph)

        enddo
!
!       redirect
        philm = philm_fmm
!
!
!===================================================================================
!     NON FMM                                                                      |
!===================================================================================
!
      else
!              
!       loop over atoms
        do isph = 1,nsph
!       
!         phi = A_oo g
          call mkrvec(     isph, zero, glm, philm(    :,isph), xlm, x, basloc, vplm, vcos, vsin )
!!!       call mkrvec_fmm( isph, zero, glm, philm_fmm(:,isph), xlm, x, basloc, vplm, vcos, vsin )
!          
        enddo
!
!!!        s1 = zero ; s2 = zero ; s3 = zero
!!!        do isph = 1,nsph
!!!          do j = 1,nbasis
!!!!
!!!            s1 = s1 + ( philm(    j,isph) - philm_fmm(j,isph) )**2
!!!            s2 = s2 + ( philm(    j,isph)                     )**2
!!!            s3 = s3 + ( philm_fmm(j,isph)                     )**2
!!!
!!!          enddo
!!!        enddo
!!!!
!!!        write(*,1002) sqrt(s2), sqrt(s3), sqrt(s1/s2)
!!! 1002   format(' FMM : norm,norm_fmm,error = ',3(e12.5,2x))     
!
      endif
!      
!===================================================================================
!
!
!     $omp parallel do default(shared) private(isph)
!      
!      
!     STEP 2 : solve A_eps W = phi
!     -----------------------------
!    
!fl NEW
      call pcm(.false., .false., .true., xx, philm, wlm_new)
!!!      call prtsph('solution to the pcm equations - new:', nsph, 0, wlm_new)
!     initialize
!fl
      allocate (prec(nbasis,nbasis,nsph), precm1(nbasis,nbasis,nsph) , stat=istatus)
      if ( istatus.ne.0 ) then
        write(*,*)'iefpcm : allocation failed ! [2]'
        stop
      endif
      prec(:,:,:) = zero ; precm1(:,:,:) = zero
!      
!     loop over atoms
      do isph = 1,nsph
!
!       compute inverse of diagonal block A_eps,ii
        call mkprec( isph, .true., prec(:,:,isph), precm1(:,:,isph) )
!        
      end do
!
      wlm = zero
!                                   n
!     STEP 2.1 : Jacobi method for W
!     -------------------------------
!
!     main loop
      if ( iprint.gt.1 ) then
        write(iout,1000)
        write(iout,1010)
 1000   format('   first loop: computing V(eps)')
 1010   format('   it        error        err-00')
      endif
!      
!     Jacobi iteration
      do it = 1,nitmax
!      
!       initialize residual to zero
        wlm(:,:) = zero
!
!
    !!!$omp parallel do default(shared) private(isph,xlm,x,basloc,vplm,vcos,vsin)
!
!                                    n
!       STEP 2.2 : compute residual R
!       ------------------------------
!
!
!===================================================================================
!       FMM                                                                        |
!===================================================================================
!
        if (use_fmm) then
!
!                        m'  4 pi l'     l'+1           m'
!         build [ Phi_j ]  = -------  r_j     [ vold_j ]
!                        l'  2l' + 1                    l'
!         ------------------------------------------------
          do jsph = 1,nsph         
!
!           initialize r_j^(l+1) 
            tt = one
!
            do l = 0,lmax
!
              ind = l*l + l + 1
!
!             update r_j^(l+1)
              tt = tt*rsph(jsph)

              do m = -l,l

                phi_j(ind+m,jsph) = tt * vold(ind+m,jsph) / facl(ind+m) * dble(l)

              enddo
            enddo
          enddo
!
!         call to FMM
!         -----------
          call fmm( lmax, ns, phi_j, xs, ys, zs, nt, xt(1:nt), yt(1:nt), zt(1:nt), fmm_vec(1:nt) )
!          
!         expand vold at integration points
!         ---------------------------------
          voldgrid(:,:) = zero
          do its=1,ngrid

            call ylmbas( grid(:,its), basloc_fmm(:,its), vplm, vcos, vsin )

            do isph = 1,nsph
!            
              do l = 0,lmax
!
                ind = l*l + l + 1
!
                do m = -l,l
!
                  voldgrid(its,isph) = voldgrid(its,isph) + basloc_fmm(ind+m,its)*vold(ind+m,isph)
!
                enddo
              enddo
            enddo
          enddo
!          
!         explode result of fmm and multiply by U_i^n
!         -------------------------------------------
          xx(:,:) = zero
          nt = 0 
          do isph = 1,nsph
            do its = 1,ngrid
!
!             non-zero contribution from target point
              if ( ui(its,isph).gt.zero ) then
!
!               advance target point index
                nt = nt + 1
!                
!               retrive and multiply
                xx(its,isph) = ui(its,isph)*( fmm_vec(nt) - 2.d0*pi*voldgrid(its,isph) )
!                
              endif
            enddo
          enddo
!
!         integrate against SH, add action of identity term
!         -------------------------------------------------
          fep = two*pi*(eps+one)/(eps-one)
          if ( eps.eq.zero )  fep = two*pi
!
          do isph = 1,nsph

            call intrhs( isph, xx(:,isph), wlm_fmm(:,isph) )
            wlm_fmm(:,isph) = fep * vold(:,isph) - wlm_fmm(:,isph)

          enddo
!
!         redirect
          wlm = wlm_fmm
!
!
!===================================================================================
!       NON FMM                                                                    |
!===================================================================================
!
        else
!
!         loop over atoms
          do isph = 1,nsph
!
            call mkrvec(     isph, eps, vold, wlm(:,isph), xlm, x, basloc, vplm, vcos, vsin )
!!!            call mkrvec_fmm( isph, eps, vold, wlm_fmm(:,isph), xlm, x, basloc, vplm, vcos, vsin )
!          
          enddo
!
!!!          s1 = zero ; s2 = zero ; s3 = zero
!!!          do isph = 1,nsph
!!!            do j = 1,nbasis
!!!!
!!!              s1 = s1 + ( wlm(    j,isph) - wlm_fmm(j,isph) )**2
!!!              s2 = s2 + ( wlm(    j,isph)                   )**2
!!!              s3 = s3 + ( wlm_fmm(j,isph)                   )**2
!!!
!!!            enddo
!!!          enddo
!!!!
!!!          write(*,*) 'it = ',it
!!!          write(*,1002) sqrt(s2), sqrt(s3), sqrt(s1/s2)
!        
        endif
!        
!===================================================================================
!
!
!       R^n = Phi - A_eps * R^n-1 
        wlm(:,:) = philm(:,:) - wlm(:,:)
!
!                                n          
!       STEP 2.3 : solve for  W_i  = A_eps,ii^-1 ( ... )
!       ------------------------------------------------
!
!       loop over atoms
        do isph = 1,nsph
!
!         apply inverse
          call DGEMV( 'N', nbasis, nbasis, one, precm1(:,:,isph), nbasis, wlm(:,isph), 1, zero, xlm, 1 )
!
!         update : W_i^n = Phi_i - A_eps,ii W^n-1 
          wlm(:,isph) = vold(:,isph) + xlm(:)
!          
        end do
!
!               n-1    n
!       update W    = W
        vold = vold - wlm
!
!
!       STEP 3.3 : check for convergence
!       --------------------------------
!
        err(:) = zero
!
!       loop over atoms
        do isph = 1,nsph
!        
!         accumulate
          err(:) = err(:) + vold(:,isph)**2
!          
        end do
!
        err = sqrt(err/dble(nsph))
!
!
!!!!===================================================================
!!!!       OLD VERSION
!!!!       compute rms- and max-norm of v_old
!!!        call rmsvec( nbasis*nsph, vold, vrms, vmax )
!!!!===================================================================
!
!
!!!!===================================================================
!!!!       OLD VERSION
!!!        if ( vrms.le.tredis .and. ndiis.gt.0 )  dodiis = .true.
!!!!===================================================================
!
!
!===================================================================
!       NEW VERSION
        if ( ndiis.gt.0 )  dodiis = .true.
!===================================================================
!
!
!       diis extrapolation
        if ( dodiis ) then
!                
          xdiis(:,:,nmat) = wlm
          ediis(:,:,nmat) = vold
!          
          call diis( nbasis*nsph, nmat, xdiis, ediis, bmat, wlm )
!          
        end if
!        
!        
!===================================================================
!       NEW VERSION
!!!        vold = vold - wlm
        vrmsold = vrms
        call rmsvec( nbasis*nsph, vold, vrms, vmax )
!===================================================================
!
!        
        if (iprint.gt.1) then
          write(iout,1020) it, vrms, err(1)
 1020     format(1x,i4,f14.8,121d12.4)
        endif
!
!       convergence has been achieved
        if ( vrms.lt.tol )  goto 999
!        
!       update
        vold = wlm
!        
      end do
!
!     convergence has not been reached
      stop ' convergence failure!!!'
!
  999 continue
!  
!!!      call prtsph('solution to the PCM equations - old:', nsph, 0, wlm)


!     overwrite with new
      wlm = wlm_new

!
!     compute charge distribution and energy
!
!
!!!!     initialize
!!!      sigma_g = zero
!!!!
!!!!     solve  L sigma = W , compute energy
!!!      call itsolv2( .false., .true., wlm, psi, sigma_g, ene )
!!!!
!!!!     save phi_eps for computing forces
!!!      phi_eps = wlm
!
!
!     check solution
!     --------------
!
      if (iprint.gt.1) then
!              
        s1 = zero ; s2 = zero ; s3 = zero ; Awlm = zero
!
        do isph = 1,nsph
!       
!         action of ( A_eps )_ij 
          call mkrvec( isph, eps, wlm, Awlm, xlm, x, basloc, vplm, vcos, vsin )
!        
          do j=1,nbasis
!          
            s1 = s1 + ( philm(j,isph) - Awlm(j) )**2
            s2 = s2 + ( philm(j,isph)           )**2
            s3 = s3 + ( wlm(  j,isph)           )**2
!            
          enddo
        enddo
!
        write(*,*) 'ddPCM : '
        write(*,*) ' '
        write(*,1001) sqrt(s3) , sqrt(s1) / sqrt(s2)
!!!        write(*,1002) sqrt(s2)
 1001   format(' || Phi_eps || , || A_oo Phi - A_eps Phi_eps || / || A_oo Phi || =  ',2(e12.5,2x) )     
 1002   format(' || A_oo Phi || = ',e12.5)
!
        if ( abs(sqrt(s1)/sqrt(s2)) .gt. 1.E-06 ) then
          write(*,*) 'solution failed!'
          write(*,*)'PAUSE - type any key to continue'
          read(*,*)
        endif        
!      
      endif        
!
      if (iprint.gt.1) then
      write(iout,2000)
 2000 format('   first loop has converged.',/,'   second loop: solving ddCOSMO equations for V(eps)')
      endif
!
!
!===================================================================================
! ddCOSMO                                                                          |
!===================================================================================
!
!     initialize
      sigma_g = zero
!
!     solve  L sigma = W , compute energy
      call cosmo(.false., .false., phi, wlm, psi, sigma_g, ene)
      esolv = ene
!     call itsolv2( .false., .true., wlm, psi, sigma_g, ene )
!
!     save phi_eps for computing forces
      phi_eps = wlm
!
!
      deallocate(prec,precm1 , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'iefpcm : deallocation failed ! [1]'
        stop
      endif
!
!     free the memory
      deallocate( err, xlm, x, vold, philm, wlm, wlm_new, glm, basloc, vplm, vcos, vsin, &
                  xdiis, ediis, bmat, stat=istatus )
!                 xdiis, ediis, bmat, prec, precm1 , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'iefpcm : deallocation failed ! [2]'
        stop
      endif
!      
      return
!
!
endsubroutine iefpcm
!----------------------------------------------------------------------------------------
!
!
!
!
!
!-------------------------------------------------------------------------------
! Purpose : solve PCM equation A_eps W = Phi
!-------------------------------------------------------------------------------
subroutine ADJpcm( philm, wlm )
!
      use  ddcosmo
!      
      implicit none
!
      real*8, dimension(nbasis,nsph), intent(in ) :: philm
      real*8, dimension(nbasis,nsph), intent(out) :: wlm
!      
!     local arrays:
      real*8, allocatable :: vold(:,:)
!
!     scratch arrays:
      real*8, allocatable :: basloc(:), vplm(:), vcos(:), vsin(:), xlm(:), x(:)
!
!     diis arrays:
      real*8, allocatable :: xdiis(:,:,:), ediis(:,:,:), bmat(:)
!
      integer :: it, isph, nmat, lenb, istatus
      real*8  :: ene, vrms, vmax, tol,vrmsold, xx(1)
      real*8, allocatable :: err(:), ddiag(:)
      logical :: dodiis
!
      integer, parameter :: nitmax=300
      real*8,  parameter :: ten=10.0d0, tredis=1.d-2
!
!-------------------------------------------------------------------------------
!
!     set solver tolerance
      tol = ten**(-iconv)
!      
!     set up DIIS iterative solver
      dodiis = .false.
      nmat   = 1
      lenb   = ndiis + 1
!
!     allocate memory
      allocate( err(   nbasis)          , &
                xlm(   nbasis)          , &
                x(     ngrid)           , &
                vold(  nbasis,nsph)     , &
                basloc(nbasis)          , &  ! spherical harmonics (SH)
                vplm(  nbasis)          , &  ! service array for SH
                vcos(  lmax+1)          , &  ! service array for SH
                vsin(  lmax+1)          , &  ! service array for SH
                xdiis(nbasis,nsph,ndiis), &
                ediis(nbasis,nsph,ndiis), &
                bmat(lenb*lenb)         , &
                prec(nbasis,nbasis,nsph), precm1(nbasis,nbasis,nsph) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'ADJpcm : allocation failed !'
        stop
      endif
!
!     initial residual : R^0  = phi - A_eps * x^0 = phi
      vold(:,:) = philm(:,:)
!      
!     STEP 2 : solve A_eps W = phi
!     -----------------------------
!    
!     loop over atoms
      do isph = 1,nsph
!
!       compute inverse of diagonal block A_eps^T,ii
        call ADJprec( isph, .true., prec(:,:,isph), precm1(:,:,isph) )
!        
      enddo
!
!
!                                   n
!     STEP 2.1 : Jacobi method for W
!     -------------------------------
!
!     main loop
      if ( .not. iquiet ) then
        write(iout,1000)
        write(iout,1010)
      endif
 1000 format('   first loop: computing V(eps)')
 1010 format('   it        error        err-00')
!      
!     Jacobi iteration
      do it = 1,nitmax
!      
!       initialize residual to zero
        wlm(:,:) = zero
!
      !$omp parallel do default(shared) private(isph,xlm,x,basloc,vplm,vcos,vsin)
!
!                                    n
!       STEP 2.2 : compute residual R
!       ------------------------------
!
!       loop over atoms
        do isph = 1,nsph
!
          call ADJvec( isph, eps, vold, wlm(:,isph), xlm, x, basloc, vplm, vcos, vsin )
!          
        end do
!
!       R^n = Phi - A_eps * R^n-1 
        wlm(:,:) = philm(:,:) - wlm(:,:)
!
!                                n          
!       STEP 2.3 : solve for  W_i  = A_eps,ii^-1 ( ... )
!       ------------------------------------------------
!
!       loop over atoms
        do isph = 1,nsph
!
!         apply inverse
          call DGEMV( 'N', nbasis, nbasis, one, precm1(:,:,isph), nbasis, wlm(:,isph), 1, zero, xlm, 1 )
!
!         update : W_i^n = Phi_i - A_eps,ii W^n-1 
          wlm(:,isph) = vold(:,isph) + xlm(:)
!          
        end do
!
!               n-1    n
!       update W    = W
        vold = vold - wlm
!
!
!       STEP 3.3 : check for convergence
!       --------------------------------
!
        err(:) = zero
!
!       loop over atoms
        do isph = 1,nsph
!        
!         accumulate
          err(:) = err(:) + vold(:,isph)**2
!          
        end do
!
        err = sqrt(err/dble(nsph))
!
!!!!===================================================================
!!!!       OLD VERSION
!!!!       compute rms- and max-norm of v_old
!!!        vrms=zero ; vmax=zero
!!!        call rmsvec( nbasis*nsph, vold, vrms, vmax )
!!!!===================================================================
!
!!!!===================================================================
!!!!       OLD VERSION
!!!        if ( vrms.le.tredis .and. ndiis.gt.0 )  dodiis = .true.
!!!!===================================================================
!
!
!===================================================================
!       NEW VERSION
        if ( ndiis.gt.0 )  dodiis = .true.
!===================================================================
!
!
!       diis extrapolation
        if ( dodiis ) then
!                
          xdiis(:,:,nmat) = wlm
          ediis(:,:,nmat) = vold
!          
          call diis( nbasis*nsph, nmat, xdiis, ediis, bmat, wlm )
!          
        end if
!        
!        
!===================================================================
!       NEW VERSION
        vrmsold = vrms
        call rmsvec( nbasis*nsph, vold, vrms, vmax )
!===================================================================
!
!        
        if ( .not. iquiet ) then
          write(iout,1020) it, vrms, err(1)
        endif
 1020   format(1x,i4,f14.8,121d12.4)
!
!       convergence has been achieved
        if ( vrms.lt.tol )  goto 999
!        
!       update
        vold = wlm
!        
      end do
!
!     convergence has not been reached
      stop ' convergence failure!!!'
!
  999 continue
!  
      if ( .not. iquiet ) then
        write(iout,2000)
      endif
 2000 format('   first loop has converged.',/,'   second loop: solving ddCOSMO equations for V(eps)')
!
!
!
!     free the memory
      deallocate( err, xlm, x, vold, basloc, vplm, vcos, vsin, &
                  xdiis, ediis, bmat, prec, precm1 , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*)'ADJpcm : deallocation failed !'
        stop
      endif
!
!
endsubroutine ADJpcm
!----------------------------------------------------------------------------------------
!
!
!
!
!----------------------------------------------------------------------------------------
! Purpose : compute 
!
!   sum  A_ij v_j
!    j
!
! where A is the PCM matrix, through the following steps . 
!
! 1. Off-diagonal terms :
!
!                        4\pi l           l+1
!   x(n) =   sum    sum  ------ ( t_n^ij )    Y_l^m( s_n^ij ) [ v_j ]_l^m
!          j \ne i  l,m  2l + 1
!
! 2. Add in diagonal term :
!
!                         2\pi
!   x(n) = x(n) + sum  - ------ Y_l^m( s_n^ij ) [ v_i ]_l^m
!                 l,m    2l + 1
!
! 3. Integrate :
!
!   [ dv ]_l^m = sum  w_n U_n^i Y_l^m( s_n ) x(n)
!                 n
!
! 4. Add in identity term :
!
!                     \eps+1
!   [ dv ]_l^m = 2\pi ------ [ v_i ]_l^m - [ dv ]_l^m
!                     \eps-1
!
!
! Remark : when eps_s=0, the eps=oo case is triggered.
!----------------------------------------------------------------------------------------
!
subroutine mkrvec( isph, eps_s, vlm, dvlm, xlm, x, basloc, vplm, vcos, vsin )
!
      use  ddcosmo , only : nbasis, nsph, ngrid, lmax, csph, rsph, grid, basis, ui, facl, &
                            one, pi, zero, pt5, two, dtslm, dtslm2, intrhs, ylmbas, ext1, &
                            do_diag
!      
      implicit none
      integer,                         intent(in   ) :: isph
      real*8,                          intent(in   ) :: eps_s
      real*8,  dimension(nbasis,nsph), intent(in   ) :: vlm
      real*8,  dimension(nbasis),      intent(inout) :: dvlm
      real*8,  dimension(ngrid),       intent(inout) :: x
      real*8,  dimension(nbasis),      intent(inout) :: xlm, basloc, vplm
      real*8,  dimension(lmax+1),      intent(inout) :: vcos, vsin
!
      integer :: its, jsph
      real*8  :: vij(3), sij(3), fep
      real*8  :: vvij, tij, stslm, stslm2, stslm3
!
!----------------------------------------------------------------------------------------
!
!     compute multiplicative coefficient
      fep = two*pi*(eps_s+one)/(eps_s-one)
      if ( eps_s.eq.zero )  fep = two*pi
!
!     initialize
      x(:) = zero
!
!     loop over integration points
      do its = 1,ngrid
! 
!       positive contribution
        if ( ui(its,isph).gt.zero ) then
!
!         loop over spheres
          do jsph = 1,nsph
!
!
!           action of A_ij
!           ==============
!
            if ( jsph.ne.isph ) then
!
!             compute t_n^ij = | r_i + \rho_i s_n - r_j | / \rho_j
              vij  = csph(:,isph) + rsph(isph)*grid(:,its) - csph(:,jsph)
              vvij = sqrt(dot_product(vij,vij))
              tij  = vvij/rsph(jsph) 
!
!             compute s_n^ij = ( r_i + \rho_i s_n - r_j ) / | ... |
              sij  = vij/vvij
!              
!             compute Y_l^m( s_n^ij )
              call ylmbas( sij, basloc, vplm, vcos, vsin )
!              
!             point is INSIDE j-sphere [ EXTENSION ]
!             --------------------------------------
              if ( tij.lt.one ) then
!
!               STEP 1 : x(n) = x(n) + ...
!
!               extension of potential
                select case(ext1)
!
!               t^l extension
                case(0)
                x(its) = x(its) + dtslm2( tij, vlm(:,jsph), basloc )
!
!               constant extension
                case(1)
                x(its) = x(its) + dtslm2( one, vlm(:,jsph), basloc )

!               t^-l extension
                case(2)
                x(its) = x(its) + dtslm(  tij, vlm(:,jsph), basloc )
!
               endselect
!                      
!             point is OUTSIDE j-sphere
!             -------------------------
              else
!                      
!               STEP 1 : x(n) = x(n) + ... 
                x(its) = x(its) + dtslm2( tij, vlm(:,jsph), basloc )
!                      
              endif
!                      
!                      
!           action of A_ii [ excluding identity term ]
!           ==============
!                      
            elseif ( do_diag ) then
!
!             compute x_l^m = - 2pi / (2l+1) [ v_i ]_l^m
              xlm(:) = -pt5 * vlm(:,isph) / facl(:)
!
!             STEP 2 : x(n) = x(n) + sum_l,m  Y_l^m( s_n ) x_l^m
              x(its) = x(its) + dot_product( basis(:,its), xlm(:) )
!              
!              
            endif
          enddo
        endif
      enddo
!
!     STEP 2
      x(:) = x(:) * ui(:,isph)
!
!     STEP 3 : integrate against SH
      call intrhs( isph, x, dvlm )
!
!
!     STEP 4 : add action of identity term
      if ( do_diag ) then
!              
        dvlm(:) = fep * vlm(:,isph) - dvlm(:)
!        
      else
!              
        dvlm(:) = -dvlm(:)
!        
      endif
!
!
endsubroutine mkrvec
!-------------------------------------------------------------------------------
!
!
!
!----------------------------------------------------------------------------------------
! Purpose : action of A_i:^eps , i.e., sum_j A_ij^eps v_j .
!
! Let's drop ^eps for semplicity. Then :
!
!   dvlm_i =   sum   A_ij vl'm'_j + A_ii vl'm'_i =
!            j \ne i
!
!                            4pi l'                              l'+1
!          = -  sum     sum  ------ sum w_n  Y_l^m(s_n)  U_i^n  t      Y_l'^m'(s_ijn)  vl'm'_j
!             j \ne i  l',m' 2l'+1   n
!
!                  eps+1                2pi
!            + 2pi ----- vlm_i +  sum  -----  sum w_n  Y_l^m(s_n)  U_i^n  Y_l'^m'(s_n)  vl'm'_i
!                  eps-1         l',m' 2l'+1   n      
!
!                                                         4pi l'   l'+1
!          = - sum  w_n  Y_l^m(s_n)  U_i^n   sum     sum  ------  t      Y_l'^m'(s_ijn)  vl'm'_j
!               n                          j \ne i  l',m' 2l'+1
! 
!                                    |------------------------- x1(n) -------------------------|
!
!                                                     2pi
!            - sum  w_n  Y_l^m(s_n)  U_i^n   sum   - -----  Y_l'^m'(s_n)  vl'm'_i
!                n                          l',m'    2l'+1
!
!                                    |----------------- x2(n) ------------------|
!
!                  eps+1         
!            + 2pi ----- vlm_i 
!                  eps-1        
!
!                                                            eps+1         
!          = - sum  w_n  Y_l^m(s_n)  ( x1(n) + x2(n) ) + 2pi ----- vlm_i
!               n                                            eps-1        
!
! Remark : when eps_s=0, the eps=oo case is triggered.
!----------------------------------------------------------------------------------------
!
subroutine mkrvec_fmm( isph, eps_s, vlm, dvlm, xlm, x, basloc, vplm, vcos, vsin )
!
      use  ddcosmo , only : nbasis, nsph, ngrid, lmax, csph, rsph, grid, basis, ui, facl, &
                            one, pi, zero, pt5, two, dtslm, dtslm2, intrhs, ylmbas, ext1
!      
      implicit none
      integer,                         intent(in   ) :: isph
      real*8,                          intent(in   ) :: eps_s
      real*8,  dimension(nbasis,nsph), intent(in   ) :: vlm
      real*8,  dimension(nbasis),      intent(inout) :: dvlm
      real*8,  dimension(ngrid),       intent(inout) :: x
      real*8,  dimension(nbasis),      intent(inout) :: xlm, basloc, vplm
      real*8,  dimension(lmax+1),      intent(inout) :: vcos, vsin
!
      integer :: its, jsph, m, l, ns, nt, ind
      real*8  :: fep, xs(nsph),ys(nsph),zs(nsph),xt(ngrid),yt(ngrid),zt(ngrid)
      real*8  :: phi(nbasis,nsph),tt,fmm_vec(ngrid)

      real*8 :: vgrid(ngrid)
!
!----------------------------------------------------------------------------------------
!
!     compute f( \eps )
      fep = two*pi*(eps_s+one)/(eps_s-one)
      if ( eps_s.eq.zero )  fep = two*pi
!
!     initialize
      x(:) = zero
!
!     sources for FMM
      xs(:) = csph(1,:)
      ys(:) = csph(2,:)
      zs(:) = csph(3,:)
      ns = nsph
!
!     targets for FMM
      nt = 0
      do its = 1,ngrid
!
!       non-zero contribution from target point
        if ( ui(its,isph).gt.zero ) then
!
!         increment number of target points
          nt = nt + 1
!
!         store target point
          xt(nt) = csph(1,isph) + rsph(isph)*grid(1,its)
          yt(nt) = csph(2,isph) + rsph(isph)*grid(2,its)
          zt(nt) = csph(3,isph) + rsph(isph)*grid(3,its)
!
        endif
      enddo
!      
!                    m'  4 pi l'     l'+1        m'
!     build [ Phi_j ]  = -------  r_j     [ v_j ]
!                    l'  2l' + 1                 l'
!
      do jsph = 1,nsph         
!
!       initialize r_j^(l+1) 
        tt = one
!
        do l = 0,lmax
!
          ind = l*l + l + 1
!
!         update r_j^(l+1)
          tt = tt*rsph(jsph)

          do m = -l,l

            phi(ind+m,jsph) = tt * vlm(ind+m,jsph) / facl(ind+m) * dble(l)

          enddo
        enddo
      enddo
!
!     call to FMM
      fmm_vec = 0.d0
      call fmm( lmax, ns, phi, xs, ys, zs, nt, xt(1:nt), yt(1:nt), zt(1:nt), fmm_vec(1:nt) )
!
      vgrid(:) = zero
      do its=1,ngrid

        call ylmbas( grid(:,its), basloc, vplm, vcos, vsin )

        do l = 0,lmax
!
          ind = l*l + l + 1
!
          do m = -l,l

            vgrid(its) = vgrid(its) + basloc(ind+m)*vlm(ind+m,isph)

          enddo
        enddo
      enddo
!
!     expand result of fmm and multiply by U_i^n
      nt = 0
      do its = 1,ngrid
!
!       non-zero contribution from target point
        if ( ui(its,isph).gt.zero ) then
!
!         advance target point index
          nt = nt + 1
!          
!         retrive and multiply
          x(its) = ui(its,isph)*( fmm_vec(nt) - 2.d0*pi*vgrid(its) )
!          
        endif
      enddo
!
!     integrate against SH
      call intrhs( isph, x, dvlm )
!
!
!     add action of identity term
!     ===========================
!
      dvlm(:) =  fep * vlm(:,isph) - dvlm(:)
!
!
endsubroutine mkrvec_fmm
!-------------------------------------------------------------------------------
!
!
!
!
!
!-------------------------------------------------------------------------------
! Purpose : compute
!
!   sum (A^T)_ij v_j
!    j
!
! where A is the PCM matrix, through the following steps .
!
! 1. Compute :
!
!   ss = sum  Y_l^m( s_n ) [ v_j ]_l^m
!        l,m                       
!
! 2. Compute :
!
!           4\pi l                         -(l+1)
!   f_l^m = ------ Y_l^m( s_n^ji ) (t_n^ji)      
!           2l + 1
!
! 3. Off-diagonal terms :
!
!  x_l^m = x_l^m - U_n^j * f_l^m * ss
!
! 4. Compute :
!
!   ss = sum  Y_l^m( s_n ) [ v_i ]_l^m
!        l,m                       
!
! 5. Add in diagonal term :
!
!                    2\pi
!   x_l^m = x_l^m - ------ U_n^i * Y_l^m( s_n ) * ss
!                   2l + 1
!
! 6. Accumulate over integration points :
!
!  [ dv ]_l^m = [ dv ]_l^m + w_n * x_l^m
!
! 7. Add in identity term :
!
!                    \eps+1
!  [ dv ]_l^m = 2\pi ------ [ v_i ]_l^m - [ dv ]_l^m
!                    \eps-1
!
!-------------------------------------------------------------------------------
!
subroutine ADJvec( isph, eps_s, vlm, dvlm, xlm, x, basloc, vplm, vcos, vsin )
!
      use  ddcosmo , only : nbasis, nsph, ngrid, lmax, csph, rsph, grid, basis, &
                            ui, facl, two, pi, one, zero, pt5, w, ylmbas, ext1, &
                            do_diag, tylm
!      
      implicit none
      integer,                         intent(in   ) :: isph
      real*8,                          intent(in   ) :: eps_s
      real*8,  dimension(nbasis,nsph), intent(in   ) :: vlm
      real*8,  dimension(nbasis),      intent(inout) :: dvlm
      real*8,  dimension(ngrid),       intent(inout) :: x
      real*8,  dimension(nbasis),      intent(inout) :: xlm, basloc, vplm
      real*8,  dimension(lmax+1),      intent(inout) :: vcos, vsin
!
      real*8  :: vji(3), sji(3), fep, vvji, tji, tt, ss, flm(nbasis)
      integer :: n, jsph
!
!-------------------------------------------------------------------------------
!
!     compute multiplicative coefficient
      fep = two*pi*(eps_s+one)/(eps_s-one)
      if ( eps_s.eq.zero )  fep = two*pi
!
!     initialize [ UNUSED !!! ] workspace 
      x = zero
!
!     initialize
      dvlm = zero
!
!     loop over integration points
      do n = 1,ngrid
!
!       initialize
        xlm = zero

!       loop over spheres
        do jsph = 1,nsph
!          
!         non-null contribution from integration point 
          if ( ui(n,jsph).gt.zero ) then
!
!
!           action of (A^T)_ij
!           ==================
!
            if ( jsph.ne.isph ) then
!
!             compute t_n^ji
              vji  = csph(:,jsph) + rsph(jsph)*grid(:,n) - csph(:,isph)
              vvji = sqrt( dot_product( vji,vji ) )
              tji  =  vvji/rsph(isph)
!              
!             compute s_n^ji
              sji = vji/vvji
!              
!             contract over l, m
              ss = dot_product( basis(:,n), vlm(:,jsph) )
!              
!             compute Y_l^m( s_n^ji )
              call ylmbas( sji, basloc, vplm, vcos, vsin )
!              
!             point is INSIDE i-sphere [ EXTENSION ]
!             --------------------------------------
              if ( tji.lt.one ) then
!
!               extension of potential
                select case(ext1)
!
!               t^-l extension
                case(0)
                call tylm( tji, basloc, flm )
!                
!               constant extension
                case(1)
                call tylm( one, basloc, flm )
!
!               (1/t)^-l extension
                case(2)
                call tylm( one/tji, basloc, flm )
!
                endselect
!                
!
!             point is OUTSIDE i-sphere
!             -------------------------
              else 
!
                call tylm( tji, basloc, flm )
!                
              endif
!
!             accumulate over j
              xlm(:) = xlm(:) + ui(n,jsph) * flm(:) * ss
!                      
!                      
!           action of (A^T)_ii [ excluding identity term ]
!           ==================
!                      
            elseif ( do_diag ) then
!
!             contract over l, m
              ss = dot_product( basis(:,n), vlm(:,isph) )
!              
!             accumulate x_l^m = x_l^m - U_i^n * Y_l^m(s_n) * ss * 2\pi / (2l + 1)
              xlm(:) = xlm(:) - ui(n,isph) * basis(:,n) * ss * pt5/facl(:)
!              
            endif
          endif
        enddo
!        
!       accumulate over n
        dvlm(:) = dvlm(:) + w(n) * xlm(:)
!
      enddo
!
!
!     add action of identity term
!     ===========================
!
      if ( do_diag ) then 
!              
        dvlm(:) = fep*vlm(:,isph) - dvlm(:)
!        
      else
!              
        dvlm(:) = - dvlm(:)
!        
      endif
!
!
endsubroutine ADJvec
!-------------------------------------------------------------------------------
!
!
!
!
!
!
!-------------------------------------------------------------------------------
! Purpose : compute A_eps,ii^-1.
!
! Remark : this is a bit more than a preconditioner... in fact, it is the
!          exact (numerical) inverse. 
!
!-------------------------------------------------------------------------------
!
subroutine mkprec( isph, doinv, p, pm1 )
!
      use  ddcosmo
!     use  dd_utils , only : nbasis, ngrid, basis, lmax, ui, w, eps
!      
      implicit none
      integer,                          intent(in   ) :: isph
      logical,                          intent(in   ) :: doinv
      real*8, dimension(nbasis,nbasis), intent(inout) :: p
      real*8, dimension(nbasis,nbasis), intent(inout) :: pm1
!
      integer :: l, m, ind, l1, m1, ind1, info, its, istatus
      real*8  :: f, f1
!
      integer, allocatable :: ipiv(:)
      real*8,  allocatable :: work(:)
!
!-------------------------------------------------------------------------------
!
      allocate( ipiv(nbasis) , work(nbasis*nbasis) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'mkprec : allocation failed !'
        stop
      endif
!
!     initialize
      p(:,:) = zero
!      
!      
!     STEP 1
!     ------
!
!     loop over grid points
      do its = 1,ngrid
!      
!       1st loop over SH degree
        do l = 0,lmax
!
!         index associated to Y_l^0
          ind = l*l + l + 1
!
!         1st loop over SH order
          do m = -l,l
!          
            f = w(its) * basis(ind+m,its) * ui(its,isph)
!            
!           2nd loop over SH degree
            do l1 = 0,lmax
!
!             index associated to Y_l1^0
              ind1 = l1*l1 + l1 + 1
!
!             2nd loop over SH order
              do m1 = -l1, l1
!              
                f1 = two * pi / (two*dble(l1) + one) * basis(ind1+m1,its)
!                
!               accumulate
                p(ind+m,ind1+m1) = p(ind+m,ind1+m1) + f * f1
!                
              end do
            end do
          end do
        end do
      end do
!
!
!     STEP 2 : diagonal part
!     ----------------------
!
!     lopp over SH degree
      do l = 0,lmax
!
!       index associated to Y_l^0
        ind = l*l + l + 1
!
!       loop over SH order
        do m = -l,l
!        
          f = two * pi * (eps+one) / (eps-one)
!
!         accumulate
          p(ind+m,ind+m) = p(ind+m,ind+m) + f
!          
        end do
      end do
!      
!
!     STEP 3 : invert preconditioner
!     ------------------------------
!
      if ( doinv ) then
!              
        pm1(:,:) = p(:,:)
!        
!       compute LU factorization
        call DGETRF( nbasis, nbasis, pm1, nbasis, ipiv, info )
!
!       check
        if ( info.ne.0 ) then
          write(6,*) 'mkprec, DGETRF : info = ', info
        end if
!        
!       invert factorization
        call DGETRI( nbasis, pm1, nbasis, ipiv, work, nbasis*nbasis, info )
!        
!       check
        if ( info.ne.0 ) then
          write(6,*) 'mkprec, DGETRI : info = ', info
        end if
!
      end if
!      
      deallocate ( work, ipiv, stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'mkprec : deallocation failed !'
        stop
      endif
!      
      return
!      
!      
endsubroutine mkprec
!-------------------------------------------------------------------------------
!
!
!
!
!-------------------------------------------------------------------------------
subroutine ADJprec( isph, doinv, p, pm1 )
!
      use  ddcosmo
!     use  dd_utils , only : nbasis, ngrid, basis, lmax, ui, w, eps
!      
      implicit none
      integer,                          intent(in   ) :: isph
      logical,                          intent(in   ) :: doinv
      real*8, dimension(nbasis,nbasis), intent(inout) :: p
      real*8, dimension(nbasis,nbasis), intent(inout) :: pm1
!
      integer :: l, m, ind, l1, m1, ind1, info, its, istatus, iiprint
      real*8  :: f, f1
!
      integer, allocatable :: ipiv(:)
      real*8,  allocatable :: work(:)
!
!-------------------------------------------------------------------------------
!
      allocate( ipiv(nbasis) , work(nbasis*nbasis) , stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'ADJprec : allocation failed !'
        stop
      endif
!
!     initialize
      p(:,:) = zero
!      
!      
!     STEP 1
!     ------
!
!     loop over grid points
      do its = 1,ngrid
!      
!       1st loop over SH degree
        do l = 0,lmax
!
!         index associated to Y_l^0
          ind = l*l + l + 1
!
!         1st loop over SH order
          do m = -l,l
!          
            f = two*pi/(two*dble(l) + one) * w(its) * basis(ind+m,its) * ui(its,isph)
!            
!           2nd loop over SH degree
            do l1 = 0,lmax
!
!             index associated to Y_l1^0
              ind1 = l1*l1 + l1 + 1
!
!             2nd loop over SH order
              do m1 = -l1, l1
!                
!               accumulate
                p(ind+m,ind1+m1) = p(ind+m,ind1+m1) + f * basis(ind1+m1,its)
!                
              end do
            end do
          end do
        end do
      end do
!
!
!     STEP 2 : diagonal part
!     ----------------------
!
!     loop over SH degree
      do l = 0,lmax
!
!       index associated to Y_l^0
        ind = l*l + l + 1
!
!       loop over SH order
        do m = -l,l
!        
          f = two * pi * (eps+one) / (eps-one)
!
!         accumulate
          p(ind+m,ind+m) = p(ind+m,ind+m) + f
!          
        end do
      end do
!
!     print to screen
      iiprint=0
      if ( iiprint.ne.0 ) then
      write(*,1000) isph
 1000 format(' Preconditioner for isph = ',i2)     
      write(*,*)''
!      
      write(*,*) 'M = '
      do l= 1,nbasis
!      
        write(*,"(4x,300(e12.5,2x))") ( p(l,l1), l1=1,nbasis )
!
      enddo
      write(*,*)''
      endif
!
!      
!
!     STEP 3 : invert preconditioner
!     ------------------------------
!
      if ( doinv ) then
!              
        pm1(:,:) = p(:,:)
!        
!       compute LU factorization
        call DGETRF( nbasis, nbasis, pm1, nbasis, ipiv, info )
!
!       check
        if ( info.ne.0 ) then
          write(6,*) 'mkprec, DGETRF : info = ', info
        end if
!        
!       invert factorization
        call DGETRI( nbasis, pm1, nbasis, ipiv, work, nbasis*nbasis, info )
!        
!       check
        if ( info.ne.0 ) then
          write(6,*) 'mkprec, DGETRI : info = ', info
        end if
!
      end if
!      
      deallocate ( work, ipiv, stat=istatus )
      if ( istatus.ne.0 ) then
        write(*,*) 'mkprec : deallocation failed !'
        stop
      endif
!      
      return
!      
!      
endsubroutine ADJprec

!!!subroutine pause
!!!      implicit none
!!!      write(*,*)'PAUSE - press any key'
!!!      read(*,*)
!!!      return
!!!endsubroutine pause
