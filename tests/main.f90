!> @copyright (c) 2020-2020 RWTH Aachen. All rights reserved.
!!
!! ddX software
!!
!! @file tests/main.f90
!! Main ddx driver.
!!
!! @version 1.0.0
!! @author Aleksandr Mikhalev
!! @date 2020-12-17

program main
use dd_core
use dd_operators
use dd_solvers
use dd_cosmo
use dd_pcm
implicit none

character(len=255) :: fname
type(dd_data_type) :: dd_data
integer :: iprint, nproc, lmax, pmax, ngrid, iconv, igrad, n, force, fmm, model
integer :: nngmax=200, niter, ndiis=25, info
logical :: ok
real(dp) :: eps, eta, tol, se=zero, kappa
real(dp), allocatable :: x(:), y(:), z(:), rvdw(:), charge(:)
real(dp), allocatable :: phi(:), gradphi(:, :), psi(:, :), xs(:, :)
real(dp), allocatable :: g(:, :), rhs(:, :)
real(dp), parameter :: toang=0.52917721092d0, tokcal=627.509469d0
real(dp), parameter :: tobohr=1d0/toang

integer :: i, j

! Read input file name
call getarg(1, fname)
write(*, *) "Reading input file ", fname
open (unit=100,file=fname,form='formatted',access='sequential')
!
! scalar parameters. the variables are defined in the ddcosmo module and are common to
! all the ddcosmo routines (no need to declare them if ddcosmo.mod is loaded.)
!
read(100,*) iprint      ! printing flag
read(100,*) nproc       ! number of openmp threads
read(100,*) lmax        ! max angular momentum of spherical harmonics basis
read(100,*) pmax        ! max degree of harmonics for the FMM
read(100,*) ngrid       ! number of lebedev points
read(100,*) iconv       ! 10^(-iconv) is the convergence threshold for the iterative solver
read(100,*) igrad       ! whether to compute (1) or not (0) forces
read(100,*) eps         ! dielectric constant of the solvent
read(100,*) eta         ! regularization parameter
!
read(100,*) n           ! number of atoms
!
allocate(x(n), y(n), z(n), rvdw(n), charge(n))
!
! we also read from the same file the charges, coordinates and vdw radii.
! in this example, the coordinates and radii are read in angstrom and
! converted in bohr before calling ddinit.
!
do i = 1, n
  read(100,*) charge(i), x(i), y(i), z(i), rvdw(i)
end do
x    = x*tobohr
y    = y*tobohr
z    = z*tobohr
rvdw = rvdw*tobohr
!
close (100)

model=2
force=0
fmm=0
kappa=0d0
call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pmax, pmax, &
    & iprint, nngmax, se, eta, eps, kappa, dd_data, info)
allocate(phi(dd_data % ncav), gradphi(3, dd_data % ncav), &
    & psi(dd_data % nbasis,n))
call mkrhs(n, charge, x, y, z, dd_data % ncav, dd_data % ccav, phi, gradphi, &
    & dd_data % nbasis, psi)
tol = 10d0 ** (-iconv)
niter = 200
call ddpcm(dd_data, phi, psi, tol, ndiis, niter)
deallocate(phi, gradphi, psi)
deallocate(x, y, z, rvdw, charge)
call ddfree(dd_data)

end program main

