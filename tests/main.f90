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
use dd_operators
implicit none

type(dd_data_type) :: dd_data
integer :: iprint, nproc, lmax, pmax, ngrid, iconv, igrad, n, force, fmm, model
integer :: nngmax=200
real(kind=rp) :: eps, eta
real(kind=rp), allocatable :: x(:), y(:), z(:), rvdw(:), charge(:)
real(kind=rp), allocatable :: phi(:), psi(:,:)
real(kind=rp), parameter :: toang=0.52917721092d0, tokcal=627.509469d0
real(kind=rp), parameter :: tobohr=1d0/toang

integer :: i

open (unit=100,file='Input.txt',form='formatted',access='sequential')
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
!
! call the initialization routine. this routine allocates memory, computes some
! quantities for internal use and creates and fills an array ccav(3,ncav) with
! the coordinates of the grid points at which the user needs to compute the potential.
! ncav is the number of external grid points and nbasis the number of spherical
! harmonics functions used for the expansion of the various ddcosmo quantities;
! both are computed by ddinit and defined as common variables in ddcosmo.mod.

model=1
force=0
fmm=0
call ddinit(n, x, y, z, rvdw, model, lmax, ngrid, force, fmm, pmax, pmax, &
    & iprint, nngmax, eta, dd_data)
allocate (phi(dd_data % ncav), psi(dd_data % nbasis,n))
!call mkrhs(n, charge, x, y, z, dd_data % ncav, dd_data % ccav, phi, &
!    & dd_data % nbasis, psi)
deallocate(phi, psi)
deallocate(x, y, z, rvdw, charge)
end program main

