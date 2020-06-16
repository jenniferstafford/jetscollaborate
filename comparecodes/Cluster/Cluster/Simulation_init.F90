!!****if* source/Simulation/SimulationMain/Sedov/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init(integer(IN) :: myPE)
!!
!! DESCRIPTION
!!
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for Sedov Spherical Explosion 
!!  problem.
!!
!! ARGUMENTS
!!
!!   myPE -   current processor number
!!
!! PARAMETERS
!!
!!***

subroutine Simulation_init(myPE)

  use Simulation_data 
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  implicit none

#include "constants.h"
#include "Flash.h"

  integer, intent(in) :: myPE
  integer :: i, ii
  real :: tmp

  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('xmin',sim_xMin)
  call RuntimeParameters_get('xmax',sim_xMax)
  call RuntimeParameters_get('ymin',sim_yMin)
  call RuntimeParameters_get('ymax',sim_yMax)
  call RuntimeParameters_get('zmin',sim_zMin)
  call RuntimeParameters_get('zmax',sim_zMax)
  call RuntimeParameters_get('eos_singleSpeciesA', sim_eos_singleSpeciesA)
  call RuntimeParameters_get('killdivb',sim_killdivb)
  call RuntimeParameters_get('Bfield0', sim_Bfield0)
  call RuntimeParameters_get('saturated', sim_saturatedConduction)
!  call RuntimeParameters_get('T_in',T_in)
!  call RuntimeParameters_get('T_out',T_out)
!  call RuntimeParameters_get('rho_in',rho_in)
!  call RuntimeParameters_get('Mvir',Mvir)
  call RuntimeParameters_get('M0',M0)
  call RuntimeParameters_get('rs',rs)
  call RuntimeParameters_get('rc',rc)

  call PhysicalConstants_get('Boltzmann', kb)
  call PhysicalConstants_get('Newton', newton)
  call PhysicalConstants_get('proton mass', mp)

!  h = 0.7e0
!  Mvir = Mvir/h
!  c  = 6.0* ((Mvir/(1.0e14/h))**(-0.2e0))
!  mc = log(1.0+c) - c/(1.0+c)
!  M0 = (Mvir*1.9889225e33)/(2.0*mc)
!  rvir = 4.95e24/h *(Mvir*h/5.0e14)**(1.0/3.0)
!  rs = rvir/c


  Msun = 1.9889225E33
  Mpc  = 3.0856775807E24
  km   = 1.0e5 
  kpc  = 3.0856775807E21

  rs = rs * kpc
  rc = rc * kpc
  M0 = M0 * Msun


  DeltaT = T_out-T_in


!  open(unit=1,file='galaxies.dat',status='old')
!
!  read(1,*)ntot_tmp
!
!  allocate(Positions_Ascii_tmp(3,ntot_tmp))
!  allocate(Velocities_Ascii_tmp(3,ntot_tmp))
!  allocate(Masses_Ascii_tmp(ntot_tmp))
!
!  do i=1,ntot_tmp
!
!     read(1,*) Positions_Ascii_tmp(1,i), Positions_Ascii_tmp(2,i), Positions_Ascii_tmp(3,i), &
!               Velocities_Ascii_tmp(1,i), Velocities_Ascii_tmp(2,i), Velocities_Ascii_tmp(3,i), &
!               Masses_Ascii_tmp(i)
!
!  enddo
!  close (1)
!
!  do i=1,ntot_tmp
!
!     Positions_Ascii_tmp(1,i)   = Positions_Ascii_tmp(1,i)    * Mpc 
!     Positions_Ascii_tmp(2,i)   = Positions_Ascii_tmp(2,i)    * Mpc
!     Positions_Ascii_tmp(3,i)   = Positions_Ascii_tmp(3,i)    * Mpc
!     Velocities_Ascii_tmp(1,i)  = Velocities_Ascii_tmp(1,i)   * km
!     Velocities_Ascii_tmp(2,i)  = Velocities_Ascii_tmp(2,i)   * km
!     Velocities_Ascii_tmp(3,i)  = Velocities_Ascii_tmp(3,i)   * km
!!    Masses_Ascii_tmp(i)        = Masses_Ascii_tmp(i)         * Msun
!     Masses_Ascii_tmp(i)        = 5.0e11 * Msun
!
!  enddo
!
!
!! filtering out the particles that fall outside the domain
!
!  ii = 0
!  do i=1,ntot_tmp
!     if ( &
!            (abs(Positions_Ascii_tmp(1,i)) .lt. sim_xMax) .or. &
!            (abs(Positions_Ascii_tmp(1,i)) .lt. sim_xMax) .or. &
!            (abs(Positions_Ascii_tmp(1,i)) .lt. sim_xMax)           ) then
!        ii = ii + 1
!     endif
!  enddo
!
!  ntot = ii   !actual number of particles used in the beginning of the simulation (some may escape) in the domain
!
!
!  write(*,*)'ntot = ', ntot
!!  stop
!
!
!  allocate(Positions_Ascii(3,ntot))
!  allocate(Velocities_Ascii(3,ntot))
!  allocate(Masses_Ascii(ntot))
!
!  ii = 0
!  do i=1,ntot_tmp
!     if ( &
!            (abs(Positions_Ascii_tmp(1,i)) .lt. sim_xMax) .or. &
!            (abs(Positions_Ascii_tmp(1,i)) .lt. sim_xMax) .or. &
!            (abs(Positions_Ascii_tmp(1,i)) .lt. sim_xMax)           ) then
!
!         ii = ii + 1
!         Positions_Ascii(1,ii)   = Positions_Ascii_tmp(1,i)
!         Positions_Ascii(2,ii)   = Positions_Ascii_tmp(2,i) 
!         Positions_Ascii(3,ii)   = Positions_Ascii_tmp(3,i)   
!         Velocities_Ascii(1,ii)  = Velocities_Ascii_tmp(1,i)  
!         Velocities_Ascii(2,ii)  = Velocities_Ascii_tmp(2,i)  
!         Velocities_Ascii(3,ii)  = Velocities_Ascii_tmp(3,i)  
!         Masses_Ascii(ii)        = Masses_Ascii_tmp(i)        
!
!     endif
!  enddo
!
!!  write(*,*)'particle data', Positions_Ascii, Velocities_Ascii, Masses_Ascii
!!  stop
!
!  deallocate(Positions_Ascii_tmp)
!  deallocate(Velocities_Ascii_tmp)
!  deallocate(Masses_Ascii_tmp)

!  open(unit=1,file='bx.dat',status='old',form='unformatted')
  open(unit=1,file='bx.dat',status='old')
  read(1,*)ntotal, tbx
  allocate(bx(ntotal,ntotal,ntotal))
  read(1,*) bx
  close (1)

!  open(unit=1,file='by.dat',status='old',form='unformatted')
  open(unit=1,file='by.dat',status='old')
  read(1,*)ntotal, tby
  allocate(by(ntotal,ntotal,ntotal))
  read(1,*) by
  close (1)

!  open(unit=1,file='bz.dat',status='old',form='unformatted')
  open(unit=1,file='bz.dat',status='old')
  read(1,*)ntotal, tbz
  allocate(bz(ntotal,ntotal,ntotal))
  read(1,*) bz
  close (1)

! write(*,*) bx(3,1,1), bx(8,2,120), bx(4,50,4)
! stop

  !reading in the ICs for the gas
  open(unit=1,file='A2199_gas_profile_linear.dat',status='old')
  read(1,*) ngas, dr
  allocate(rhogas(ngas))
  allocate(tempgas(ngas))
  do i = 1, ngas
     read(1,*) tmp, rhogas(i), tempgas(i) 
  enddo
  close (1)

end subroutine Simulation_init

