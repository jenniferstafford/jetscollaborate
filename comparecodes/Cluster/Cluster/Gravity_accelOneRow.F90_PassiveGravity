!!****if* source/physics/Gravity/GravityMain/Poisson/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow
!!
!!
!! SYNOPSIS
!!
!!  Gravity_accelOneRow(integer(2),intent(in):: pos, 
!!                      integer, intent(in) :: sweepDir, 
!!                      integer, intent(in) :: blockID, 
!!                      integer, intent(in) :: numCells, 
!!                      real(numCells),intent(out) :: grav, 
!!                      integer, intent(in),optional :: potentialIndex)
!!                      
!!                      
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration on all mesh blocks.  Either a single component
!!  of the acceleration or all three can be computed.
!!
!!  This routine computes the gravitational acceleration for a row
!!  of zones in a specified direction in a given block. First-order
!!  finite-volume differencing is used everywhere.  Formulae based
!!  on long stencils (usually high-order) may produce differences
!!  at the block boundaries for siblings as hydro solver may require
!!  several valid guard cells (e.g., PPM with parabolic
!!  interpolation for force terms needs 3 valid guard cells). Not
!!  providing such valid data may result in violation of conservation. 
!!
!! ARGUMENTS
!!
!!  pos     -       Row indices transverse to the sweep direction
!!  sweepDir   -       The sweep direction:  test against sweep_x,
!!                                 sweep_y, and sweep_z
!!  blockID   -     The local identifier of the block to work on
!!  grav()   -       Array to receive result
!!  numCells -       Number of cells to update in grav array
!!  potentialIndex      -  if specified,  Variable # to take as potential.
!!                         Default is GPOT_VAR for the potential stored in the
!!                         gpot slot of unk, which should correspond to the
!!                         potential at the current timestep.
!!
!!
!!***


subroutine Gravity_accelOneRow (pos, sweepDir, blockID, numCells, grav, potentialIndex)

  use Simulation_data

  use Grid_interface, ONLY : Grid_getBlkPhysicalSize, Grid_getBlkPtr, &
    Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_getCellCoords
  implicit none

#include "Flash.h"
#include "constants.h"

  integer, dimension(2), intent(in) :: pos
  integer, intent(in)               :: sweepDir, blockID,  numCells
  real, intent(inout)               :: grav(numCells)
  integer, intent(IN),optional        :: potentialIndex
  real            :: blockSize(MDIM)
  real, POINTER, DIMENSION(:,:,:,:) :: solnVec
  
  integer         :: ii, iimin, iimax, lb
  real            :: gpot(numCells), delxinv
  real, parameter :: onesixth = 1.e0/6.e0
  integer         :: potVar

  real :: g, x
  
!==========================================================================

#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_KHI_GC) :: zCenter
  real,dimension(GRID_JHI_GC) :: yCenter
  real,dimension(GRID_IHI_GC) :: xCenter
#else
  real,allocatable,dimension(:) ::xCenter,yCenter,zCenter
#endif
  real :: dr32, tmpdr32

  integer :: sizeX,sizeY,sizez

  integer :: j,k
  logical :: gcell = .true.

!==================================================
  
  
  j=pos(1)
  k=pos(2)
#ifndef FIXEDBLOCKSIZE
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX=blkLimitsGC(HIGH,IAXIS)
  sizeY=blkLimitsGC(HIGH,JAXIS)
  sizeZ=blkLimitsGC(HIGH,KAXIS)
  allocate(xCenter(sizeX))
  allocate(yCenter(sizeY))
  allocate(zCenter(sizeZ))
#else
  sizeX=GRID_IHI_GC
  sizeY=GRID_JHI_GC
  sizeZ=GRID_KHI_GC
#endif
  zCenter = 0.
  yCenter = 0.
  if (NDIM == 3) call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zCenter, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, yCenter, sizeY)
  call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, xCenter, sizeX)

  if (sweepDir .eq. SWEEP_X) then               ! x-component

     tmpdr32 = yCenter(j)*yCenter(j) + zCenter(k)*zCenter(k) 

     do ii = 1, numCells

        dr32 = sqrt(xCenter(ii)*xCenter(ii) + tmpdr32)
        x = dr32/rs
        !g = -2.0*newton*M0/rs/rs
        !g = g*( log(1.0+x)/x/x - 1.0/(x*(1.0+x)) )

        g = -2.0*newton*M0/(((rs-rc)**2.0)*(dr32+rs))  / (rs*rs*x*x)
        g = g * (-rs*(rs-rc)*dr32 + (dr32+rs)*( rs*(rs-2.0*rc)*log(1.0+x)+rc*rc*log(1.0+dr32/rc) ) )

        grav(ii) = g * xCenter(ii)/dr32

     enddo


  else if (sweepDir .eq. SWEEP_Y) then          ! y-component

     tmpdr32 = xCenter(j)*xCenter(j) + zCenter(k)*zCenter(k) 

     do ii = 1, numCells
        
        dr32 = sqrt(yCenter(ii)*yCenter(ii) + tmpdr32)
        x = dr32/rs
        !g = -2.0*newton*M0/rs/rs
        !g = g*( log(1.0+x)/x/x - 1.0/(x*(1.0+x)) )

        g = -2.0*newton*M0/(((rs-rc)**2.0)*(dr32+rs))  / (rs*rs*x*x)
        g = g * (-rs*(rs-rc)*dr32 + (dr32+rs)*( rs*(rs-2.0*rc)*log(1.0+x)+rc*rc*log(1.0+dr32/rc) ) )

        grav(ii) = g * yCenter(ii)/dr32

     enddo

  else if (sweepDir .eq. SWEEP_Z) then          ! z-component

     tmpdr32 = xCenter(j)*xCenter(j) + yCenter(k)*yCenter(k) 

     do ii = 1, numCells
        
        dr32 = sqrt(zCenter(ii)*zCenter(ii) + tmpdr32)           
        x = dr32/rs
        !g = -2.0*newton*M0/rs/rs
        !g = g*( log(1.0+x)/x/x - 1.0/(x*(1.0+x)) )

        g = -2.0*newton*M0/(((rs-rc)**2.0)*(dr32+rs))  / (rs*rs*x*x)
        g = g * (-rs*(rs-rc)*dr32 + (dr32+rs)*( rs*(rs-2.0*rc)*log(1.0+x)+rc*rc*log(1.0+dr32/rc) ) )

        grav(ii) = g * zCenter(ii)/dr32

     enddo

  endif

!==============================================================================
#ifndef FIXEDBLOCKSIZE
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
#endif


  
  call Grid_releaseBlkPtr(blockID, solnVec)
  
  return
   
end subroutine Gravity_accelOneRow


