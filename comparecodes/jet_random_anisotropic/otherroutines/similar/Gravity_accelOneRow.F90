!!****if* source/physics/Gravity/GravityMain/PointMass/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow 
!!
!! SYNOPSIS
!!
!!  call Gravity_accelOneRow(integer(IN)  :: pos(2),
!!                           integer(IN)  :: sweepDir,
!!                           integer(IN)  :: blockID,
!!                           integer(IN)  :: numCells,
!!                           real(INOUT)  :: grav(numCells),
!!                           integer(IN),optional :: potentialIndex,
!!                           integer(IN),optional :: extraAccelVars(MDIM))
!!
!! DESCRIPTION
!!
!!  This routine computes the gravitational acceleration for a row
!!  of cells in a specified direction in a given block.
!!
!! ARGUMENTS
!!
!!  pos      :  Row indices transverse to the sweep direction
!!  sweepDir :    The sweep direction:  allowed values are 
!!              SWEEP_X, SWEEP_Y, and SWEEP_Z. These values are defined
!!              in constants.h.
!!  blockID  :  The local identifier of the block to work on
!!  numCells :  Number of cells to update in grav()
!!  grav()   :   Array to receive result
!!  potentialIndex :  optional, not applicable in pointmass gravity
!!  extraAccelVars :  optional, ignored in this implementation
!! 
!!***

subroutine Gravity_accelOneRow (pos, sweepDir, blockID, numCells, grav, &
                                potentialIndex, extraAccelVars)

!=======================================================================

  use Gravity_data, ONLY: useGravity
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use Simulation_data

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: sweepDir,blockID,numCells
  integer, dimension(2),INTENT(in) ::pos
  real, dimension(numCells),INTENT(inout) :: grav
  integer,intent(IN),optional :: potentialIndex
  integer,intent(IN),OPTIONAL :: extraAccelVars(MDIM)

!==========================================================================


#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_KHI_GC) :: zCenter
  real,dimension(GRID_JHI_GC) :: yCenter
  real,dimension(GRID_IHI_GC) :: xCenter
#else
  real,allocatable,dimension(:) ::xCenter,yCenter,zCenter
  integer, dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
#endif
  real :: r, rt, gasConst

  integer :: sizeX,sizeY,sizez

  integer :: ii,j,k,nozzle=1
  logical :: gcell = .true.

!==============================================================================

  if (.NOT.useGravity .or. sim_densityProfile=="uniform") return
  call PhysicalConstants_get("ideal gas constant", gasConst)
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
  call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, xCenter, sizeX)
  xCenter = xCenter - sim(nozzle)%pos(1)
  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, yCenter, sizeY)
     yCenter = yCenter - sim(nozzle)%pos(2)
  endif
  if (NDIM == 3) then 
     call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zCenter, sizeZ)
     zCenter = zCenter - sim(nozzle)%pos(3)
  endif
  

  if (sweepDir .eq. SWEEP_X) then                       ! x-component

     rt = yCenter(j)*yCenter(j) + zCenter(k)*zCenter(k) 

     do ii = 1, numCells
        r=sqrt(xCenter(ii)*xCenter(ii) + rt)
        if (r .lt. sim_rCut) then
           grav(ii) = gasConst/sim_mu*xCenter(ii)&
                      ! Density derivative
                      *(-3.*sim_densityBeta/(1.+r*r/sim_rCore**2)/sim_rCore**2&
                      ! Temperature profile
                      *sim_Tout*(1.0+(r/sim_rCoreT)**3)&
                      /(sim_Tout/sim_Tcore+(r/sim_rCoreT)**3)&
                      ! Temperature derivative
                      +3.0*r*sim_Tout*(sim_Tout/sim_Tcore-1.0)*(sim_rCoreT)**3&
                      /(sim_Tout/sim_Tcore*(sim_rCoreT)**3+r**3)**2)
        else
           grav(ii) = 0.
        endif
           !grav(ii) = -2.*1.5*sim_densityBeta*sim_pAmbient/sim_rhoCore/&
           !     (1.0 + r2/sim_rCore**2)*xCenter(ii)/sim_rCore**2
     enddo
     
  else if (sweepDir .eq. SWEEP_Y) then          ! y-component

     rt = xCenter(j)*xCenter(j) + zCenter(k)*zCenter(k) 

     do ii = 1, numCells
        r=sqrt(yCenter(ii)*yCenter(ii) + rt)
        if (r .lt. sim_rCut) then
           grav(ii) = gasConst/sim_mu*yCenter(ii)&
                      ! Density derivative
                      *(-3.*sim_densityBeta/(1.+r*r/sim_rCore**2)/sim_rCore**2&
                      ! Temperature profile
                      *sim_Tout*(1.0+(r/sim_rCoreT)**3)&
                      /(sim_Tout/sim_Tcore+(r/sim_rCoreT)**3)&
                      ! Temperature derivative
                      +3.0*r*sim_Tout*(sim_Tout/sim_Tcore-1.0)*(sim_rCoreT)**3&
                      /(sim_Tout/sim_Tcore*(sim_rCoreT)**3+r**3)**2)
        else
           grav(ii) = 0.
        endif
           !grav(ii) = -2.*1.5*sim_densityBeta*sim_pAmbient/sim_rhoCore/&
           !     (1.0 + r2/sim_rCore**2)*yCenter(ii)/sim_rCore**2
     enddo
     
  else if (sweepDir .eq. SWEEP_Z) then          ! z-component

     rt = xCenter(j)*xCenter(j) + yCenter(k)*yCenter(k) 

     do ii = 1, numCells
        r=sqrt(zCenter(ii)*zCenter(ii) + rt)
        if (r .lt. sim_rCut) then
           grav(ii) = gasConst/sim_mu*zCenter(ii)&
                      ! Density derivative
                      *(-3.*sim_densityBeta/(1.+r*r/sim_rCore**2)/sim_rCore**2&
                      ! Temperature profile
                      *sim_Tout*(1.0+(r/sim_rCoreT)**3)&
                      /(sim_Tout/sim_Tcore+(r/sim_rCoreT)**3)&
                      ! Temperature derivative
                      +3.0*r*sim_Tout*(sim_Tout/sim_Tcore-1.0)*(sim_rCoreT)**3&
                      /(sim_Tout/sim_Tcore*(sim_rCoreT)**3+r**3)**2)
        else
           grav(ii)=0.
        endif
           !grav(ii) = -2.*1.5*sim_densityBeta*sim_pAmbient/sim_rhoCore/&
           !     (1.0 + r2/sim_rCore**2)*zCenter(ii)/sim_rCore**2
     enddo

  endif

!==============================================================================
#ifndef FIXEDBLOCKSIZE
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
#endif

  return

end subroutine Gravity_accelOneRow
