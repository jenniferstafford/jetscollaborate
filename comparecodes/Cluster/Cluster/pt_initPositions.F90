!!****if* source/Simulation/SimulationMain/sb/pt_initPositions
!!
!! NAME
!!
!!  pt_initPositions
!!
!! SYNOPSIS
!!
!!  pt_initPositions(integer, INTENT(in)  :: blockid)
!!
!! DESCRIPTION
!!    An override for the SB simulation. 
!!
!!
!! ARGUMENTS
!!
!!   blockid : ID of block in current processor
!!
!!
!!***

subroutine pt_initPositions (blockID, success)

  use Simulation_data

  use Grid_interface, ONLY : Grid_getBlkPhysicalSize, &
    Grid_getBlkCenterCoords
  use Particles_data, ONLY: pt_numLocal, particles, pt_maxPerProc
  use Grid_data, ONLY : gr_myPE

#include "Flash.h"
#include "constants.h"

  implicit none

  integer, INTENT(in) :: blockID
  logical, intent(out) :: success
  integer       :: i, p
  logical       :: IsInBlock
  real          :: xpos, ypos, zpos, bxLower, byLower, bzLower, bxUpper, byUpper, bzUpper
  real          :: xvel, yvel, zvel, pm, blockSize(3), blockCenter(3)

  real :: pscale, vscale, mscale

!-------------------------------------------------------------------------------

! Particle slot number (incremented and saved between calls)

  p = pt_numLocal

!-------------------------------------------------------------------------------

    !
    ! unit conversion

    pscale = 1.0
    vscale = 1.0
    mscale = 1.0


! Get locations of block faces.

  call Grid_getBlkPhysicalSize(blockID, blockSize)
  call Grid_getBlkCenterCoords(blockID, blockCenter)
  bxLower    = blockCenter(1) - 0.5*blockSize(1)
  bxUpper    = blockCenter(1) + 0.5*blockSize(1)
  if (NDIM >= 2) then
     byLower = blockCenter(2) - 0.5*blockSize(2)
     byUpper = blockCenter(2) + 0.5*blockSize(2)
  endif
  if (NDIM == 3) then
     bzLower = blockCenter(3) - 0.5*blockSize(3)
     bzUpper = blockCenter(3) + 0.5*blockSize(3)
  endif

! Loop over both particles and compute their positions.

  do i = 1, ntot       !

       xpos = pscale*Positions_Ascii(1,i)
       ypos = pscale*Positions_Ascii(2,i)
       zpos = pscale*Positions_Ascii(3,i)
       !
       xvel = vscale*Velocities_Ascii(1,i)
       yvel = vscale*Velocities_Ascii(2,i)
       zvel = vscale*Velocities_Ascii(3,i)
       !
       pm = mscale*Masses_Ascii(i)

! Check to see if the particle lies within this block.

     IsInBlock = (xpos >= bxLower) .and. (xpos < bxUpper)
     if (NDIM >= 2) &
          IsInBlock = IsInBlock .and. ((ypos >= byLower) .and. (ypos < byUpper))
     if (NDIM == 3) &
          IsInBlock = IsInBlock .and. ((zpos >= bzLower) .and. (zpos < bzUpper))

! If yes, and adequate particle buffer space is available, initialize it.

     if (IsInBlock) then
        p = p + 1
        if (p > pt_maxPerProc) then
           call Driver_abortFlash &
                ("InitParticlePositions:  Exceeded max # of particles/processor!")
        endif

! Particle current block number.
        particles(BLK_PART_PROP,p) = real(blockID)
! Particle mass.

        particles(MASS_PART_PROP,p) = pm

! Particle position and velocity.
        particles(POSX_PART_PROP,p) = xpos
        particles(VELX_PART_PROP,p) = xvel
        particles(POSY_PART_PROP,p) = ypos
        particles(VELY_PART_PROP,p) = yvel
        particles(POSZ_PART_PROP,p) = zpos
        particles(VELZ_PART_PROP,p) = zvel

     endif
  enddo
    
! Set the particle database local number of particles.

  pt_numLocal = p

  success = .true.
  return

!-------------------------------------------------------------------------------
  
end subroutine pt_initPositions
