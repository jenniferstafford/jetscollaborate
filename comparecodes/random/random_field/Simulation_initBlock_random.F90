!!
!! NAME
!!
!!  Simulation_a_cubic
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: direction
!!                       real(IN) :: x,y,z)
!!
!!
!!
!! DESCRIPTION
!!
!!  Calculates vector potential A at position
!!  x,y,z using cubic interpolation on root grid of
!!  A set up in Simulation_init, assuming period
!!  boundary conditions on A
!!
!! ARGUMENTS
!!
!!  diection -         component of A to be calculated
!!  	     	       1:x, 2:y, 3:z
!!  x,y,z -            location on grid
!!***


REAL FUNCTION Simulation_a_cubic(direction,xx,yy,zz)

  use Simulation_data

#include "Simulation.h"

  IMPLICIT NONE

  real, dimension(64) :: invec, xyzvec
  real, intent(IN) :: xx,yy,zz
  real :: x,y,z,deltax,deltay,deltaz,dx,dy,dz
  integer, intent(IN) :: direction
  integer :: i0,j0,k0,i,j,k,xindex,yindex,zindex

  ! periodic outside set boundaries  
  x = sim_xmin + MODULO(xx - sim_xmin,(sim_xmax-sim_xmin)*DBLE(NXAVEC)/DBLE(NXAVEC-1))
  y = sim_ymin + MODULO(yy - sim_ymin,(sim_ymax-sim_ymin)*DBLE(NYAVEC)/DBLE(NYAVEC-1))
  z = sim_zmin + MODULO(zz - sim_zmin,(sim_zmax-sim_zmin)*DBLE(NZAVEC)/DBLE(NZAVEC-1))
  
  ! find indices for point within periodic block
  i0 = MODULO(floor((x - sim_xmin)/(sim_xmax-sim_xmin)*(NXAVEC-1)) + 1,NXAVEC+1)
  j0 = MODULO(floor((y - sim_ymin)/(sim_ymax-sim_ymin)*(NYAVEC-1)) + 1,NYAVEC+1)
  k0 = MODULO(floor((z - sim_zmin)/(sim_zmax-sim_zmin)*(NZAVEC-1)) + 1,NZAVEC+1)
  
  ! step size of a-grid in each direction
  deltax = sim_ax_coords(2)-sim_ax_coords(1)
  deltay = sim_ay_coords(2)-sim_ay_coords(1)
  deltaz = sim_az_coords(2)-sim_az_coords(1)

  ! normalized offset for cubic interpolation
  dx = (x - sim_ax_coords(i0))/deltax - 0.5d0
  dy = (y - sim_ay_coords(j0))/deltay - 0.5d0
  dz = (z - sim_az_coords(k0))/deltaz - 0.5d0

  ! fill interpolation vectors
  do i=1,4
     do j=1,4
        do k=1,4
           xindex=MODULO(i0 + 2 - i,NXAVEC)+1
           yindex=MODULO(j0 + 2 - j,NYAVEC)+1
           zindex=MODULO(k0 + 2 - k,NZAVEC)+1
           SELECT CASE(direction)
              CASE(1)
                 invec(16*(i-1) + 4*(j-1) + k) = sim_Ax(xindex,yindex,zindex)
              CASE(2)
                 invec(16*(i-1) + 4*(j-1) + k) = sim_Ay(xindex,yindex,zindex)
              CASE(3)
                 invec(16*(i-1) + 4*(j-1) + k) = sim_Az(xindex,yindex,zindex)
              END SELECT
           xyzvec(16*(i-1) + 4*(j-1) + k) = dx**DBLE(4 - i)*dy**DBLE(4 - j)*dz**DBLE(4 - k)
        enddo
     enddo
  enddo

  ! interpolation step
  Simulation_a_cubic=dot_product(xyzvec,matmul(sim_cubic_matrix,invec))

END FUNCTION Simulation_a_cubic


!!****if* source/Simulation/SimulationMain/magnetoHD/Random_field_test/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID, 
!!                       integer(IN) :: myPE)
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.
!!
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!  myPE -             my processor number
!!***


subroutine Simulation_initBlock(blockID)

  !use Tree
  use Simulation_data
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getBlkPtr, Grid_releaseBlkPtr, &
    Grid_getCellCoords, Grid_putPointData

  use Logfile_interface, ONLY : Logfile_stamp


  implicit none

#include "constants.h"
#include "Flash.h"
#include "Simulation.h"

  real :: Simulation_a_cubic
  integer :: istat, sizeX, sizeY, sizeZ, i, j, k, ix, jy, kz, ii, jj, kk
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  real,allocatable,dimension(:) :: xCoordf,yCoordf,zCoordf
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData
  real :: delx,dely,delz
  logical :: gcell = .true.
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC

  integer,intent(IN) :: blockID ! removed myPE
  real,pointer :: solnData(:,:,:,:)
  real :: rho_zone, velx_zone, vely_zone, velz_zone, pres_zone, &
       ener_zone, ekin_zone, eint_zone, &
       bx_zone, by_zone, bz_zone, bx_face_1, by_face_1, bz_face_1, &
       bx_face_2,by_face_2,bz_face_2,magp_zone

  logical, save :: firstCall = .TRUE.

  real, allocatable, dimension(:,:,:) :: ax_block, ay_block, az_block
  
  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX),stat=istat)
  allocate(xCoordf(sizeX+1),stat=istat)
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY),stat=istat)
  allocate(yCoordf(sizeY+1),stat=istat)
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ),stat=istat)
  allocate(zCoordf(sizeZ+1),stat=istat)

  allocate(ax_block(sizeX,sizeY+1,sizeZ+1),stat=istat)
  allocate(ay_block(sizeX+1,sizeY,sizeZ+1),stat=istat)
  allocate(az_block(sizeX+1,sizeY+1,sizeZ),stat=istat)

  call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
  call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)
  call Grid_getCellCoords(KAXIS, blockId, FACES, gcell, zCoordf, sizeZ+1)
  call Grid_getCellCoords(JAXIS, blockId, FACES, gcell, yCoordf, sizeY+1)
  call Grid_getCellCoords(IAXIS, blockId, FACES, gcell, xCoordf, sizeX+1)

!===============================================================================

  if (firstCall) then

     firstCall = .FALSE.

  endif

  call Grid_getBlkPtr(blockID, solnData, CENTER)
  call Grid_getBlkPtr(blockID,facexData,FACEX)
  call Grid_getBlkPtr(blockID,faceyData,FACEY)
  call Grid_getBlkPtr(blockID,facezData,FACEZ)

  ! Cubic interpolation of A_x pontential on x-edge of block coordinates
  do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)+1
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)+1
        do i=  blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)
           ax_block(i,j,k) = Simulation_a_cubic(1,xCoord(i),yCoordf(j),zCoordf(k))
        enddo
     enddo
  enddo

  ! Cubic interpolation of A_y pontential on y-edge of block coordinates
  do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)+1
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        do i=  blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)+1
           ay_block(i,j,k) = Simulation_a_cubic(2,xCoordf(i),yCoord(j),zCoordf(k))
        enddo
     enddo
  enddo

  ! Cubic interpolation of A_z pontential on z-edge of block coordinates
  do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)+1
        do i=  blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)+1
           az_block(i,j,k) = Simulation_a_cubic(3,xCoordf(i),yCoordf(j),zCoord(k))
        enddo
     enddo
  enddo

  do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        do i=  blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)
           
           delx=xCoord(blkLimitsGC(LOW,IAXIS)+1)-&
                xCoord(blkLimitsGC(LOW,IAXIS))
           dely=yCoord(blkLimitsGC(LOW,JAXIS)+1)-&
                yCoord(blkLimitsGC(LOW,JAXIS))
           delz=zCoord(blkLimitsGC(LOW,KAXIS)+1)-&
                zCoord(blkLimitsGC(LOW,KAXIS))

           bx_zone=0.
           by_zone=0.
           bz_zone=0.
           magp_zone=0.
           bx_face_1=0.
           by_face_1=0.
           bz_face_1=0.
           bx_face_2=0.
           by_face_2=0.
           bz_face_2=0.

           ! B=curl(A)
           bx_face_1 = (az_block(i,j+1,k) - az_block(i,j,k))/dely - (ay_block(i,j,k+1) - ay_block(i,j,k))/delz
           bx_face_2 = (az_block(i+1,j+1,k) - az_block(i+1,j,k))/dely - (ay_block(i+1,j,k+1) - ay_block(i+1,j,k))/delz
           by_face_1 = (ax_block(i,j,k+1) - ax_block(i,j,k))/delz - (az_block(i+1,j,k) - az_block(i,j,k))/delx
           by_face_2 = (ax_block(i,j+1,k+1) - ax_block(i,j+1,k))/delz - (az_block(i+1,j+1,k) - az_block(i,j+1,k))/delx
           bz_face_1 = (ay_block(i+1,j,k) - ay_block(i,j,k))/delx - (ax_block(i,j+1,k) - ax_block(i,j,k))/dely
           bz_face_2 = (ay_block(i+1,j,k+1) - ay_block(i,j,k+1))/delx - (ax_block(i,j+1,k+1) - ax_block(i,j,k+1))/dely
           
           bx_zone=(0.5 * (bx_face_1 + bx_face_2))
           by_zone=(0.5 * (by_face_1 + by_face_2))
           bz_zone=(0.5 * (bz_face_1 + bz_face_2))

           magp_zone=0.5*(bx_zone**2 + by_zone**2 + bz_zone**2)
           
	   rho_zone=1
	   pres_zone=1
	   velx_zone=0.0
	   vely_zone=0.0
	   velz_zone=0.0
	   ekin_zone=0.0
           eint_zone = pres_zone/rho_zone*1.5
	   ener_zone=ekin_zone+eint_zone

	   if ((xCoord(i).gt.1.874).and.(yCoord(j).gt.1.874).and.(zCoord(k).gt.1.874)) then
	      rho_zone=10
	   endif
	   
           ! store the variables in the block's unk data
           solnData(DENS_VAR,i,j,k) = rho_zone
           solnData(PRES_VAR,i,j,k) = pres_zone
           solnData(ENER_VAR,i,j,k) = ener_zone
#ifdef EINT_VAR
           solnData(EINT_VAR,i,j,k) = eint_zone
#endif
           solnData(GAMC_VAR,i,j,k) = 5./3.
           solnData(GAME_VAR,i,j,k) = 5./3.
           
           solnData(VELX_VAR,i,j,k) = velx_zone
           solnData(VELY_VAR,i,j,k) = vely_zone
           solnData(VELZ_VAR,i,j,k) = velz_zone
           
           solnData(MAGX_VAR,i,j,k) = bx_zone
           solnData(MAGY_VAR,i,j,k) = by_zone
           solnData(MAGZ_VAR,i,j,k) = bz_zone
           solnData(MAGP_VAR,i,j,k) = magp_zone
           solnData(DIVB_VAR,i,j,k) = 0.

#if NSPECIES > 0
           solnData(SPECIES_BEGIN,i,j,k) =  1.0-1.0e-10
           solnData(SPECIES_BEGIN+1:SPECIES_END,i,j,k) = 1.0e-10
#endif

        enddo
     enddo
  enddo

  ! fill face-centered fields on extended block
  do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        do i=  blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)+1

           facexData(MAG_FACE_VAR,i,j,k) = (az_block(i,j+1,k) - az_block(i,j,k))/dely - (ay_block(i,j,k+1) - ay_block(i,j,k))/delz

        enddo
     enddo
  enddo

  do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)+1
        do i=  blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)

           faceyData(MAG_FACE_VAR,i,j,k)= (ax_block(i,j,k+1) - ax_block(i,j,k))/delz - (az_block(i+1,j,k) - az_block(i,j,k))/delx

        enddo
     enddo
  enddo
  do k = blkLimitsGC(LOW, KAXIS), blkLimitsGC(HIGH, KAXIS)+1
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        do i=  blkLimitsGC(LOW, IAXIS), blkLimitsGC(HIGH, IAXIS)

           facezData(MAG_FACE_VAR,i,j,k)=(ay_block(i+1,j,k) - ay_block(i,j,k))/delx - (ax_block(i,j+1,k) - ax_block(i,j,k))/dely

        enddo
     enddo
  enddo

  call Grid_releaseBlkPtr(blockID,solnData, CENTER)
  call Grid_releaseBlkPtr(blockID,facexData,FACEX)
  call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
  call Grid_releaseBlkPtr(blockID,facezData,FACEZ)

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)
  deallocate(xCoordf)
  deallocate(yCoordf)
  deallocate(zCoordf)
  deallocate(ax_block)
  deallocate(ay_block)
  deallocate(az_block)

  return

end subroutine Simulation_initBlock
