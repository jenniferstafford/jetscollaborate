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


REAL FUNCTION getDensity(distance)

  use Simulation_data

  
#include "Simulation.h"

  IMPLICIT NONE

  real :: densityBG, tempBG
  real, intent(IN) :: distance


  
 ! Background
       if (sim_densityProfile =="betacore") then
          ! beta model for the density profile
          densityBG = max(sim_rhoCore*(1.0 + (distance/sim_rCore)**2)**(-1.5*sim_densityBeta),&
                          sim_rhoFloor)

       else
          ! uniform background density and pressure
          densityBG = sim_rhoCore
       endif
       
       getDensity=densityBG
      
  
END FUNCTION getDensity



REAL FUNCTION getPressure(distance)

  use Simulation_data

  
#include "Simulation.h"

  IMPLICIT NONE

  real :: densityBG, tempBG
  real, intent(IN) :: distance


  
  ! Background
   densityBG = getDensity(distance)
       if (sim_densityProfile =="betacore") then
        
          
          ! isothermal two-temperature atmosphere
          tempBG = sim_Tout*(1.0+(distance/sim_rCoreT)**3)&
                       /(sim_Tout/sim_Tcore+(distance/sim_rCoreT)**3)
       else
          ! uniform background density and pressure
          tempBG = sim_Tcore
       endif
       
       getPressure=tempBG*densityBG*gasConst/sim_mu
      
  
END FUNCTION getPressure




!!****if* source/Simulation/SimulationMain/magnetoHD/MHD_Jet/Simulation_initBlock
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


!!REORDER(4): solnData

subroutine Simulation_initBlock(blockID)


  !use Tree
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
       Grid_getBlkPtr, Grid_releaseBlkPtr, &
       Grid_getCellCoords, Grid_putPointData

    
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get
  use Logfile_interface, ONLY : Logfile_stamp
  use Eos_interface, ONLY : Eos_wrapped
  use Simulation_data
  use Driver_data, ONLY : dr_simTime, dr_dtInit
  use Hydro_data, ONLY : hy_bref


  implicit none

#include "constants.h"
#include "Flash.h"
#include "Simulation.h"

  integer,intent(IN) :: blockID
  real,pointer,dimension(:,:,:,:) :: solnData, &
        solnFaceXData, solnFaceYData, solnFaceZData


  real :: rho_zone, velx_zone, vely_zone, velz_zone, temp_zone, pres_zone, &
       ener_zone, ekin_zone, eint_zone, gasConst, &
       bx_zone, by_zone, bz_zone, bx_face_1, by_face_1, bz_face_1, &
       bx_face_2,by_face_2,bz_face_2,magp_zone
  real :: densityBG, tempBG!, rhoCut
  real :: vel, fac


  integer :: i,j,k, istat
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  logical :: gcell = .true.

  integer :: nozzle=1
  real :: radius, length, sig, distance, theta
  real, dimension(3) :: cellvec, plnvec, jetvec, rvec, phivec, velvec, voutvec
  real :: r2, bf, rout, rmix

  
  real :: Simulation_a_cubic
  integer :: istat, sizeX, sizeY, sizeZ, i, j, k, ix, jy, kz, ii, jj, kk
  real,allocatable,dimension(:) :: xCoordf,yCoordf,zCoordf
  real, pointer, dimension(:,:,:,:) :: facexData,faceyData,facezData
  real :: delx,dely,delz


  integer,intent(IN) :: blockID ! removed myPE


  logical, save :: firstCall = .TRUE.

  real, allocatable, dimension(:,:,:) :: ax_block, ay_block, az_block



  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(xCoord(sizeX),stat=istat)
  allocate(xCoordf(sizeX+1),stat=istat)
  allocate(yCoord(sizeY),stat=istat)
  allocate(yCoordf(sizeY+1),stat=istat)
  allocate(zCoord(sizeZ),stat=istat)
  allocate(zCoordf(sizeZ+1),stat=istat)

  
  allocate(ax_block(sizeX,sizeY+1,sizeZ+1),stat=istat)
  allocate(ay_block(sizeX+1,sizeY,sizeZ+1),stat=istat)
  allocate(az_block(sizeX+1,sizeY+1,sizeZ),stat=istat)

  
  call Grid_getCellCoords(IAXIS,blockID,CENTER,gcell, xCoord, sizeX)
  call Grid_getCellCoords(JAXIS,blockID,CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(KAXIS,blockID,CENTER,gcell, zCoord, sizeZ)
  call Grid_getCellCoords(KAXIS, blockId, FACES, gcell, zCoordf, sizeZ+1)
  call Grid_getCellCoords(JAXIS, blockId, FACES, gcell, yCoordf, sizeY+1)
  call Grid_getCellCoords(IAXIS, blockId, FACES, gcell, xCoordf, sizeX+1)


  

  if (firstCall) then

     firstCall = .FALSE.

  endif
  
  call Grid_getBlkPtr(blockID,solnData,CENTER)
  call Grid_getBlkPtr(blockID,solnFaceXData,FACEX)
  call Grid_getBlkPtr(blockID,solnFaceYData,FACEY)
  call Grid_getBlkPtr(blockID,solnFaceZData,FACEZ)




  
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

           
           
	  
	   velx_zone=0.0
	   vely_zone=0.0
	   velz_zone=0.0
	   ekin_zone=0.0
           eint_zone = pres_zone/rho_zone*1.5
	   ener_zone=ekin_zone+eint_zone

    fac=0.0
          densityBG = getDensity(distance)
       
       
          solnData(DENS_VAR,i,j,k) = densityBG
          solnData(PRES_VAR,i,j,k) = getPressure(distance)
	   
           ! store the variables in the block's unk data
           
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


  call Eos_wrapped(MODE_DENS_PRES, blkLimitsGC, blockID, CENTER)
  solnData(ENER_VAR,:,:,:) = solnData(EINT_VAR,:,:,:)+&
                        0.5*(solnData(VELX_VAR,:,:,:)**2 +&
                             solnData(VELY_VAR,:,:,:)**2 +&
                             solnData(VELZ_VAR,:,:,:)**2)
  
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





  ! Initial conditions are spatially uniform.

  call PhysicalConstants_get("ideal gas constant", gasConst)
  
  rho_zone = sim_rhoCore
  temp_zone = sim_Tcore

  velx_zone = 0.0
  vely_zone = sim_windVel
  velz_zone = 0.0


  ! Compute the gas energy and set the gamma-values needed for
  ! the equation of state.
  !ekin_zone = 0.5 * (velx_zone**2 + vely_zone**2 + velz_zone**2)

  !eint_zone = pres_zone / (sim_gamma-1.) / rho_zone
  !ener_zone = eint_zone + ekin_zone
  !ener_zone = max(ener_zone, sim_smallP)


#if NSPECIES > 0
  solnData(SPECIES_BEGIN,:,:,:) =  1.0-(NSPECIES-1)*sim_smallx
  solnData(SPECIES_BEGIN+1:SPECIES_END,:,:,:) =     sim_smallx
#endif

  ! store the variables in the block's unk data
  solnData(DENS_VAR,:,:,:) = rho_zone
  solnData(TEMP_VAR,:,:,:) = temp_zone
  solnData(GAMC_VAR,:,:,:) = sim_gamma
  solnData(GAMC_VAR,:,:,:) = sim_gamma
  solnData(VELX_VAR,:,:,:) = velx_zone
  solnData(VELY_VAR,:,:,:) = vely_zone
  solnData(VELZ_VAR,:,:,:) = velz_zone

#ifdef EINT_VAR
  !solnData(EINT_VAR,:,:,:) = eint_zone
  !solnData(EINT_VAR,:,:,:) = solnData(PRES_VAR,:,:,:)/solnData(DENS_VAR,:,:,:)&
  !                                   /(solnData(GAMC_VAR,:,:,:)-1.0)
  solnData(EINT_VAR,:,:,:) = gasConst*solnData(TEMP_VAR,:,:,:)&
                             /(solnData(GAMC_VAR,:,:,:)-1.0)/sim_mu
#endif
  !solnData(ENER_VAR,:,:,:) = ener_zone
  solnData(ENER_VAR,:,:,:) = solnData(EINT_VAR,:,:,:)+&
                        0.5*(solnData(VELX_VAR,:,:,:)**2 +&
                             solnData(VELY_VAR,:,:,:)**2 +&
                             solnData(VELZ_VAR,:,:,:)**2)

  solndata(MAGX_VAR,:,:,:) = 0.0
  solndata(MAGY_VAR,:,:,:) = 0.0
  solndata(MAGZ_VAR,:,:,:) = sim_bzAmbient
  solndata(MAGP_VAR,:,:,:) = max((sim_bzAmbient/hy_bref)**2/2.0, sim_smallP)

  solnFaceXdata(MAG_FACE_VAR,:,:,:) = 0.0
  solnFaceYdata(MAG_FACE_VAR,:,:,:) = 0.0
  solnFaceZdata(MAG_FACE_VAR,:,:,:) = sim_bzAmbient

  !rhoCut = sim_rhoCore*(1.0 + (sim_rCut/sim_rCore)**2)**(-1.5*sim_densityBeta)

  bf = sim(nozzle)%rFeatherOut
  r2 = sim(nozzle)%radius
  rout = r2 + bf
  rmix = rout + sim(nozzle)%rFeatherMix

  ! Initialize the nozzle and environment
  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
   do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
    do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
       cellvec = (/ xCoord(i), yCoord(j), zCoord(k) /)
       call hy_uhd_jetNozzleGeometry(nozzle,cellvec,radius,length,distance,&
                                     sig,theta,jetvec,rvec,plnvec,phivec)
       ! inside the nozzle
     !  if ((radius.le.rmix)&
      !     .and.(abs(length).le.2.0*(sim(nozzle)%length+sim(nozzle)%zFeather))) then
       !   fac = taper(nozzle, radius, 0.5*length, 1.0, 1.0, 0.0)
          !vel = sim(nozzle)%velocity&
          !      *(0.5*(1.0+cos(PI*(max(0.0, min(1.0, (radius-r2)/bf)))))&
          !      *(1.0-sim(nozzle)%outflowR)+sim(nozzle)%outflowR ) &
          !      *sin(PI/2.0*min(abs(length),0.5*sim(nozzle)%length)*sig/sim(nozzle)%length/0.5)
          !voutvec = sim(nozzle)%outflowR*sim(nozzle)%velocity*plnvec&
          !          !*coshat(radius-0.5*(r2+2.0*bf), 0.5*(r2+bf), bf, 1.0)
          !          *0.5*(1.0+cos(PI*( max(-1.0, min(0.0,(radius-rout)/bf)) )))

          !velvec = vel*jetvec + voutvec &
          !         + sim(nozzle)%linVel + cross(sim(nozzle)%angVel,rvec*distance)
          !solnData(VELX_VAR:VELZ_VAR,i,j,k) = velvec*fac + solnData(VELX_VAR:VELZ_VAR,i,j,k)*(1.0-fac)
      ! endif
       ! cylindrical initial cavity

       
     !  if (sim(nozzle)%initGeometry == 'cylindrical') then
      !    if ((radius.le.rmix)&
       !       .and.(abs(length).le.2.0*(sim(nozzle)%length+sim(nozzle)%zFeather))) then
             ! inside the extended nozzle and feather
        !     fac = taper(nozzle, radius, 0.5*length, 1.0, 1.0, 0.0)
             !fac = 1.0
         ! else
          !   fac = 0.0
          !endif
       ! spherical initial cavity
     !  else
      !    if (distance.le.2.0*max( rmix, (sim(nozzle)%length+sim(nozzle)%zFeather) ) ) then
       !      fac = taperSph(nozzle, 0.5*distance, 1.0, 0.0)
        !  else
         !    fac = 0.0
          !endif
      ! endif

       ! Background

       fac=0.0
       densityBG = getDensity(distance)
       
       
       solnData(DENS_VAR,i,j,k) = densityBG
       solnData(PRES_VAR,i,j,k) = getPressure(distance)
       solnData(JET_SPEC,i,j,k) = max(sim_smallx, fac)
       solnData(ISM_SPEC,i,j,k) = max(sim_smallx, 1.0-fac)

    enddo
   enddo
  enddo

  call Eos_wrapped(MODE_DENS_PRES, blkLimitsGC, blockID, CENTER)
  solnData(ENER_VAR,:,:,:) = solnData(EINT_VAR,:,:,:)+&
                        0.5*(solnData(VELX_VAR,:,:,:)**2 +&
                             solnData(VELY_VAR,:,:,:)**2 +&
                             solnData(VELZ_VAR,:,:,:)**2)
  call Grid_releaseBlkPtr(blockID, solnData, CENTER)
  call Grid_releaseBlkPtr(blockID,solnFaceXData,FACEX)
  call Grid_releaseBlkPtr(blockID,solnFaceYData,FACEY)
  call Grid_releaseBlkPtr(blockID,solnFaceZData,FACEZ)
 

  return
end subroutine Simulation_initBlock



