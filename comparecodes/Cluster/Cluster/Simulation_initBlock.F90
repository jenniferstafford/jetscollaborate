!!****if* source/Simulation/SimulationMain/sb/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!! 
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer :: blockId, 
!!                       integer :: myPE)
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid data (density, pressure, velocity, etc.) for
!!  a specified block.  
!!
!!
!! ARGUMENTS
!!
!!  blockId -        The number of the block to initialize
!!  myPE -          current processor number
!!
!!
!!***

subroutine Simulation_initBlock (blockId, myPE)

  use Simulation_data
  use Eos_data
  use Eos_interface
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
      Grid_getCellCoords, Grid_putPointData  
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_interface, ONLY : Driver_abortFlash

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"

  integer,intent(IN) ::  blockId
  integer,intent(IN) ::  myPE
  
  integer  ::  i, j, k, n, jLo, jHi
  integer  ::  ii, jj, kk
  real     ::  sumRho, sumP, sumVX, sumVY, sumVZ
  real     ::  vel, diagonal
  real     ::  xx, dxx, yy, dyy, zz, dzz, frac
  real     ::  xf, yf, zf, rxf, ryf, rzf, r0
  real     ::  vx, vy, vz, p, e, ek
  logical  ::  validGeom
  integer  ::  istat
  integer  ::  nxb, nyb, nzb, lrefine_max

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  real,allocatable,dimension(:) :: xCoordf,yCoordf,zCoordf
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  integer,dimension(MDIM) :: axis

  real :: tempZone, eintZone, ekinZone
  real, pointer, dimension(:,:,:,:) :: solnData, facexData, faceyData, facezData

  logical :: gcell = .true.

  logical :: first

  real :: integral, dx, g, d_integral, x, y, r, Bfield0, Btmp, xm, r_ii
!  real, allocatable, dimension(:) :: T, rho
!  integer :: nmax


! get the coordinate information for the current block from the database

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)
  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
  allocate(xCoord(sizeX),stat=istat)
  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
  allocate(yCoord(sizeY),stat=istat)
  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
  allocate(zCoord(sizeZ),stat=istat)


  allocate(xCoordf(sizeX+1),stat=istat)  ! plus 1 ???
  allocate(yCoordf(sizeY+1),stat=istat)
  allocate(zCoordf(sizeZ+1),stat=istat)


  nxb = blkLimits(HIGH,IAXIS) - blkLimits(LOW,IAXIS) + 1
  nyb = blkLimits(HIGH,JAXIS) - blkLimits(LOW,JAXIS) + 1
  nzb = blkLimits(HIGH,KAXIS) - blkLimits(LOW,KAXIS) + 1

  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, CENTER,gcell, yCoord, sizeY)
  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)


  if (NDIM == 3) call Grid_getCellCoords&
                      (KAXIS, blockId, FACES, gcell, zCoordf, sizeZ+1)
  if (NDIM >= 2) call Grid_getCellCoords&
                      (JAXIS, blockId, FACES, gcell, yCoordf, sizeY+1)
  call Grid_getCellCoords(IAXIS, blockId, FACES, gcell, xCoordf, sizeX+1)


  call Grid_getBlkPtr(blockId,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_getBlkPtr(blockID,facexData,FACEX)
     call Grid_getBlkPtr(blockID,faceyData,FACEY)
     if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

! setup 1D HSE solution from r=0 to r=xmax. Use this below to look up (i.e., essentially interpolate)
! the solution within r=0 and r=sqrt(3*(2*xmax)^2) = 2*sqrt(3) * xmax = 3.461 * xmax
! note that 0>xmin<x<xmax>0 and that we need to define rho & T outside the box as well
! as we require isolated boundary conditions

!  call RuntimeParameters_get('lrefine_max', lrefine_max)
!  nmax = 1000 * nxb * (2**(lrefine_max-1))
!
!  allocate(T(nmax))
!  allocate(rho(nmax))

  xm = abs(sim_xMin)

!  dx = 4.0*sim_xMax/dble(nmax)  !4.0>3.461
!  integral = 0.0
!
!  do i = 1, nmax
!    r = dx*dble(i)
!    x = r/rs
!    g = -2.0*newton*M0/rs/rs
!    g = g*( log(1.0+x)/x/x - 1.0/(x*(1.0+x)) )
!    y = r/sim_xMax
!    if (y .lt. 1.0) then
!       T(i) = T_in + 3.0*DeltaT*y*y -2.0*DeltaT*y*y*y
!    else
!       T(i) = T_out
!    endif
!    d_integral = sim_eos_singleSpeciesA*mp*g/(kb*T(i))*dx
!    integral = integral + d_integral
!    rho(i) = rho_in * T_in/T(i)	* exp(integral)
!  enddo


  dx = dr

  first = .true.
  do i = 1, ngas
    r = dx*dble(i)
    ! normalizing magnetic field
    if (r > 0.2*0.1*sim_xMax .and. first) then
        Bfield0 = sqrt(2.0)*sqrt( sim_Bfield0 * eos_gasConstant*tempgas(i)*rhogas(i)/sim_eos_singleSpeciesA )
	!Bfield0 = sqrt(2.0)*sqrt( sim_Bfield0 * eos_gasConstant*T(i)*rho(i)/sim_eos_singleSpeciesA )
        !!Bfield0 = Bfield0 * sqrt(4.0e0*PI)
	!r0 = 0.1*sim_xMax
	first = .false.
    endif
  enddo 


  do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
     zz = zCoord(k)
     zf = zCoordf(k)
     do j = blkLimitsGC(LOW, JAXIS), blkLimitsGC(HIGH, JAXIS)
        yy = yCoord(j)
	yf = yCoordf(j)
        do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH, IAXIS)
           xx = xCoord(i)
	   xf = xCoordf(i)

              ! interpolating the ICs to get T[] and rho[] at equal intervals of dx = dr
	      ! this assumes that the cluster center is at (0,0,0)
	      r  = sqrt(xx**2.0 + yy**2.0 + zz**2.0)
	      !ii = int((r-0.5*dx)/dx) + 1
              ii = int(r/dx) + 1
              !r_ii = ii*dr - 0.5*dr
              r_ii = (ii-1)*dr 

              solnData(DENS_VAR,i,j,k) = rhogas(ii) + (rhogas(ii+1)-rhogas(ii))*( (r-r_ii)/dr )
	      tempZone = tempgas(ii) + (tempgas(ii+1)-tempgas(ii))*( (r-r_ii)/dr )

!              solnData(PRES_VAR,i,j,k) = eos_gasConstant*tempZone*rho(ii)/sim_eos_singleSpeciesA
              solnData(PRES_VAR,i,j,k) = eos_gasConstant*tempZone*solnData(DENS_VAR,i,j,k)/sim_eos_singleSpeciesA

              solnData(VELX_VAR,i,j,k) = 0.0
              solnData(VELY_VAR,i,j,k) = 0.0
              solnData(VELZ_VAR,i,j,k) = 0.0

              ekinZone = 0.5 * dot_product(solnData(VELX_VAR:VELZ_VAR,i,j,k),&
                                           solnData(VELX_VAR:VELZ_VAR,i,j,k))

              eintZone = solnData(PRES_VAR,i,j,k)/(sim_gamma-1.0)/solnData(DENS_VAR,i,j,k)

              solnData(ENER_VAR,i,j,k) = eintZone + ekinZone
              solnData(EINT_VAR,i,j,k) = eintZone

              solnData(GAMC_VAR,i,j,k) = sim_gamma
              solnData(GAME_VAR,i,j,k) = sim_gamma

#if NFACE_VARS > 0
               solnData(DIVB_VAR,i,j,k) = 0.0
!               solnData(MAGX_VAR,i,j,k) = Bfield0*xx/r /((r/r0)**2.0) !sim_Bfield0
!               solnData(MAGY_VAR,i,j,k) = Bfield0*yy/r /((r/r0)**2.0)
!               solnDatA(MAGZ_VAR,i,j,k) = Bfield0*zz/r /((r/r0)**2.0)

               ! bx, by, bz are the arrays of the size = ( nxb*2^(lrefmin-1) + 2 )^3 from IDL

               call TDint(bx,ntotal,xm,xx,yy,zz,Btmp)
               solnData(MAGX_VAR,i,j,k) = Btmp*Bfield0
               call TDint(by,ntotal,xm,xx,yy,zz,Btmp)
               solnData(MAGY_VAR,i,j,k) = Btmp*Bfield0
               call TDint(bz,ntotal,xm,xx,yy,zz,Btmp)
               solnDatA(MAGZ_VAR,i,j,k) = Btmp*Bfield0

               solnData(MAGP_VAR,i,j,k) = 0.5*dot_product(solnData(MAGX_VAR:MAGZ_VAR,i,j,k),&
                                                          solnData(MAGX_VAR:MAGZ_VAR,i,j,k))

#endif

              ! Cell face-centered variables for StaggeredMesh scheme
#if NFACE_VARS > 0
              if (sim_killdivb) then

!                 rxf = sqrt(xf**2.0 + yy**2.0 + zz**2.0)
!                 ryf = sqrt(xx**2.0 + yf**2.0 + zz**2.0)
!                 rzf = sqrt(xx**2.0 + yy**2.0 + zf**2.0)
!                 facexData(MAG_FACE_VAR,i,j,k)= Bfield0*xf/rxf /((rxf/r0)**2.0)
!                 faceyData(MAG_FACE_VAR,i,j,k)= Bfield0*yf/ryf /((ryf/r0)**2.0)
!                 facezData(MAG_FACE_VAR,i,j,k)= Bfield0*zf/rzf /((rzf/r0)**2.0)

                 call TDint(bx,ntotal,xm,xf,yy,zz,Btmp)
                 facexData(MAG_FACE_VAR,i,j,k)= Btmp*Bfield0
                 call TDint(by,ntotal,xm,xx,yf,zz,Btmp)
                 faceyData(MAG_FACE_VAR,i,j,k)= Btmp*Bfield0
                 call TDint(bz,ntotal,xm,xx,yy,zf,Btmp)
                 facezData(MAG_FACE_VAR,i,j,k)= Btmp*Bfield0

              endif
#endif

           enddo
        enddo
     enddo

#if NFACE_VARS > 0
      if (sim_killdivb) then

         do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)

            i=blkLimitsGC(HIGH,IAXIS)+1
            do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
               zz = zCoord(k)
               yy = yCoord(j)
               xf = xCoordf(i)
!               rxf = sqrt(xf**2.0 + yy**2.0 + zz**2.0)
!               facexData(MAG_FACE_VAR,i,j,k) = Bfield0*xf/rxf /((rxf/r0)**2.0)

               call TDint(bx,ntotal,xm,xf,yy,zz,Btmp)
               facexData(MAG_FACE_VAR,i,j,k) = Btmp*Bfield0

            enddo

            j=blkLimitsGC(HIGH,JAXIS)+1
            do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
               zz = zCoord(k)
               yf = yCoordf(j)
               xx = xCoord(i)
!               ryf = sqrt(xx**2.0 + yf**2.0 + zz**2.0)
!               faceyData(MAG_FACE_VAR,i,j,k) = Bfield0*yf/ryf /((ryf/r0)**2.0)

               call TDint(by,ntotal,xm,xx,yf,zz,Btmp)
               faceyData(MAG_FACE_VAR,i,j,k) = Btmp*Bfield0

            enddo

         enddo

         k=blkLimitsGC(HIGH,KAXIS)+1
         do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
            do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
               zf = zCoordf(k)
               yy = yCoord(j)
               xx = xCoord(i)
!               rzf = sqrt(xx**2.0 + yy**2.0 + zf**2.0)
!               facezData(MAG_FACE_VAR,i,j,k) = Bfield0*zf/rzf /((rzf/r0)**2.0)


!               write(*,*)'before',i,j,k,xx,yy,zz, xm

               call TDint(bz,ntotal,xm,xx,yy,zf,Btmp)

!               write(*,*)'after',i,j,k,xx,yy,zz, xm




               facezData(MAG_FACE_VAR,i,j,k) = Btmp*Bfield0

            enddo
         enddo

     endif

#endif

  ! Release pointer
  call Grid_releaseBlkPtr(blockId,solnData,CENTER)

#if NFACE_VARS > 0
  if (sim_killdivb) then
     call Grid_releaseBlkPtr(blockID,facexData,FACEX)
     call Grid_releaseBlkPtr(blockID,faceyData,FACEY)
     call Grid_releaseBlkPtr(blockID,facezData,FACEZ)
  endif
#endif

  deallocate(xCoord)
  deallocate(yCoord)
  deallocate(zCoord)

!  deallocate(T)
!  deallocate(rho)

  return
end subroutine Simulation_initBlock


subroutine TDint(b_arr,N,xm,xx,yy,zz,b)

implicit none

real :: b000, b100, b010, b001, b101, b011, b110, b111, b, x,y,z, xx,yy,zz, dx, xm
integer :: i,j,k, n
real :: b_arr(N,N,N)

! if any of the coordinates are outside the box, set B-field to zero. Note that the values
! in those boundary regions will be set by the BC set in flash.par (e.g., "reflect")

if (abs(xx).gt.xm .or. abs(yy).gt.xm .or. abs(zz).gt.xm) then

   b = 0.0

else

  dx = 2.0*xm/(N-2)

  i = int((xx+xm+0.5*dx)/dx) + 1
  j = int((yy+xm+0.5*dx)/dx) + 1
  k = int((zz+xm+0.5*dx)/dx) + 1

  if (i .lt.0 .or. j.lt.0 .or. k.lt.0) then

     write(*,*)'kupa',i,j,k,xx,yy,zz
     stop

  endif

  x = (xx+xm+0.5*dx)/dx -i +1
  y = (yy+xm+0.5*dx)/dx -j +1
  z = (zz+xm+0.5*dx)/dx -k +1

  b000 = b_arr(i,j,k)
  b100 = b_arr(i+1,j,k)
  b010 = b_arr(i,j+1,k)
  b001 = b_arr(i,j,k+1)
  b101 = b_arr(i+1,j,k+1)
  b011 = b_arr(i,j+1,k+1)
  b110 = b_arr(i+1,j+1,k)
!  b000 = b_arr(i,j,k) 
  b111 = b_arr(i+1,j+1,k+1)


  b = b000*(1.0-x)*(1.0-y)*(1.0-z)  + &
      b100*     x *(1.0-y)*(1.0-z)  + &
      b010*(1.0-x)*     y *(1.0-z)  + &
      b001*(1.0-x)*(1.0-y)*     z   + &
      b101*     x *(1.0-y)*     z   + &
      b011*(1.0-x)*     y *     z   + &
      b110*     x *     y *(1.0-z)  + &
      b111*     x *     y *     z 

endif

end subroutine TDint
