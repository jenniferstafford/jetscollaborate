!!****if* source/physics/Hydro/HydroMain/unsplit/hy_uhd_addThermalFluxes
!!
!! NAME
!!
!!  hy_uhd_addThermalFluxes
!!
!!
!! SYNOPSIS
!!
!!  hy_uhd_addThermalFluxes(blockID,blkLimitsGC,ix,iy,iz,Flux,visc,sweepDir,Dcff)
!!
!!  hy_uhd_addThermalFluxes(integer(IN) :: blockID,
!!                          integer(IN) :: blkLimitsGC(:,:),
!!                          integer(IN) :: ix,
!!                          integer(IN) :: iy,
!!                          integer(IN) :: iz,
!!                          real(IN)    :: Flux,
!!                          real(IN)    :: cond(:,:,:),
!!                          integer(IN) :: sweepDir)
!!
!!
!! DESCRIPTION
!!
!!  Adds thermal flux contributions to total fluxes
!!
!!
!! ARGUMENTS
!!
!!  blockID     - a local blockID
!!  blkLimitsGC - an array that holds the lower and upper indices of the section
!!                of block with the guard cells
!!  ix,iy,iz    - indices of the line along which the sweep is made
!!  Flux        - array containing MHD fluxes
!!  cond        - array containing thermal conductivity coefficients
!!  sweepDir    - direction of sweep
!!
!!***

Subroutine hy_uhd_addThermalFluxes(blockID,blkLimitsGC,ix,iy,iz,Flux,cond,sweepDir,Dcff)

  use Hydro_data,     ONLY : hy_qref
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr, Grid_getDeltas,&
	                     Grid_getBlkIndexLimits, Grid_getCellCoords

  ! SB modification -------------
  use Eos_data, ONLY : eos_gasConstant, eos_gamma, eos_singleSpeciesA
  use Conductivity_data, ONLY : SpitzerFraction, isotropic
!  use MHD_data,       ONLY : mhd_qref
  use Cosmology_interface, ONLY : Cosmology_getRedshift
  use Simulation_data, ONLY : sim_saturatedConduction
  ! --------------------

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"

  !! Argument List ----------------------------------------------------------
  integer, INTENT(IN) :: blockID,ix,iy,iz
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
  real, dimension(HY_VARINUM), intent(INOUT) :: Flux
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, &
                  GRID_JLO_GC:GRID_JHI_GC, &
                  GRID_KLO_GC:GRID_KHI_GC),intent(IN) :: cond
#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
                  intent(IN) :: cond
#endif
  integer, INTENT(IN) :: sweepDir
  real, intent(OUT) :: Dcff
  !! ----------------------------------------------------------------------

  real    :: kappa_th,idx,idy,idz
  real, dimension(MDIM) :: del
  real, pointer, dimension(:,:,:,:) :: U


  ! SB modification
  real, pointer, dimension(:,:,:,:) :: facexData, faceyData, facezData
  real :: modB, bx, by, bz, Fthx, Fthy, Fthz, FthB, lambda_mfp
  real :: gradx, grady, gradz, fx, fy, fz, signx, signy, signz
  real :: qsata, qsatb, qsat, cond_eff, gTx, gTy, gTz, mod_gT, gr, F_sat, ff

  real :: modB_im1jk, modB_im1jp1k, modB_ijp1k, modB_ijk, modB_im1jkp1,&
	  modB_ijkp1, modB_ijm1k, modB_ip1jm1k, modB_ip1jk, modB_ijm1kp1, &
	  modB_ijkm1, modB_ip1jkm1, modB_ijp1km1
  real :: bxx, bxy, bxz, al, bl, ar, br, byy, byx, byz, bzz, bzx, bzy
  real :: dtdyxl, dtdyxr, dtdyx, dtdzxl, dtdzxr, dtdzx, dtdxyl, dtdxyr, dtdxy, &
	  dtdzyl, dtdzyr, dtdzy, dtdxzl, dtdxzr, dtdxz, dtdyzl, dtdyzr, dtdyz
  real :: FthBx, FthBx_sat, FthBy, FthBy_sat, FthBz, FthBz_sat
  real :: a_tmp1, a_tmp2, b_tmp1, b_tmp2
  real :: qsatyxl, qsatyxr, qsatyx, qsatzxl, qsatzxr, qsatzx, qsatxyl, qsatxyr, &
	  qsatxy, qsatzyl, qsatzyr, qsatzy, qsatxzl, qsatxzr, qsatxz, qsatyzl, qsatyzr, qsatyz

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer :: sizeX,sizeY,sizeZ, istat
  integer,dimension(MDIM) :: axis
  logical :: gcell = .true.

  real :: Bx_mean, By_mean, Bz_mean
  real :: cv, cs, xt, ft, qq, cent, dist
!  logical :: saturated
  real :: currentRedshift, temp_max
  integer :: i,j,k

  ! ---------------------

   call Cosmology_getRedshift(currentRedshift)
!   sim_saturatedConduction = .false.  !.true.   !.false.


  !! Get deltas
  call Grid_getDeltas(blockID,del)

  idx=1./del(DIR_X)
  if (NDIM >= 2) then
     idy=1./del(DIR_Y)
     if (NDIM == 3) then
        idz=1./del(DIR_Z)
     endif
  endif

  !! Get pointer
  call Grid_getBlkPtr(blockID,U,CENTER)


!  temp_max = 1.0e8
!  do i = ix-1, ix+1 
!    do j = iy-1, iy+1
!      do k = iz-1, iz
!        if ( U(TEMP_VAR,i,j,k) .gt. temp_max) U(TEMP_VAR,i,j,k) = temp_max
!      enddo
!    enddo
!  enddo


  call Grid_getBlkPtr(blockID,facexData,FACEX)
  call Grid_getBlkPtr(blockID,faceyData,FACEY)
  if (NDIM == 3) call Grid_getBlkPtr(blockID,facezData,FACEZ)


  cv = eos_gasConstant/(eos_gamma-1.0)/eos_singleSpeciesA
  cs = sqrt(eos_gamma*(eos_gamma-1.0)*cv*U(TEMP_VAR,ix,iy,iz))





!  Dcff  = cond(ix,iy,iz)/(U(DENS_VAR,ix,iy,iz)*cv)
!
!  sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
!  allocate(xCoord(sizeX),stat=istat)
!  sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
!  allocate(yCoord(sizeY),stat=istat)
!  sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
!  allocate(zCoord(sizeZ),stat=istat)
!  call Grid_getCellCoords(IAXIS, blockId, CENTER, gcell, xCoord, sizeX)
!  call Grid_getCellCoords(JAXIS, blockId, CENTER, gcell, yCoord, sizeY)
!  call Grid_getCellCoords(KAXIS, blockId, CENTER, gcell, zCoord, sizeZ)
!
!  cent = 1.974833651E+26 / 2.0
!  dist = sqrt( (zCoord(iz)-cent)**2.0 + (yCoord(iy)-cent)**2.0 + (xCoord(ix)-cent)**2.0 )



  select case(sweepDir)

! ***************************************************************
! *******************     X-SWEEP       *************************
! ***************************************************************

  case(DIR_X)

! --- Computing a unit vector along a magnetic field line ---
! --- bx, by, bz are actually the directional cosines     ---

  By_mean = 0.25*( faceyData(MAG_FACE_VAR,ix,iy,iz) + &
                   faceyData(MAG_FACE_VAR,ix,iy+1,iz) + &
                   faceyData(MAG_FACE_VAR,ix-1,iy,iz) + &
                   faceyData(MAG_FACE_VAR,ix-1,iy+1,iz) )

  Bz_mean = 0.25*( facezData(MAG_FACE_VAR,ix,iy,iz) + &
                   facezData(MAG_FACE_VAR,ix,iy,iz+1) + &
                   facezData(MAG_FACE_VAR,ix-1,iy,iz) + &
                   facezData(MAG_FACE_VAR,ix-1,iy,iz+1) )

  modB = facexData(MAG_FACE_VAR,ix,iy,iz)**2.0 +  & 
         By_mean**2.0 + Bz_mean**2.0
  modB = sqrt(modB)

  bxx = facexData(MAG_FACE_VAR,ix,iy,iz)/modB
  bxy = By_mean/modB
  bxz = Bz_mean/modB


  if (isotropic) then
     bxx = 1.0
     bxy = 0.0
     bxz = 0.0
  endif 


! --- (unsaturated) conductive flux calculation ---

  al = U(TEMP_VAR,ix-1,iy,iz)-U(TEMP_VAR,ix-1,iy-1,iz)
  bl = U(TEMP_VAR,ix-1,iy+1,iz)-U(TEMP_VAR,ix-1,iy,iz)

  ar = U(TEMP_VAR,ix,iy,iz)-U(TEMP_VAR,ix,iy-1,iz)
  br = U(TEMP_VAR,ix,iy+1,iz)-U(TEMP_VAR,ix,iy,iz)

  call MC_limiter(al,bl,dtdyxl)
  call MC_limiter(ar,br,dtdyxr)
  call MC_limiter(dtdyxl,dtdyxr,dtdyx)
!  call minmod_local(al,bl,dtdyxl)
!  call minmod_local(ar,br,dtdyxr)
!  call minmod_local(dtdyxl,dtdyxr,dtdyx)

  al = U(TEMP_VAR,ix-1,iy,iz)-U(TEMP_VAR,ix-1,iy,iz-1)
  bl = U(TEMP_VAR,ix-1,iy,iz+1)-U(TEMP_VAR,ix-1,iy,iz)

  ar = U(TEMP_VAR,ix,iy,iz)-U(TEMP_VAR,ix,iy,iz-1)
  br = U(TEMP_VAR,ix,iy,iz+1)-U(TEMP_VAR,ix,iy,iz)

  call MC_limiter(al,bl,dtdzxl)
  call MC_limiter(ar,br,dtdzxr)
  call MC_limiter(dtdzxl,dtdzxr,dtdzx)
!  call minmod_local(al,bl,dtdzxl)
!  call minmod_local(ar,br,dtdzxr)
!  call minmod_local(dtdzxl,dtdzxr,dtdzx)

  
  ! projecting (unsaturated) flux onto the magnetic field line
  !  (this is for the x-sweep ONLY)

  FthBx = -bxx*idx*(U(TEMP_VAR,ix,iy,iz)-U(TEMP_VAR,ix-1,iy,iz)) &
          -bxy*idy*dtdyx -bxz*idz*dtdzx 
  if (cond(ix,iy,iz)*cond(ix-1,iy,iz) .eq. 0.0) then
     cond_eff = 0.0
  else
     cond_eff = 2.0*cond(ix,iy,iz)*cond(ix-1,iy,iz) / (cond(ix,iy,iz)+cond(ix-1,iy,iz))
  endif
  !cond_eff = cond_eff / mhd_qref
  FthBx = cond_eff * FthBx


if (sim_saturatedConduction) then

  qsata = 1.03e13 * (U(TEMP_VAR,ix,iy,iz)**1.5) * U(DENS_VAR,ix,iy,iz)
  qsatb = 1.03e13 * (U(TEMP_VAR,ix-1,iy,iz)**1.5) * U(DENS_VAR,ix-1,iy,iz)
  qsat = 2.0*(qsata*qsatb)/(qsata+qsatb) * SpitzerFraction

  gTx =        U(TEMP_VAR,ix  ,iy,  iz  )-U(TEMP_VAR,ix-1,iy  ,iz  )
  gTy = 0.25*( U(TEMP_VAR,ix-1,iy,  iz  )-U(TEMP_VAR,ix-1,iy-1,iz  ) + &
               U(TEMP_VAR,ix-1,iy+1,iz  )-U(TEMP_VAR,ix-1,iy,  iz  ) + &
               U(TEMP_VAR,ix,  iy,  iz  )-U(TEMP_VAR,ix,  iy-1,iz  ) + &
               U(TEMP_VAR,ix,  iy+1,iz  )-U(TEMP_VAR,ix,  iy,  iz  ) )
  gTz = 0.25*( U(TEMP_VAR,ix-1,iy,  iz  )-U(TEMP_VAR,ix-1,iy,  iz-1) + &
               U(TEMP_VAR,ix-1,iy,  iz+1)-U(TEMP_VAR,ix-1,iy,  iz  ) + &
               U(TEMP_VAR,ix,  iy,  iz  )-U(TEMP_VAR,ix,  iy,  iz-1) + &
               U(TEMP_VAR,ix,  iy,  iz+1)-U(TEMP_VAR,ix,  iy,  iz  )       )

  gr = bxx*gTx + bxy*gTy + bxz*gTz
  mod_gT = sqrt( gTx*gTx + gTy*gTy + gTz*gTz )
  F_sat = - gr/mod_gT * qsat

  ! to be conservative we want to activate saturated fluxes as soon as possible
  ! this means that we want to use the arithmetic average as harmonic average 
  ! would naturally weigh lambda_mfp towards smaller values for which the flux would be 
  ! in the unsaturated regime  
  lambda_mfp = 1.19e-17 * 0.5* ( (U(TEMP_VAR,ix,iy,iz)**2.0) / U(DENS_VAR,ix,iy,iz) + &
				 (U(TEMP_VAR,ix-1,iy,iz)**2.0) / U(DENS_VAR,ix-1,iy,iz) )
  gr = abs(idx*gr)

  xt = 1.0 / ( 4.2*lambda_mfp*gr/(0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix-1,iy,iz))) )
  ft = tanh(xt - (1.0/xt))
  ft = 0.5*(1.0+ft)


!  qq = del(1)*gr/(min(U(TEMP_VAR,ix,iy,iz),U(TEMP_VAR,ix-1,iy,iz)))
!  if (qq .gt. 2.0) then ! huge gradient over one zone
!     qq=1.0
!  else
!     qq=0.0
!  endif


!  if ((U(DENS_VAR,ix,iy,iz)/1.67e-24) .le. 1.0e-6 .or. (U(DENS_VAR,ix-1,iy,iz)/1.67e-24) .le. 1.0e-6) F_sat = 0.0
!  if (((U(DENS_VAR,ix,iy,iz)/1.67e-24) .le. 1.0e-7 .or. & 
!       (U(DENS_VAR,ix-1,iy,iz)/1.67e-24) .le. 1.0e-7) .or. currentRedshift .gt. 4.0) F_sat = 0.0

!  if (((U(DENS_VAR,ix,iy,iz)/1.67e-24) .le. 1.0e-7 .or. & 
!       (U(DENS_VAR,ix-1,iy,iz)/1.67e-24) .le. 1.0e-7)) F_sat = 0.0
!  if ( U(TEMP_VAR,ix,iy,iz) .le. 1.0e3 .or. U(TEMP_VAR,ix-1,iy,iz) .le. 1.0e3) F_sat = 0.0

  FthBx = FthBx*ft + F_sat*(1.0-ft)


!  FthBx = FthBx/( 1.0 + 1/xt )

!  FthBx = FthBx/( 1.0 + qq/xt )

!  FthBx = -   min(1.0e8, 0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix-1,iy,iz)))*cond_eff* gr*del(1)/mod_gT / & 
!            ( min(1.0e8, 0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix-1,iy,iz)))/gr + 4.2*lambda_mfp)
!  if (  0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix-1,iy,iz)) .gt. 1.0e6 ) FthBx = 0.0
!

!  FthBx = FthBx/(1.0+del(1)*gr/(0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix-1,iy,iz))))
!  FthBx = FthBx/(1.0+del(1)*gr/(min(U(TEMP_VAR,ix,iy,iz),U(TEMP_VAR,ix-1,iy,iz))))

!  if (  0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix-1,iy,iz))/gr .lt. 1.0e-2 * 4.2*lambda_mfp) then
!     write(*,*)'aaaa', 0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix-1,iy,iz))/gr, 4.2*lambda_mfp, 0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix-1,iy,iz))
!     stop
!  endif

  Dcff  = cond(ix,iy,iz)/(U(DENS_VAR,ix,iy,iz)*cv)
!!print*,'in hy_uhd_addThermalFluxes 1, Dcff=',Dcff,U(DENS_VAR,ix,iy,iz),cv
!  if (Dcff .gt. 1.0001e32) then
!     write(*,*)'Dcff too large in x-sweep', Dcff, cond(ix,iy,iz), U(DENS_VAR,ix,iy,iz), U(TEMP_VAR,ix,iy,iz) 
!     stop
!  endif

  Dcff  =  Dcff*ft + (cs*del(1)/0.4)*(1.0-ft)
!!print*,'in hy_uhd_addThermalFluxes 2, Dcff=',Dcff,ft,cs,del(1)
!   Dcff  =  Dcff*ft + 4.0*(cs*del(1)/0.4)*(1.0-ft)
! if (U(DENS_VAR,ix,iy,iz) .le. 1.0e-6 .or. U(DENS_VAR,ix-1,iy,iz) .le. 1.0e-6) Dcff = 1.0e-20
!  if (((U(DENS_VAR,ix,iy,iz)/1.67e-24) .le. 1.0e-7 .or. & 
!       (U(DENS_VAR,ix-1,iy,iz)/1.67e-24) .le. 1.0e-7) .or. currentRedshift .gt. 4.0) Dcff = 1.0e-20
!  if ( U(TEMP_VAR,ix,iy,iz) .le. 1.0e3 .or. U(TEMP_VAR,ix-1,iy,iz) .le. 1.0e3) Dcff = 1.0e-20

!  if (((U(DENS_VAR,ix,iy,iz)/1.67e-24) .le. 1.0e-7 .or. & 
!       (U(DENS_VAR,ix-1,iy,iz)/1.67e-24) .le. 1.0e-7)) Dcff = 1.0e-20
!  if ( U(TEMP_VAR,ix,iy,iz) .le. 1.0e3 .or. U(TEMP_VAR,ix-1,iy,iz) .le. 1.0e3) Dcff = 1.0e-20

endif



! ***************************************************************
! *******************     Y-SWEEP       *************************
! ***************************************************************

  case(DIR_Y)

! --- Computing a unit vector along a magnetic field line ---
! --- bx, by, bz are actually the directional cosines     ---

  Bx_mean = 0.25*( facexData(MAG_FACE_VAR,ix,iy,iz) + &
                   facexData(MAG_FACE_VAR,ix+1,iy,iz) + &
                   facexData(MAG_FACE_VAR,ix,iy-1,iz) + &
                   facexData(MAG_FACE_VAR,ix+1,iy-1,iz) )

  Bz_mean = 0.25*( facezData(MAG_FACE_VAR,ix,iy,iz) + &
                   facezData(MAG_FACE_VAR,ix,iy,iz+1) + &
                   facezData(MAG_FACE_VAR,ix,iy-1,iz) + &
                   facezData(MAG_FACE_VAR,ix,iy-1,iz+1) )

  modB = faceyData(MAG_FACE_VAR,ix,iy,iz)**2.0 +  & 
         Bx_mean**2.0 + Bz_mean**2.0
  modB = sqrt(modB)

  byy = faceyData(MAG_FACE_VAR,ix,iy,iz)/modB
  byx = Bx_mean/modB
  byz = Bz_mean/modB


  if (isotropic) then
     byy = 1.0
     byx = 0.0
     byz = 0.0
  endif


! --- (unsaturated) conductive flux calculation ---

  al = U(TEMP_VAR,ix,iy-1,iz)-U(TEMP_VAR,ix-1,iy-1,iz)
  bl = U(TEMP_VAR,ix+1,iy-1,iz)-U(TEMP_VAR,ix,iy-1,iz)

  ar = U(TEMP_VAR,ix,iy,iz)-U(TEMP_VAR,ix-1,iy,iz)
  br = U(TEMP_VAR,ix+1,iy,iz)-U(TEMP_VAR,ix,iy,iz)

  call MC_limiter(al,bl,dtdxyl)
  call MC_limiter(ar,br,dtdxyr)
  call MC_limiter(dtdxyl,dtdxyr,dtdxy)
!  call minmod_local(al,bl,dtdxyl)
!  call minmod_local(ar,br,dtdxyr)
!  call minmod_local(dtdxyl,dtdxyr,dtdxy)

  al = U(TEMP_VAR,ix,iy-1,iz)-U(TEMP_VAR,ix,iy-1,iz-1)
  bl = U(TEMP_VAR,ix,iy-1,iz+1)-U(TEMP_VAR,ix,iy-1,iz)

  ar = U(TEMP_VAR,ix,iy,iz)-U(TEMP_VAR,ix,iy,iz-1)
  br = U(TEMP_VAR,ix,iy,iz+1)-U(TEMP_VAR,ix,iy,iz)

  call MC_limiter(al,bl,dtdzyl)
  call MC_limiter(ar,br,dtdzyr)
  call MC_limiter(dtdzyl,dtdzyr,dtdzy)
!  call minmod_local(al,bl,dtdzyl)
!  call minmod_local(ar,br,dtdzyr)
!  call minmod_local(dtdzyl,dtdzyr,dtdzy)

  
  ! projecting (unsaturated) flux onto the magnetic field line
  !  (this is for the y-sweep ONLY)

  FthBy = -byy*idy*(U(TEMP_VAR,ix,iy,iz)-U(TEMP_VAR,ix,iy-1,iz)) &
          -byx*idx*dtdxy -byz*idz*dtdzy 
  if (cond(ix,iy,iz)*cond(ix,iy-1,iz) .eq. 0.0) then
     cond_eff = 0.0
  else
     cond_eff = 2.0*cond(ix,iy,iz)*cond(ix,iy-1,iz) / (cond(ix,iy,iz)+cond(ix,iy-1,iz))
  endif
  !cond_eff = cond_eff / mhd_qref
  FthBy = cond_eff * FthBy


if (sim_saturatedConduction) then

  qsata = 1.03e13 * (U(TEMP_VAR,ix,iy,iz)**1.5) * U(DENS_VAR,ix,iy,iz)
  qsatb = 1.03e13 * (U(TEMP_VAR,ix,iy-1,iz)**1.5) * U(DENS_VAR,ix,iy-1,iz)
  qsat = 2.0*(qsata*qsatb)/(qsata+qsatb) * SpitzerFraction

  gTy =        U(TEMP_VAR,ix  ,iy,  iz  )-U(TEMP_VAR,ix  ,iy-1,iz  )
  gTx = 0.25*( U(TEMP_VAR,ix  ,iy-1,iz  )-U(TEMP_VAR,ix-1,iy-1,iz  ) + &
               U(TEMP_VAR,ix+1,iy-1,iz  )-U(TEMP_VAR,ix,  iy-1,iz  ) + &         
               U(TEMP_VAR,ix  ,iy,  iz  )-U(TEMP_VAR,ix-1,iy,  iz  ) + &
               U(TEMP_VAR,ix+1,iy,  iz  )-U(TEMP_VAR,ix,  iy,  iz  ) )
  gTz = 0.25*( U(TEMP_VAR,ix,  iy-1,iz  )-U(TEMP_VAR,ix,  iy-1,iz-1) + &
               U(TEMP_VAR,ix,  iy-1,iz+1)-U(TEMP_VAR,ix,  iy-1,iz  ) + &
               U(TEMP_VAR,ix,  iy,  iz  )-U(TEMP_VAR,ix,  iy,  iz-1) + &
               U(TEMP_VAR,ix,  iy,  iz+1)-U(TEMP_VAR,ix,  iy,  iz  )       )

  gr = byx*gTx + byy*gTy + byz*gTz
  mod_gT = sqrt( gTx*gTx + gTy*gTy + gTz*gTz )
  F_sat = - gr/mod_gT * qsat

  ! to be conservative we want to activate saturated fluxes as soon as possible
  ! this means that we want to use the arithmetic average as harmonic average 
  ! would naturally weigh lambda_mfp towards smaller values for which the flux would be 
  ! in the unsaturated regime  
  lambda_mfp = 1.19e-17 * 0.5* ( (U(TEMP_VAR,ix,iy,iz)**2.0) / U(DENS_VAR,ix,iy,iz) + &
				 (U(TEMP_VAR,ix,iy-1,iz)**2.0) / U(DENS_VAR,ix,iy-1,iz) )
  gr = abs(idy*gr)

  xt = 1.0 / ( 4.2*lambda_mfp*gr/(0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix,iy-1,iz))) )
  ft = tanh(xt - (1.0/xt))
  ft = 0.5*(1.0+ft)


!  qq = del(1)*gr/(min(U(TEMP_VAR,ix,iy,iz),U(TEMP_VAR,ix,iy-1,iz)))
!  if (qq .gt. 2.0) then ! huge gradient over one zone
!     qq=1.0
!  else
!     qq=0.0
!  endif


!  if ((U(DENS_VAR,ix,iy,iz)/1.67e-24) .le. 1.0e-6 .or. (U(DENS_VAR,ix,iy-1,iz)/1.67e-24) .le. 1.0e-6) F_sat = 0.0
!  if (((U(DENS_VAR,ix,iy,iz)/1.67e-24) .le. 1.0e-7 .or. & 
!       (U(DENS_VAR,ix,iy-1,iz)/1.67e-24) .le. 1.0e-7) .or. currentRedshift .gt. 4.0) F_sat = 0.0
!  if ( U(TEMP_VAR,ix,iy,iz) .le. 1.0e3 .or. U(TEMP_VAR,ix,iy-1,iz) .le. 1.0e3) F_sat = 0.0

!  if (((U(DENS_VAR,ix,iy,iz)/1.67e-24) .le. 1.0e-7 .or. & 
!       (U(DENS_VAR,ix,iy-1,iz)/1.67e-24) .le. 1.0e-7)) F_sat = 0.0
!  if ( U(TEMP_VAR,ix,iy,iz) .le. 1.0e3 .or. U(TEMP_VAR,ix,iy-1,iz) .le. 1.0e3) F_sat = 0.0

  FthBy = FthBy*ft + F_sat*(1.0-ft)


!  FthBy = FthBy/( 1.0 + 1/xt )
!  FthBy = FthBy/( 1.0 + qq/xt )

!  FthBy = -   min(1.0e8, 0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix,iy-1,iz)))*cond_eff* gr*del(1)/mod_gT / & 
!            ( min(1.0e8, 0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix,iy-1,iz)))/gr + 4.2*lambda_mfp)
!  if (  0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix,iy-1,iz)) .gt. 1.0e6 ) FthBy = 0.0
!

!  FthBy = FthBy/(1.0+del(1)*gr/(0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix,iy-1,iz))))
!  FthBy = FthBy/(1.0+del(1)*gr/(min(U(TEMP_VAR,ix,iy,iz),U(TEMP_VAR,ix,iy-1,iz))))

  Dcff  = cond(ix,iy,iz)/(U(DENS_VAR,ix,iy,iz)*cv)

!  if (Dcff .gt. 1.0001e32) then
!     write(*,*)'Dcff too large in y-sweep', Dcff, cond(ix,iy,iz), U(DENS_VAR,ix,iy,iz)/1.67e-24, cv
!     stop
!  endif

  Dcff  =  Dcff*ft + (cs*del(1)/0.4)*(1.0-ft)
!  Dcff  =  Dcff*ft + 4.0*(cs*del(1)/0.4)*(1.0-ft)

!  if (U(DENS_VAR,ix,iy,iz) .le. 1.0e-6 .or. U(DENS_VAR,ix,iy-1,iz) .le. 1.0e-6) Dcff = 1.0e-20
!  if (((U(DENS_VAR,ix,iy,iz)/1.67e-24) .le. 1.0e-7 .or. & 
!       (U(DENS_VAR,ix,iy-1,iz)/1.67e-24) .le. 1.0e-7) .or. currentRedshift .gt. 4.0) Dcff = 1.0e-20
!  if ( U(TEMP_VAR,ix,iy,iz) .le. 1.0e3 .or. U(TEMP_VAR,ix,iy-1,iz) .le. 1.0e3) Dcff = 1.0e-20

!  if (((U(DENS_VAR,ix,iy,iz)/1.67e-24) .le. 1.0e-7 .or. & 
!       (U(DENS_VAR,ix,iy-1,iz)/1.67e-24) .le. 1.0e-7)) Dcff = 1.0e-20
!  if ( U(TEMP_VAR,ix,iy,iz) .le. 1.0e3 .or. U(TEMP_VAR,ix,iy-1,iz) .le. 1.0e3) Dcff = 1.0e-20


endif


! ***************************************************************
! *******************     Z-SWEEP       *************************
! ***************************************************************

  case(DIR_Z)

! --- Computing a unit vector along a magnetic field line ---
! --- bx, by, bz are actually the directional cosines     ---

  Bx_mean = 0.25*( facexData(MAG_FACE_VAR,ix,iy,iz) + &
                   facexData(MAG_FACE_VAR,ix+1,iy,iz) + &
                   facexData(MAG_FACE_VAR,ix,iy,iz-1) + &
                   facexData(MAG_FACE_VAR,ix+1,iy,iz-1) )

  By_mean = 0.25*( faceyData(MAG_FACE_VAR,ix,iy,iz) + &
                   faceyData(MAG_FACE_VAR,ix,iy+1,iz) + &
                   faceyData(MAG_FACE_VAR,ix,iy,iz-1) + &
                   faceyData(MAG_FACE_VAR,ix,iy+1,iz-1) )

  modB = facezData(MAG_FACE_VAR,ix,iy,iz)**2.0 +  & 
         By_mean**2.0 + Bx_mean**2.0
  modB = sqrt(modB)

  bzz = facezData(MAG_FACE_VAR,ix,iy,iz)/modB
  bzx = Bx_mean/modB
  bzy = By_mean/modB


  if (isotropic) then
     bzz = 1.0
     bzx = 0.0
     bzy = 0.0
  endif


! --- (unsaturated) conductive flux calculation --- 

  al = U(TEMP_VAR,ix,iy,iz-1)-U(TEMP_VAR,ix-1,iy,iz-1)
  bl = U(TEMP_VAR,ix+1,iy,iz-1)-U(TEMP_VAR,ix,iy,iz-1)

  ar = U(TEMP_VAR,ix,iy,iz)-U(TEMP_VAR,ix-1,iy,iz)
  br = U(TEMP_VAR,ix+1,iy,iz)-U(TEMP_VAR,ix,iy,iz)

  call MC_limiter(al,bl,dtdxzl)
  call MC_limiter(ar,br,dtdxzr)
  call MC_limiter(dtdxzl,dtdxzr,dtdxz)
!  call minmod_local(al,bl,dtdxzl)
!  call minmod_local(ar,br,dtdxzr)
!  call minmod_local(dtdxzl,dtdxzr,dtdxz)

  al = U(TEMP_VAR,ix,iy,iz-1)-U(TEMP_VAR,ix,iy-1,iz-1)
  bl = U(TEMP_VAR,ix,iy+1,iz-1)-U(TEMP_VAR,ix,iy,iz-1)

  ar = U(TEMP_VAR,ix,iy,iz)-U(TEMP_VAR,ix,iy-1,iz)
  br = U(TEMP_VAR,ix,iy+1,iz)-U(TEMP_VAR,ix,iy,iz)

  call MC_limiter(al,bl,dtdyzl)
  call MC_limiter(ar,br,dtdyzr)
  call MC_limiter(dtdyzl,dtdyzr,dtdyz)
!  call minmod_local(al,bl,dtdyzl)
!  call minmod_local(ar,br,dtdyzr)
!  call minmod_local(dtdyzl,dtdyzr,dtdyz)

  
  ! projecting (unsaturated) flux onto the magnetic field line
  !  (this is for the z-sweep ONLY)

  FthBz = -bzz*idz*(U(TEMP_VAR,ix,iy,iz)-U(TEMP_VAR,ix,iy,iz-1)) &
          -bzx*idx*dtdxz -bzy*idy*dtdyz 
  if (cond(ix,iy,iz)*cond(ix,iy,iz-1) .eq. 0.0) then
     cond_eff = 0.0
  else
     cond_eff = 2.0*cond(ix,iy,iz)*cond(ix,iy,iz-1) / (cond(ix,iy,iz)+cond(ix,iy,iz-1))
  endif
  !cond_eff = cond_eff / mhd_qref
  FthBz = cond_eff * FthBz



if (sim_saturatedConduction) then

  qsata = 1.03e13 * (U(TEMP_VAR,ix,iy,iz)**1.5) * U(DENS_VAR,ix,iy,iz)
  qsatb = 1.03e13 * (U(TEMP_VAR,ix,iy,iz-1)**1.5) * U(DENS_VAR,ix,iy,iz-1)
  qsat = 2.0*(qsata*qsatb)/(qsata+qsatb) * SpitzerFraction

  gTz =        U(TEMP_VAR,ix  ,iy,  iz  )-U(TEMP_VAR,ix  ,iy,  iz-1)
  gTx = 0.25*( U(TEMP_VAR,ix,  iy,  iz-1)-U(TEMP_VAR,ix-1,iy,  iz-1) + &
               U(TEMP_VAR,ix+1,iy,  iz-1)-U(TEMP_VAR,ix,  iy,  iz-1) + &
               U(TEMP_VAR,ix,  iy,  iz  )-U(TEMP_VAR,ix-1,iy,  iz  ) + &
               U(TEMP_VAR,ix+1,iy,  iz  )-U(TEMP_VAR,ix,  iy,  iz  )       )
  gTy = 0.25*( U(TEMP_VAR,ix  ,iy,  iz-1)-U(TEMP_VAR,ix,  iy-1,iz-1) + &
               U(TEMP_VAR,ix,  iy+1,iz-1)-U(TEMP_VAR,ix,  iy,  iz-1) + &         
               U(TEMP_VAR,ix  ,iy,  iz  )-U(TEMP_VAR,ix,  iy-1,iz  ) + &
               U(TEMP_VAR,ix,  iy+1,iz  )-U(TEMP_VAR,ix,  iy,  iz  ) )

  gr = bzx*gTx + bzy*gTy + bzz*gTz
  mod_gT = sqrt( gTx*gTx + gTy*gTy + gTz*gTz )
  F_sat = - gr/mod_gT * qsat

  ! to be conservative we want to activate saturated fluxes as soon as possible
  ! this means that we want to use the arithmetic average as harmonic average 
  ! would naturally weigh lambda_mfp towards smaller values for which the flux would be 
  ! in the unsaturated regime  
  lambda_mfp = 1.19e-17 * 0.5* ( (U(TEMP_VAR,ix,iy,iz)**2.0) / U(DENS_VAR,ix,iy,iz) + &
				 (U(TEMP_VAR,ix,iy,iz-1)**2.0) / U(DENS_VAR,ix,iy,iz-1) )
  gr = abs(idz*gr)

  xt = 1.0 / ( 4.2*lambda_mfp*gr/(0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix,iy,iz-1))) )
  ft = tanh(xt - (1.0/xt))
  ft = 0.5*(1.0+ft)


!  qq = del(1)*gr/(min(U(TEMP_VAR,ix,iy,iz),U(TEMP_VAR,ix,iy,iz-1)))
!  if (qq .gt. 2.0) then ! huge gradient over one zone
!     qq=1.0
!  else
!     qq=0.0
!  endif


!  if ((U(DENS_VAR,ix,iy,iz)/1.67e-24) .le. 1.0e-6 .or. (U(DENS_VAR,ix,iy,iz-1)/1.67e-24) .le. 1.0e-6) F_sat = 0.0
!  if (((U(DENS_VAR,ix,iy,iz)/1.67e-24) .le. 1.0e-7 .or. & 
!       (U(DENS_VAR,ix,iy,iz-1)/1.67e-24) .le. 1.0e-7) .or. currentRedshift .gt. 4.0) F_sat = 0.0
!  if ( U(TEMP_VAR,ix,iy,iz) .le. 1.0e3 .or. U(TEMP_VAR,ix,iy,iz-1) .le. 1.0e3) F_sat = 0.0

!  if (((U(DENS_VAR,ix,iy,iz)/1.67e-24) .le. 1.0e-7 .or. & 
!       (U(DENS_VAR,ix,iy,iz-1)/1.67e-24) .le. 1.0e-7)) F_sat = 0.0
!  if ( U(TEMP_VAR,ix,iy,iz) .le. 1.0e3 .or. U(TEMP_VAR,ix,iy,iz-1) .le. 1.0e3) F_sat = 0.0

  FthBz = FthBz*ft + F_sat*(1.0-ft)


! FthBz = FthBz/( 1.0 + 1/xt )
!  FthBz = FthBz/( 1.0 + qq/xt )

!  FthBz = -   min(1.0e8, 0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix,iy,iz-1)))*cond_eff* gr*del(1)/mod_gT / & 
!            ( min(1.0e8, 0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix,iy,iz-1)))/gr + 4.2*lambda_mfp)
!  if (  0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix,iy,iz-1)) .gt. 1.0e6 ) FthBz = 0.0
!

!  FthBz = FthBz/(1.0+del(1)*gr/(0.5*(U(TEMP_VAR,ix,iy,iz)+U(TEMP_VAR,ix,iy,iz-1))))
!  FthBz = FthBz/(1.0+del(1)*gr/(min(U(TEMP_VAR,ix,iy,iz),U(TEMP_VAR,ix,iy,iz-1))))

  Dcff  = cond(ix,iy,iz)/(U(DENS_VAR,ix,iy,iz)*cv)

!  if (Dcff .gt. 1.0001e32) then
!     write(*,*)'Dcff too large in z-sweep', Dcff, cond(ix,iy,iz), U(DENS_VAR,ix,iy,iz)/1.67e-24, cv
!     stop
!  endif

  Dcff  =  Dcff*ft + (cs*del(1)/0.4)*(1.0-ft)
!  Dcff  =  Dcff*ft + 4.0*(cs*del(1)/0.4)*(1.0-ft)

!  if (U(DENS_VAR,ix,iy,iz) .le. 1.0e-6 .or. U(DENS_VAR,ix,iy,iz-1) .le. 1.0e-6) Dcff = 1.0e-20
!  if (((U(DENS_VAR,ix,iy,iz)/1.67e-24) .le. 1.0e-7 .or. & 
!       (U(DENS_VAR,ix,iy,iz-1)/1.67e-24) .le. 1.0e-7) .or. currentRedshift .gt. 4.0) Dcff = 1.0e-20
!  if ( U(TEMP_VAR,ix,iy,iz) .le. 1.0e3 .or. U(TEMP_VAR,ix,iy,iz-1) .le. 1.0e3) Dcff = 1.0e-20

!  if (((U(DENS_VAR,ix,iy,iz)/1.67e-24) .le. 1.0e-7 .or. & 
!       (U(DENS_VAR,ix,iy,iz-1)/1.67e-24) .le. 1.0e-7)) Dcff = 1.0e-20
!  if ( U(TEMP_VAR,ix,iy,iz) .le. 1.0e3 .or. U(TEMP_VAR,ix,iy,iz-1) .le. 1.0e3) Dcff = 1.0e-20

endif



end select


! --- projecting the conductive flux along the field line onto the sweep directions

  select case(sweepDir)
  case(DIR_X)

!      if (dist >= 0.9*cent)  FthBx = 0.0

      Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX)    + FthBx*bxx 

#if NDIM >= 2
  case(DIR_Y)

!      if (dist >= 0.9*cent)  FthBy = 0.0

      Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX)    + FthBy*byy 

#endif

#if NDIM == 3
  case(DIR_Z)

!      if (dist >= 0.9*cent)  FthBz = 0.0

      Flux(F05ENER_FLUX) = Flux(F05ENER_FLUX)    + FthBz*bzz  

#endif

  end select


  call Grid_releaseBlkPtr(blockID,U,CENTER)
!!print*,'in hy_uhd_addThermalFluxes, Dcff in dir=',Dcff,sweepDir
end subroutine hy_uhd_addThermalFluxes




subroutine MC_limiter(a,b,c)

real :: a,b,c,ab_tmp1,ab_tmp2,ab_tmp3

call minmod_local(a,b,ab_tmp1)
ab_tmp2 = 2.0e0*ab_tmp1
ab_tmp3 = 0.5e0*(a+b)
call minmod_local(ab_tmp2,ab_tmp3,c)

end subroutine MC_limiter



subroutine minmod_local(a,b,c)

real :: a,b,c

if (a > 0.0e0 .and. b > 0.0e0) c = min(a,b)
if (a < 0.0e0 .and. b < 0.0e0) c = max(a,b)
if (a*b <= 0.0e0)              c = 0.0

end subroutine minmod_local
