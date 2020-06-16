!!****if* source/physics/Hydro/HydroMain/unsplit/MHD_StaggeredMesh/hy_uhd_getFaceFlux
!!
!! NAME
!!
!!  hy_uhd_getFaceFlux
!!
!! SYNOPSIS
!!
!!  hy_uhd_getFaceFlux( integer(IN) :: blockID,
!!                      integer(IN) :: blkLimits(2,MDIM)
!!                      integer(IN) :: blkLimitsGC(2,MDIM), 
!!                      integer(IN) :: datasize(3),
!!                      real(OUT)   :: xflux(:,:,:,:),
!!                      real(OUT)   :: yflux(:,:,:,:),
!!                      real(OUT)   :: zflux(:,:,:,:),
!!                      real(OUT)   :: vint(:,:,:,:) )
!!
!! ARGUMENTS
!!
!!  blockID           - a current block ID
!!  blkLimits         - an array that holds the lower and upper indices of the section
!!                      of block without the guard cells
!!  blkLimitsGC       - an array that holds the lower and upper indices of the section
!!                      of block with the guard cells
!!  datasize          - data size for boundary extrapolated data, bdryData
!!  xflux,yflux,zflux - face fluxes at each {x,y,z} direction
!!  vint              - interface velocity
!!
!! DESCRIPTION
!!
!!  This routine computes high-order Godunov fluxes at cell interface centers 
!!  for each spatial direction using a choice of Riemann solvers.
!!  Choices of Riemann solvers are Roe-type, HLL(E), HLLC, HLLD, and Lax-Friedrichs solvers.
!!
!!*** 

Subroutine hy_uhd_getFaceFlux ( blockID,blkLimits,blkLimitsGC, &
                                datasize,xflux,yflux,zflux,vint)

  use Hydro_data,                    ONLY : hy_nref,hy_kref,hy_mref,&
                                            hy_useResistivity,      &
                                            hy_maxMagDiff,          &
                                            hy_shockInstabilityFix, &
                                            hy_RiemannSolver,       &
                                            hy_useDiffuse,          &
                                            hy_useViscosity,        &
                                            hy_useConductivity,     &
                                            hy_useMagneticResistivity

  use hy_uhd_interface,              ONLY : hy_uhd_addViscousFluxes,  &
                                            hy_uhd_addThermalFluxes,  &
                                            hy_uhd_addResistiveFluxes,&
                                            hy_uhd_Roe, &
                                            hy_uhd_LF,  &
                                            hy_uhd_HLL, &
                                            hy_uhd_HLLC,&
                                            hy_uhd_HLLD,&
                                            hy_uhd_HLL7,&
                                            hy_uhd_shockInstabilityFix
  use Grid_interface,                ONLY : Grid_getBlkPtr, Grid_releaseBlkPtr
  use Eos_interface,                 ONLY : Eos
  use Conductivity_interface,        ONLY : Conductivity
  use Viscosity_interface,           ONLY : Viscosity
  use MagneticResistivity_interface, ONLY : MagneticResistivity
  use Simulation_data,               ONLY : sim_saturatedConduction

  implicit none

#include "constants.h"
#include "Flash.h"
#include "Eos.h"
#include "UHD.h"

  !! Arguments type declaration ------------------------------
  integer, intent(IN)  :: blockID
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM),     intent(IN) :: datasize

#ifdef FIXEDBLOCKSIZE
  real, DIMENSION(NFLUXES,      &
       GRID_ILO_GC:GRID_IHI_GC, &
       GRID_JLO_GC:GRID_JHI_GC, &
       GRID_KLO_GC:GRID_KHI_GC),&
       intent(OUT) :: xflux, yflux, zflux
  real, DIMENSION(MDIM,         &
       GRID_ILO_GC:GRID_IHI_GC, &
       GRID_JLO_GC:GRID_JHI_GC, &
       GRID_KLO_GC:GRID_KHI_GC),&
       intent(OUT) :: vint
#else
  real, DIMENSION(NFLUXES,                             &
       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
       intent(OUT) :: xflux, yflux, zflux
  real, DIMENSION(MDIM,                                &
       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
       intent(OUT) :: vint
#endif
  !! ---------------------------------------------------------

  integer :: i0,imax,j0,jmax,k0,kmax,jbeg,jend,kbeg,kend, i,j,k
  real, pointer, dimension(:,:,:,:) :: scratch, U
  real, dimension(HY_VARINUM+4) :: VL,VR 
  !(DENS,VELX,VELY,VELZ,MAGX,MAGY,MAGZ,PRES + GAMC,GAME,EINT,TEMP)

!#ifdef MAGNETIC_RESISTIVITY
  real, dimension(NSPECIES) :: speciesArr
  real :: maxMagResist,maxKinVisc,maxCond,maxDcff
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, &
                  GRID_JLO_GC:GRID_JHI_GC, &
                  GRID_KLO_GC:GRID_KHI_GC) &
                  :: magResist,viscDynamic,viscKinematic,cond,dcff,&
	             Dcffx, Dcffy, Dcffz

#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) &
                  :: magResist,viscDynamic,viscKinematic,cond,dcff,&
	             Dcffx, Dcffy, Dcffz
#endif
!#endif

#ifdef FIXEDBLOCKSIZE
  real, dimension(MDIM,2,&
       GRID_ILO_GC:GRID_IHI_GC, &
       GRID_JLO_GC:GRID_JHI_GC, &
       GRID_KLO_GC:GRID_KHI_GC) :: LambdaFix
#else
  real, dimension(MDIM,2, &
       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)) &
       :: LambdaFix
#endif
  real,dimension(EOS_NUM)::eosData
  logical,dimension(EOS_VARS+1:EOS_NUM)::eosMask=.false.
  real, dimension(NSPECIES) :: eosMf
  integer :: vecLen
  integer :: interp_eosMode = MODE_DENS_PRES

  real :: maxDcffx, maxDcffy, maxDcffz


  !! initialize
  vecLen=1
  eosMf =1.

  !! Call scratch for storing fluxes
  call Grid_getBlkPtr(blockID,scratch,SCRATCH)
  call Grid_getBlkPtr(blockID,U,CENTER)

  i0   = blkLimits(LOW, IAXIS)
  imax = blkLimits(HIGH,IAXIS)
  j0   = blkLimits(LOW, JAXIS)
  jmax = blkLimits(HIGH,JAXIS)
  k0   = blkLimits(LOW, KAXIS)
  kmax = blkLimits(HIGH,KAXIS)

  if (NDIM == 1) then
     jbeg= 4
     jend=-2
     kbeg= 4
     kend=-2
  elseif (NDIM == 2) then
     jbeg= j0
     jend= jmax
     kbeg= 4
     kend=-2
  else
     jbeg= j0
     jend= jmax
     kbeg= k0
     kend= kmax
  endif


  if (hy_useDiffuse) then
     !! Initialize
     magResist     = 0.0
     viscDynamic   = 0.0
     viscKinematic = 0.0
     cond          = 0.0
     dcff          = 0.0

     do k=kbeg-3,kend+3
        do j=jbeg-3,jend+3
           do i=i0-3,imax+3

              !! copy species to a temporal array
              speciesArr(:) = U(SPECIES_BEGIN:SPECIES_END,i,j,k)

              if (hy_useViscosity) then
                 !! Get viscosity
                 call Viscosity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k),&
                      speciesArr,viscDynamic(i,j,k),viscKinematic(i,j,k))
              endif

              if (hy_useConductivity) then
                 !! Get heat conductivity
                 call Conductivity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k),&
                      speciesArr,cond(i,j,k),dcff(i,j,k))
              endif

              if (hy_useMagneticResistivity) then
                 !! Get magnetic viscosity
                 call MagneticResistivity(U(TEMP_VAR,i,j,k),U(DENS_VAR,i,j,k),&
                      speciesArr,magResist(i,j,k))
              endif

           enddo
        enddo
     enddo

     !! For non-ideal fluxes
     magResist  = magResist/hy_mref
     viscDynamic= viscDynamic/hy_nref

     if (.not. sim_saturatedConduction) then
       !! For time-step calculation
       maxMagResist = maxval(magResist)
       maxKinVisc   = maxval(viscKinematic)
       maxDcff      = maxval(dcff)

       !! Maximum values for time-step calculation
       hy_maxMagDiff = max(maxMagResist,maxKinVisc,MaxDcff)
!print*,'in hy_uhd_getFaceFlux, hy_maxMagdiff=',hy_maxMagDiff
     endif

  endif





  !! Initialize flux arrays
  xflux = 0.
  yflux = 0.
  zflux = 0.

  !! Fix for shock-instabilities such as odd-even, carbuncle instabilities
  if (hy_shockInstabilityFix) then
     call hy_uhd_shockInstabilityFix(blockID,blkLimits,blkLimitsGC,LambdaFix)
  endif


  !! Compute intercell fluxes using the updated left & right states
  !! Calculate x-flux first
  do k=kbeg-3,kend+3
     do j=jbeg-3,jend+3
        do i=i0-1,imax+2
           VL(HY_DENS:HY_MAGZ)=scratch(i-1,j,k,XP01_SCRATCH_GRID_VAR:XP08_SCRATCH_GRID_VAR)
           VR(HY_DENS:HY_MAGZ)=scratch(i,  j,k,XN01_SCRATCH_GRID_VAR:XN08_SCRATCH_GRID_VAR)

           !! Left states
           eosData(EOS_DENS)=VL(HY_DENS)
           eosData(EOS_PRES)=VL(HY_PRES)
           eosData(EOS_TEMP)=0.5*(U(TEMP_VAR,i-1,j,k)+U(TEMP_VAR,i,j,k))
#ifdef YE_MSCALAR
           eosData(EOS_ABAR)=0.5*(1./U(SUMY_MSCALAR,i-1,j,k)&
                                 +1./U(SUMY_MSCALAR,i,  j,k))
           eosData(EOS_ZBAR)=0.5*(   U(YE_MSCALAR,  i-1,j,k)/U(SUMY_MSCALAR,i-1,j,k)&
                                    +U(YE_MSCALAR,  i,  j,k)/U(SUMY_MSCALAR,i,  j,k))
#endif
           call Eos(interp_eosMode,vecLen,eosData,eosMf,eosMask)
           VL(HY_GAMC) =eosData(EOS_GAMC)
           VL(HY_EINT) =eosData(EOS_EINT)
           VL(HY_GAME) =1.+eosData(EOS_PRES)/(eosData(EOS_DENS)*eosData(EOS_EINT)) ! game
           VL(HY_TEMP) =eosData(EOS_TEMP)


           !! Right states
           eosData(EOS_DENS)=VR(HY_DENS)
           eosData(EOS_PRES)=VR(HY_PRES)
           eosData(EOS_TEMP)=0.5*(U(TEMP_VAR,i,j,k)+U(TEMP_VAR,i+1,j,k))
#ifdef YE_MSCALAR
           eosData(EOS_ABAR)=0.5*(1./U(SUMY_MSCALAR,i,  j,k)&
                                 +1./U(SUMY_MSCALAR,i+1,j,k))
           eosData(EOS_ZBAR)=0.5*(   U(YE_MSCALAR,  i,  j,k)/U(SUMY_MSCALAR,i,  j,k)&
                                    +U(YE_MSCALAR,  i+1,j,k)/U(SUMY_MSCALAR,i+1,j,k))
#endif
           call Eos(interp_eosMode,vecLen,eosData,eosMf,eosMask)
           VR(HY_GAMC) =eosData(EOS_GAMC)
           VR(HY_EINT) =eosData(EOS_EINT)
           VR(HY_GAME) =1.+eosData(EOS_PRES)/(eosData(EOS_DENS)*eosData(EOS_EINT)) ! game
           VR(HY_TEMP) =eosData(EOS_TEMP)

           if (hy_RiemannSolver == ROE) then

              call hy_uhd_Roe(DIR_X,VL,VR,xflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_X,i,j,k),&
                              i,j,k,blkLimitsGC,LambdaFix)

           elseif (hy_RiemannSolver == HLL) then
              call hy_uhd_HLL(DIR_X,VL,VR,xflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_X,i,j,k))

           elseif (hy_RiemannSolver == HLLC) then
              call hy_uhd_HLLC(DIR_X,VL,VR,xflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_X,i,j,k),&
                               i,j,k,blkLimitsGC,LambdaFix)

           elseif (hy_RiemannSolver == HLLD) then
              call hy_uhd_HLLD(DIR_X,VL,VR,xflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_X,i,j,k),&
                               i,j,k,blkLimitsGC,LambdaFix)

           elseif (hy_RiemannSolver == HLL7) then
              call hy_uhd_HLL7(DIR_X,VL,VR,xflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_X,i,j,k),&
                               i,j,k,blkLimitsGC,LambdaFix)

           elseif (hy_RiemannSolver == LF) then
              call hy_uhd_LF(DIR_X,VL,VR,xflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_X,i,j,k))

           elseif (hy_RiemannSolver == HYBR) then
              if (scratch(i,j,k,FLAG_SCRATCH_GRID_VAR)==0.) then
                 call hy_uhd_Roe(DIR_X,VL,VR,xflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_X,i,j,k),&
                                 i,j,k,blkLimitsGC,LambdaFix)
              else
                 call hy_uhd_HLLD(DIR_X,VL,VR,xflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_X,i,j,k),&
                                  i,j,k,blkLimitsGC,LambdaFix)
              endif
           endif


           !! Flux for internal energy density
           !! Reference: "Simple Method to Track Pressure Accurately", S. Li, Astronum Proceeding, 2007
           if (xflux(F01DENS_FLUX,i,j,k) <= 0.) then
              xflux(F09EINT_FLUX,i,j,k) = xflux(F01DENS_FLUX,i,j,k)*VR(HY_EINT)
              xflux(F10PRES_FLUX,i,j,k) = xflux(F01DENS_FLUX,i,j,k)/VR(HY_DENS)
           else
              xflux(F09EINT_FLUX,i,j,k) = xflux(F01DENS_FLUX,i,j,k)*VL(HY_EINT)
              xflux(F10PRES_FLUX,i,j,k) = xflux(F01DENS_FLUX,i,j,k)/VL(HY_DENS)
           endif


!#ifdef MAGNETIC_RESISTIVITY
           if (hy_useDiffuse) then
              if (hy_useViscosity) then
                 call hy_uhd_addViscousFluxes(blockID,blkLimitsGC,i,j,k,&
                      xflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),viscDynamic,DIR_X)
              endif
              if (hy_useConductivity) then
                 call hy_uhd_addThermalFluxes(blockID,blkLimitsGC,i,j,k,&
                      xflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),cond,DIR_X,Dcffx(i,j,k))
              endif
              if (hy_useMagneticResistivity) then
                 call hy_uhd_addResistiveFluxes(blockID,blkLimitsGC,i,j,k,&
                      xflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),magResist,DIR_X)
              endif
           endif
!#endif

#if NDIM >= 2
#ifndef FLASH_GRID_UG
           !! Store this only for multidimensional setups
           scratch(i,j,k,FLX01_SCRATCH_GRID_VAR:FLX10_SCRATCH_GRID_VAR) = &
                xflux(F01DENS_FLUX:F10PRES_FLUX,i,j,k)
#endif
#endif
        enddo
     enddo
  enddo

#if NDIM >= 2
  !! Calculate y-flux
  do k=kbeg-3,kend+3
     do j=jbeg-1,jend+2
        do i=i0-3,imax+3
           VL(HY_DENS:HY_MAGZ)=scratch(i,j-1,k,YP01_SCRATCH_GRID_VAR:YP08_SCRATCH_GRID_VAR)
           VR(HY_DENS:HY_MAGZ)=scratch(i,j,  k,YN01_SCRATCH_GRID_VAR:YN08_SCRATCH_GRID_VAR)

           !! Left states
           eosData(EOS_DENS)=VL(HY_DENS)
           eosData(EOS_PRES)=VL(HY_PRES)
           eosData(EOS_TEMP)=0.5*(U(TEMP_VAR,i,j-1,k)+U(TEMP_VAR,i,j,k))
#ifdef YE_MSCALAR
           eosData(EOS_ABAR)=0.5*(1./U(SUMY_MSCALAR,i,j-1,k)&
                                 +1./U(SUMY_MSCALAR,i,j,  k))
           eosData(EOS_ZBAR)=0.5*(   U(YE_MSCALAR,  i,j-1,k)/U(SUMY_MSCALAR,i,j-1,k)&
                                    +U(YE_MSCALAR,  i,j,  k)/U(SUMY_MSCALAR,i,j,  k))
#endif
           call Eos(interp_eosMode,vecLen,eosData,eosMf,eosMask)
           VL(HY_GAMC) =eosData(EOS_GAMC)
           VL(HY_EINT) =eosData(EOS_EINT)
           VL(HY_GAME) =1.+eosData(EOS_PRES)/(eosData(EOS_DENS)*eosData(EOS_EINT)) ! game
           VL(HY_TEMP) =eosData(EOS_TEMP)


           !! Right states
           eosData(EOS_DENS)=VR(HY_DENS)
           eosData(EOS_PRES)=VR(HY_PRES)
           eosData(EOS_TEMP)=0.5*(U(TEMP_VAR,i,j,k)+U(TEMP_VAR,i,j+1,k))
#ifdef YE_MSCALAR
           eosData(EOS_ABAR)=0.5*(1./U(SUMY_MSCALAR,i,j,  k)&
                                 +1./U(SUMY_MSCALAR,i,j+1,k))
           eosData(EOS_ZBAR)=0.5*(   U(YE_MSCALAR,  i,j,  k)/U(SUMY_MSCALAR,i,j,  k)&
                                    +U(YE_MSCALAR,  i,j+1,k)/U(SUMY_MSCALAR,i,j+1,k))
#endif
           call Eos(interp_eosMode,vecLen,eosData,eosMf,eosMask)
           VR(HY_GAMC) =eosData(EOS_GAMC)
           VR(HY_EINT) =eosData(EOS_EINT)
           VR(HY_GAME) =1.+eosData(EOS_PRES)/(eosData(EOS_DENS)*eosData(EOS_EINT)) ! game
           VR(HY_TEMP) =eosData(EOS_TEMP)


           if (hy_RiemannSolver == ROE) then

              call hy_uhd_Roe(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Y,i,j,k),&
                              i,j,k,blkLimitsGC,LambdaFix)

           elseif (hy_RiemannSolver == HLL) then
              call hy_uhd_HLL(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Y,i,j,k))

           elseif (hy_RiemannSolver == HLLC) then
              call hy_uhd_HLLC(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Y,i,j,k),&
                               i,j,k,blkLimitsGC,LambdaFix)

           elseif (hy_RiemannSolver == HLLD) then
              call hy_uhd_HLLD(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Y,i,j,k),&
                               i,j,k,blkLimitsGC,LambdaFix)

           elseif (hy_RiemannSolver == HLL7) then
              call hy_uhd_HLL7(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Y,i,j,k),&
                               i,j,k,blkLimitsGC,LambdaFix)

           elseif (hy_RiemannSolver == LF) then
              call hy_uhd_LF(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Y,i,j,k))

           elseif (hy_RiemannSolver == HYBR) then
              if (scratch(i,j,k,FLAG_SCRATCH_GRID_VAR)==0.) then
                 call hy_uhd_Roe(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Y,i,j,k),&
                                 i,j,k,blkLimitsGC,LambdaFix)
              else
                 call hy_uhd_HLLD(DIR_Y,VL,VR,yflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Y,i,j,k),&
                                  i,j,k,blkLimitsGC,LambdaFix)
              endif
           endif


           !! Flux for internal energy density
           !! Reference: "Simple Method to Track Pressure Accurately", S. Li, Astronum Proceeding, 2007
           if (yflux(F01DENS_FLUX,i,j,k) <= 0.) then
              yflux(F09EINT_FLUX,i,j,k) = yflux(F01DENS_FLUX,i,j,k)*VR(HY_EINT)
              yflux(F10PRES_FLUX,i,j,k) = yflux(F01DENS_FLUX,i,j,k)/VR(HY_DENS)
           else
              yflux(F09EINT_FLUX,i,j,k) = yflux(F01DENS_FLUX,i,j,k)*VL(HY_EINT)
              yflux(F10PRES_FLUX,i,j,k) = yflux(F01DENS_FLUX,i,j,k)/VL(HY_DENS)
           endif


!#ifdef MAGNETIC_RESISTIVITY
           if (hy_useDiffuse) then
              if (hy_useViscosity) then
                 call hy_uhd_addViscousFluxes(blockID,blkLimitsGC,i,j,k,&
                      yflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),viscDynamic,DIR_Y)
              endif
              if (hy_useConductivity) then
                 call hy_uhd_addThermalFluxes(blockID,blkLimitsGC,i,j,k,&
                      yflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),cond,DIR_Y,Dcffy(i,j,k))
              endif
              if (hy_useMagneticResistivity) then
                 call hy_uhd_addResistiveFluxes(blockID,blkLimitsGC,i,j,k,&
                      yflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),magResist,DIR_Y)
              endif
           endif
!#endif
#ifndef FLASH_GRID_UG
           scratch(i,j,k,FLY01_SCRATCH_GRID_VAR:FLY10_SCRATCH_GRID_VAR) = &
                yflux(F01DENS_FLUX:F10PRES_FLUX,i,j,k)
#endif

        enddo
     enddo
  enddo

#if NDIM == 3    
  !! Calculate z-flux
  do k=kbeg-1,kend+2
     do j=jbeg-3,jend+3
        do i=i0-3,imax+3
           VL(HY_DENS:HY_MAGZ)=scratch(i,j,k-1,ZP01_SCRATCH_GRID_VAR:ZP08_SCRATCH_GRID_VAR)
           VR(HY_DENS:HY_MAGZ)=scratch(i,j,k  ,ZN01_SCRATCH_GRID_VAR:ZN08_SCRATCH_GRID_VAR)

           !! Left states
           eosData(EOS_DENS)=VL(HY_DENS)
           eosData(EOS_PRES)=VL(HY_PRES)
           eosData(EOS_TEMP)=0.5*(U(TEMP_VAR,i,j,k-1)+U(TEMP_VAR,i,j,k))
#ifdef YE_MSCALAR
           eosData(EOS_ABAR)=0.5*(1./U(SUMY_MSCALAR,i,j,k-1)&
                                 +1./U(SUMY_MSCALAR,i,j,k  ))
           eosData(EOS_ZBAR)=0.5*(   U(YE_MSCALAR,  i,j,k-1)/U(SUMY_MSCALAR,i,j,k-1)&
                                    +U(YE_MSCALAR,  i,j,k  )/U(SUMY_MSCALAR,i,j,k  ))

#endif
           call Eos(interp_eosMode,vecLen,eosData,eosMf,eosMask)
           VL(HY_GAMC) =eosData(EOS_GAMC)
           VL(HY_EINT) =eosData(EOS_EINT)
           VL(HY_GAME) =1.+eosData(EOS_PRES)/(eosData(EOS_DENS)*eosData(EOS_EINT)) ! game
           VL(HY_TEMP) =eosData(EOS_TEMP)


           !! Right states
           eosData(EOS_DENS)=VR(HY_DENS)
           eosData(EOS_PRES)=VR(HY_PRES)
           eosData(EOS_TEMP)=0.5*(U(TEMP_VAR,i,j,k)+U(TEMP_VAR,i,j,k+1))
#ifdef YE_MSCALAR
           eosData(EOS_ABAR)=0.5*(1./U(SUMY_MSCALAR,i,j,k  )&
                                 +1./U(SUMY_MSCALAR,i,j,k+1))
           eosData(EOS_ZBAR)=0.5*(   U(YE_MSCALAR,  i,j,k  )/U(SUMY_MSCALAR,i,j,k  )&
                                    +U(YE_MSCALAR,  i,j,k+1)/U(SUMY_MSCALAR,i,j,k+1))
#endif
           call Eos(interp_eosMode,vecLen,eosData,eosMf,eosMask)
           VR(HY_GAMC) =eosData(EOS_GAMC)
           VR(HY_EINT) =eosData(EOS_EINT)
           VR(HY_GAME) =1.+eosData(EOS_PRES)/(eosData(EOS_DENS)*eosData(EOS_EINT)) ! game
           VR(HY_TEMP) =eosData(EOS_TEMP)


           if (hy_RiemannSolver == ROE) then
              call hy_uhd_Roe(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Z,i,j,k),&
                              i,j,k,blkLimitsGC,LambdaFix)

           elseif (hy_RiemannSolver == HLL) then
              call hy_uhd_HLL(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Z,i,j,k))

           elseif (hy_RiemannSolver == HLLC) then
              call hy_uhd_HLLC(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Z,i,j,k),&
                               i,j,k,blkLimitsGC,LambdaFix)

           elseif (hy_RiemannSolver == HLLD) then
              call hy_uhd_HLLD(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Z,i,j,k),&
                               i,j,k,blkLimitsGC,LambdaFix)

           elseif (hy_RiemannSolver == HLL7) then
              call hy_uhd_HLL7(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Z,i,j,k),&
                               i,j,k,blkLimitsGC,LambdaFix)

           elseif (hy_RiemannSolver == LF) then
              call hy_uhd_LF(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Z,i,j,k))

           elseif (hy_RiemannSolver == HYBR) then
              if (scratch(i,j,k,FLAG_SCRATCH_GRID_VAR)==0.) then
                 call hy_uhd_Roe(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Z,i,j,k),&
                                 i,j,k,blkLimitsGC,LambdaFix)
              else
                 call hy_uhd_HLLD(DIR_Z,VL,VR,zflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),vint(DIR_Z,i,j,k),&
                                  i,j,k,blkLimitsGC,LambdaFix)
              endif
           endif


           !! Flux for internal energy density
           !! Reference: "Simple Method to Track Pressure Accurately", S. Li, Astronum Proceeding, 2007
           if (zflux(F01DENS_FLUX,i,j,k) <= 0.) then
              zflux(F09EINT_FLUX,i,j,k) = zflux(F01DENS_FLUX,i,j,k)*VR(HY_EINT)
              zflux(F10PRES_FLUX,i,j,k) = zflux(F01DENS_FLUX,i,j,k)/VR(HY_DENS)
           else
              zflux(F09EINT_FLUX,i,j,k) = zflux(F01DENS_FLUX,i,j,k)*VL(HY_EINT)
              zflux(F10PRES_FLUX,i,j,k) = zflux(F01DENS_FLUX,i,j,k)/VL(HY_DENS)
           endif


!#ifdef MAGNETIC_RESISTIVITY
           if (hy_useDiffuse) then
              if (hy_useViscosity) then
                 call hy_uhd_addViscousFluxes(blockID,blkLimitsGC,i,j,k,&
                      zflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),viscDynamic,DIR_Z)
              endif
              if (hy_useConductivity) then
                 call hy_uhd_addThermalFluxes(blockID,blkLimitsGC,i,j,k,&
                      zflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),cond,DIR_Z,Dcffz(i,j,k))
              endif
              if (hy_useMagneticResistivity) then
                 call hy_uhd_addResistiveFluxes(blockID,blkLimitsGC,i,j,k,&
                      zflux(F01DENS_FLUX:F08MAGZ_FLUX,i,j,k),magResist,DIR_Z)
              endif
           endif
!#endif
#ifndef FLASH_GRID_UG
           scratch(i,j,k,FLZ01_SCRATCH_GRID_VAR:FLZ10_SCRATCH_GRID_VAR) = &
                zflux(F01DENS_FLUX:F10PRES_FLUX,i,j,k)
#endif
        enddo
     enddo
  enddo
#endif

#endif

  !! Release pointer
  call Grid_releaseBlkPtr(blockID,scratch,SCRATCH)
  call Grid_releaseBlkPtr(blockID,U,CENTER)

  if (sim_saturatedConduction) then
    ! SB modification
    maxDcffx   = maxval(Dcffx)
    maxDcffy   = maxval(Dcffy)
    maxDcffz   = maxval(Dcffz) 
    maxDcff    = max(maxDcffx, maxDcffy, maxDcffz)
    hy_maxMagDiff = max(maxMagResist,maxKinVisc,MaxDcff)
  endif

  ! ------------------

End Subroutine hy_uhd_getFaceFlux
