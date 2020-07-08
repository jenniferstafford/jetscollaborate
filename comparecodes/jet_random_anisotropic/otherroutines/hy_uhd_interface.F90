!!****ih* source/physics/Hydro/HydroMain/unsplit/hy_uhd_interface
!!
!! NAME
!!  hy_uhd_interface
!!
!! SYNOPSIS
!!  use hy_uhd_interface
!!
!! DESCRIPTION
!!  This is an interface specific for the two unsplit solvers:
!!   1. MHD_StaggeredMesh
!!   2. Hydro_MusclHancock
!!  This interface defines its public interfaces.
!!
!!***

Module hy_uhd_interface

  implicit none

#include "constants.h"
#include "Flash.h"
#include "UHD.h"


  interface
     subroutine hy_uhd_avgState(VL,VR,Vavg)
       implicit none
       real, dimension(HY_VARINUM+3), intent(IN)  :: VL,VR
       real, dimension(HY_VARINUM+2), intent(OUT) :: Vavg
     end subroutine hy_uhd_avgState
  end interface



  interface
     subroutine hy_uhd_getRiemannState(blockID,blkLimits,blkLimitsGC,hy_limiter,&
                                       dt,del,gravX,gravY,gravZ,normalFieldUpdate)
       implicit none
       integer, intent(IN)   :: blockID,hy_limiter
       real,    intent(IN)   :: dt
       real,    intent(IN),dimension(MDIM) :: del
       integer, intent(IN),dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
#ifdef FIXEDBLOCKSIZE
       real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(IN) &
            :: gravX,gravY,gravZ
#else
       real, dimension(blkLimitsGC(HIGH,IAXIS),  &
                       blkLimitsGC(HIGH,JAXIS),  &
                       blkLimitsGC(HIGH,KAXIS)), &
            intent(IN) :: gravX,gravY,gravZ
#endif
       logical, intent(IN), optional :: normalFieldUpdate
     end subroutine hy_uhd_getRiemannState
 end interface



  interface
     subroutine hy_uhd_entropyFix(lambda,lambdaL,lambdaR)
       implicit none
       real, dimension(HY_WAVENUM), intent(INOUT) :: lambda
       real, dimension(HY_WAVENUM), intent(IN)    :: lambdaL,lambdaR
     end subroutine hy_uhd_entropyFix
  end interface



  interface
     subroutine hy_uhd_getFaceFlux ( blockID,blkLimits,blkLimitsGC, &
                                     datasize, xflux, yflux,zflux,vint)
       implicit none
       integer, intent(IN)  :: blockID
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
       integer, dimension(MDIM),     intent(IN) :: datasize
#ifdef FIXEDBLOCKSIZE
       real, DIMENSION(NFLUXES,              &
                       GRID_ILO_GC:GRID_IHI_GC, &
                       GRID_JLO_GC:GRID_JHI_GC, &
                       GRID_KLO_GC:GRID_KHI_GC),&
                       intent(OUT) :: xflux, yflux, zflux
       real, DIMENSION(MDIM,              &
                       GRID_ILO_GC:GRID_IHI_GC, &
                       GRID_JLO_GC:GRID_JHI_GC, &
                       GRID_KLO_GC:GRID_KHI_GC),&
                       intent(OUT) :: vint
#else
       real, DIMENSION(NFLUXES, &
                       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
                       intent(OUT) :: xflux, yflux, zflux
       real, DIMENSION(MDIM,&
                       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
                       intent(OUT) :: vint
#endif
     end subroutine hy_uhd_getFaceFlux
  end interface



  interface
     subroutine hy_uhd_dataReconstOneStep&
          (blockID,dataSize,hy_limiter,ix,iy,iz,dt,del,gravX,gravY,gravZ,&
          V0,Vxp,Vxn,Vyp,Vyn,Vzp,Vzn,Wxp,Wxn,Wyp,Wyn,Wzp,Wzn)
       implicit none
       integer,intent(IN) :: blockID,hy_limiter
       integer,intent(IN), dimension(MDIM) :: dataSize
       integer,intent(IN) :: ix,iy,iz
       real,   intent(IN) :: dt
       real,   intent(IN), dimension(MDIM) :: del
       real,   intent(IN) :: gravX,gravY,gravZ
       real,   intent(IN),  dimension(HY_VARINUM+2) :: V0,  Vxp, Vxn, Vyp, Vyn, Vzp, Vzn
       real,   intent(OUT), dimension(HY_VARINUM)   :: Wxp, Wxn, Wyp, Wyn, Wzp, Wzn
     end subroutine hy_uhd_dataReconstOnestep
  end interface



  interface
     subroutine hy_uhd_Roe(dir,Vm,Vp,Fstar,vint,ix,iy,iz,blkLimitsGC,LambdaFix)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUM+4), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM), intent(OUT):: Fstar
       real, intent(OUT) :: vint
       integer, intent(IN) :: ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
#ifdef FIXEDBLOCKSIZE
       real, dimension(MDIM,2,&
            GRID_ILO_GC:GRID_IHI_GC, &
            GRID_JLO_GC:GRID_JHI_GC, &
            GRID_KLO_GC:GRID_KHI_GC),&
            intent(IN) :: LambdaFix
#else
       real, dimension(MDIM,2,&
            blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
            intent(IN) :: LambdaFix
#endif
     end subroutine hy_uhd_Roe
  end interface



  interface
     subroutine hy_uhd_HLL(dir,Vm,Vp,Fstar,vint)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUM+4), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM),  intent(OUT) :: Fstar
       real, intent(OUT) :: vint
     end subroutine hy_uhd_HLL
  end interface



  interface
     subroutine hy_uhd_HLLC(dir,Vm,Vp,Fstar,vint,ix,iy,iz,blkLimitsGC,LambdaFix)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUM+4), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM), intent(OUT):: Fstar
       real, intent(OUT) :: vint
       integer, intent(IN) :: ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
#ifdef FIXEDBLOCKSIZE
       real, dimension(MDIM,2,&
            GRID_ILO_GC:GRID_IHI_GC, &
            GRID_JLO_GC:GRID_JHI_GC, &
            GRID_KLO_GC:GRID_KHI_GC),&
            intent(IN) :: LambdaFix
#else
       real, dimension(MDIM,2,&
            blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
            intent(IN) :: LambdaFix
#endif
     end subroutine hy_uhd_HLLC
  end interface



  interface
     subroutine hy_uhd_LF(dir,Vm,Vp,Fstar,vint)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUM+4), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM), intent(OUT):: Fstar
       real, intent(OUT) :: vint
     end subroutine hy_uhd_LF
  end interface



  interface
     subroutine hy_uhd_TVDslope&
          (dataSize,limiter_opt,dir,V0,VR,VL,leig,delbar)
       implicit none
       integer, intent(IN) :: limiter_opt, dir
       integer,dimension(MDIM),intent(IN) :: dataSize
       real,dimension(HY_VARINUM+2),intent(IN) :: V0,VR,VL
       real,dimension(HY_WAVENUM,HY_VARINUM),intent(IN) :: leig
       real,dimension(HY_VARINUM),intent(OUT) :: delbar
     end subroutine hy_uhd_TVDslope
  end interface


  interface
     subroutine hy_uhd_unsplit( myPE, blockCount, blockList, dt, dtOld )
       implicit none
       integer, INTENT(IN) :: myPE
       integer, INTENT(IN) :: blockCount  
       integer, INTENT(IN), dimension(blockCount) :: blockList
       real,    INTENT(IN) :: dt, dtOld
     end subroutine hy_uhd_unsplit
  end interface



  interface
     subroutine hy_uhd_unsplitUpdate(blockID,range_switch,dt,del,dataSize,blkLimits,&    
                                  xflux,yflux,zflux,gravX,gravY,gravZ)
       implicit none
       integer,intent(IN) :: blockID
       integer,intent(IN) :: range_switch
       real, intent(IN)   :: dt
       real, intent(IN)   :: del(MDIM)
       integer,dimension(MDIM),intent(IN) :: dataSize
       integer,intent(IN) :: blkLimits(LOW:HIGH,MDIM)
#ifdef FIXEDBLOCKSIZE
       real, dimension(NFLUXES,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), &
            INTENT(IN) :: xflux, yflux, zflux
       real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), &
            intent(IN) :: gravX,gravY,gravZ

#else
       real, dimension(NFLUXES,dataSize(IAXIS),dataSize(JAXIS),&
                               dataSize(KAXIS)), INTENT(IN) :: xflux, yflux, zflux
  real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)), intent(IN) &
       :: gravX,gravY,gravZ

#endif
     end subroutine hy_uhd_unsplitUpdate
  end interface



  interface
     subroutine hy_uhd_updateSpecies(blockID,blkLimits,blkLimitsGC,del,vint)
       implicit none
       integer, intent(IN)  :: blockID
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits,blkLimitsGC
       real, dimension(MDIM), intent(IN) :: del
#ifdef FIXEDBLOCKSIZE
       real, dimension(MDIM,GRID_ILO_GC:GRID_IHI_GC, &
                            GRID_JLO_GC:GRID_JHI_GC, &
                            GRID_KLO_GC:GRID_KHI_GC),&
                            intent(IN) :: vint
#else
       real, dimension(MDIM,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
                            intent(IN) :: vint
#endif
     end subroutine hy_uhd_updateSpecies
  end interface



  interface
     subroutine hy_uhd_energyFix(blockID,blkLimits,dt,del,eosMode)
       implicit none
       integer, intent(IN) :: blockID
       integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
       real, intent(IN) :: dt
       real, dimension(MDIM), intent(IN) :: del
       integer, intent(IN) :: eosMode
     end subroutine hy_uhd_energyFix
  end interface



  interface 
     subroutine hy_uhd_unitConvert(blockID,convertDir)
       implicit none
       integer, intent(IN) :: blockID, convertDir
     end subroutine hy_uhd_unitConvert
  end interface
  


  interface
     subroutine hy_uhd_addViscousFluxes&
          (blockID,blkLimitsGC,ix,iy,iz,Flux,visc,sweepDir)
       implicit none
       integer, INTENT(IN) :: blockID,ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
#ifdef FIXEDBLOCKSIZE
       real, dimension(GRID_ILO_GC:GRID_IHI_GC, &
                       GRID_JLO_GC:GRID_JHI_GC, &
                       GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: visc
#else
       real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
                       intent(IN) :: visc
#endif
       integer, INTENT(IN) :: sweepDir
     end subroutine hy_uhd_addViscousFluxes
    end interface



  interface
     subroutine hy_uhd_addThermalFluxes&
          (blockID,blkLimitsGC,ix,iy,iz,Flux,cond,sweepDir,Dcff)
       implicit none
       integer, INTENT(IN) :: blockID,ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
#ifdef FIXEDBLOCKSIZE
       real, dimension(GRID_ILO_GC:GRID_IHI_GC, &
                       GRID_JLO_GC:GRID_JHI_GC, &
                       GRID_KLO_GC:GRID_KHI_GC), intent(IN) :: cond
#else
       real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
                       intent(IN) :: cond
#endif
       integer, INTENT(IN) :: sweepDir
       real, INTENT(OUT) :: Dcff
     end subroutine hy_uhd_addThermalFluxes
    end interface



    interface
       subroutine hy_uhd_eigenParameters(V,dir,cons,U_normal,C_fast,C_alfn,C_slow,A_f,A_s,B_beta)
         implicit none
         real, dimension(HY_VARINUM+2), intent(IN)  :: V
         integer, intent(IN)  :: dir
         logical, intent(IN)  :: cons
         real, intent(OUT) :: U_normal,C_fast
         real, intent(OUT), optional :: C_alfn,C_slow,A_f,A_s
         real, dimension(MDIM), intent(OUT),optional :: B_beta
       end subroutine hy_uhd_eigenParameters
    end interface



    interface
       subroutine hy_uhd_eigenValue(EigValue,U_normal,C_fast,C_alfn,C_slow)
         implicit none
         real,dimension(HY_WAVENUM), intent(OUT) :: EigValue
         real,intent(IN) :: U_normal,C_fast
         real,intent(IN), optional :: C_alfn,C_slow
       end subroutine hy_uhd_eigenValue
    end interface



    interface
       subroutine hy_uhd_eigenVector&
            (LeftEigvec,RightEigvec,V,dir,cons,C_fast,C_alfn,C_slow,A_f,A_s,B_beta)
         implicit none
         real, dimension(HY_WAVENUM,HY_VARINUM), intent(OUT) :: LeftEigvec
         real, dimension(HY_VARINUM,HY_WAVENUM), intent(OUT) :: RightEigvec
         real, dimension(HY_VARINUM+2), intent(IN) :: V
         integer, intent(IN) :: dir
         logical, intent(IN) :: cons
         real, intent(IN) :: C_fast
         real, intent(IN),optional :: C_alfn,C_slow,A_f,A_s
         real, dimension(MDIM), intent(IN),optional  :: B_beta
       end subroutine hy_uhd_eigenVector
    end interface



    interface
       subroutine hy_uhd_putGravityUnsplit&
	(blockID,blkLimitsGC,dataSize,dt,dtOld,gravX,gravY,gravZ)
         implicit none
         integer, intent(IN) :: blockID
         integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimitsGC
         integer, dimension(MDIM), intent(IN) :: dataSize
	 real,    intent(IN) :: dt, dtOld
#ifdef FIXEDBLOCKSIZE
         real, dimension(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC), intent(INOUT) :: &
              gravX,gravY,gravZ
#else
         real, dimension(dataSize(IAXIS),dataSize(JAXIS),dataSize(KAXIS)), intent(INOUT) :: &
              gravX,gravY,gravZ
#endif
       end subroutine hy_uhd_putGravityUnsplit
    end interface



    interface
       subroutine hy_uhd_shockDetect(blockID)
         implicit none
         integer, INTENT(IN) :: blockID
       end subroutine hy_uhd_shockDetect
    end interface



    interface
       subroutine hy_uhd_prim2con(V,CU)
         implicit none
         real ,dimension(HY_VARINUM+2), intent(IN) :: V
         real ,dimension(HY_VARINUM),  intent(OUT) :: CU
       end subroutine hy_uhd_prim2con
    end interface



    interface
       subroutine hy_uhd_con2prim(CU,game,V)
         implicit none
         real ,dimension(HY_VARINUM), intent(IN)  :: CU
         real, intent(IN) :: game
         real ,dimension(HY_VARINUM), intent(OUT) :: V
       end subroutine hy_uhd_con2prim
    end interface



    interface
       subroutine hy_uhd_prim2flx(dir,V,F)
         implicit none
         integer, intent(IN)  :: dir
         real, dimension(HY_VARINUM+2), intent(IN) :: V
         real, dimension(HY_VARINUM),  intent(OUT) :: F
       end subroutine hy_uhd_prim2flx
    end interface



    interface
       subroutine hy_uhd_shockInstabilityFix(blockID,blkLimits,blkLimitsGC,LambdaFix)
         implicit none
         integer, intent(IN)  :: blockID
         integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC

#ifdef FIXEDBLOCKSIZE
         real, dimension(MDIM,2,&
              GRID_ILO_GC:GRID_IHI_GC, &
              GRID_JLO_GC:GRID_JHI_GC, &
              GRID_KLO_GC:GRID_KHI_GC),&
              intent(OUT):: LambdaFix
#else
         real, dimension(MDIM,2,&
              blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
              blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
              blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
              intent(OUT) :: LambdaFix
#endif
       end subroutine hy_uhd_shockInstabilityFix
    end interface



!! FOR UNSPLIT STAGGERED MESH MHD SOLVER -------------------------------------------
#ifdef FLASH_USM_MHD
  interface
     subroutine hy_uhd_addResistiveFluxes&
          (blockID,blkLimitsGC,ix,iy,iz,Flux,magVisc,sweepDir)
       implicit none
       integer, INTENT(IN) :: blockID,ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
#ifdef FIXEDBLOCKSIZE
       real, dimension(GRID_ILO_GC:GRID_IHI_GC, &
                       GRID_JLO_GC:GRID_JHI_GC, &
                       GRID_KLO_GC:GRID_KHI_GC),intent(IN) :: magVisc
#else
       real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
                       intent(IN) :: magVisc
#endif
       integer, INTENT(IN) :: sweepDir
     end subroutine hy_uhd_addResistiveFluxes
    end interface



  interface
     subroutine hy_uhd_staggeredDivb(blockID,dt,del,blkLimits,range_switch,halfTimeAdvance)
       implicit none
       integer, intent(IN) :: blockID
       real,    intent(IN) :: dt
       real,    dimension(MDIM),   intent(IN) :: del
       integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
       integer, intent(IN) :: range_switch
       logical, intent(IN) :: halfTimeAdvance
     end subroutine hy_uhd_staggeredDivb
  end interface



  interface
     subroutine hy_uhd_HLLD(dir,Vm,Vp,Fstar,vint,ix,iy,iz,blkLimitsGC,LambdaFix)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUM+4), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM), intent(OUT):: Fstar
       real, intent(OUT) :: vint
       integer, intent(IN) :: ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
#ifdef FIXEDBLOCKSIZE
       real, dimension(MDIM,2,&
            GRID_ILO_GC:GRID_IHI_GC, &
            GRID_JLO_GC:GRID_JHI_GC, &
            GRID_KLO_GC:GRID_KHI_GC),&
            intent(IN) :: LambdaFix
#else
       real, dimension(MDIM,2,&
            blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
            intent(IN) :: LambdaFix
#endif
     end subroutine hy_uhd_HLLD
  end interface


  interface
     subroutine hy_uhd_HLL7(dir,Vm,Vp,Fstar,vint,ix,iy,iz,blkLimitsGC,LambdaFix)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUM+4), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM), intent(OUT):: Fstar
       real, intent(OUT) :: vint
       integer, intent(IN) :: ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
#ifdef FIXEDBLOCKSIZE
       real, dimension(MDIM,2,&
            GRID_ILO_GC:GRID_IHI_GC, &
            GRID_JLO_GC:GRID_JHI_GC, &
            GRID_KLO_GC:GRID_KHI_GC),&
            intent(IN) :: LambdaFix
#else
       real, dimension(MDIM,2,&
            blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
            blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
            blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
            intent(IN) :: LambdaFix
#endif
     end subroutine hy_uhd_HLL7
  end interface


  interface
     subroutine hy_uhd_getElectricFields( blockID,blkLimits,blkLimitsGC,del,flx,fly,flz)
       implicit none
       integer, intent(IN)  :: blockID
       integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
       real,    intent(IN), dimension(MDIM)  :: del  
#ifdef FIXEDBLOCKSIZE
       real, DIMENSION(NFLUXES,               &
                       GRID_ILO_GC:GRID_IHI_GC,  &
                       GRID_JLO_GC:GRID_JHI_GC,  &
                       GRID_KLO_GC:GRID_KHI_GC), &
                       intent(IN) :: flx,fly,flz
#else
       real, DIMENSION(NFLUXES,               &
                       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), &
                       intent(IN) :: flx,fly,flz
#endif
     end subroutine hy_uhd_getElectricFields
  end interface


  interface
     subroutine hy_uhd_getFluxDeriv( ix,iy,iz,blkLimitsGC,&
                                     fluxType,DerivDir,   &
                                     faceFlux,            &
                                     Flux1Deriv,          &
                                     Flux2Deriv           )
       implicit none
       integer, intent(IN) :: ix,iy,iz
       integer, intent(IN) :: fluxType,DerivDir
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
#ifdef FIXEDBLOCKSIZE
       real,    dimension(NFLUXES,                  &
                          GRID_ILO_GC:GRID_IHI_GC,  &
                          GRID_JLO_GC:GRID_JHI_GC,  &
                          GRID_KLO_GC:GRID_KHI_GC), &
                          intent(IN) :: faceFlux
#else
       real,    dimension(NFLUXES,                  &
                          blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                          blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                          blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), &
                          intent(IN) :: faceFlux
#endif
       real,    intent(OUT):: Flux1Deriv,Flux2Deriv
     end subroutine hy_uhd_getFluxDeriv
  end interface

#endif !endif #ifdef FLASH_UNSPLIT_MHD

End Module hy_uhd_interface
