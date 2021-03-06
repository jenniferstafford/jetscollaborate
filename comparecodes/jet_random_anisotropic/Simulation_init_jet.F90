!!****if* source/Simulation/SimulationMain/magnetoHD/MHD_Jet/Simulation_init.F90
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the wind tunnel with a step problem
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use Simulation_jetNozzleUpdate, ONLY : sim_jetNozzleUpdate
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_data, ONLY : dr_simTime, dr_dt, dr_restart
  use Driver_interface, ONLY : Driver_getMype
  use Grid_data, ONLY : gr_minCellSize
  use IO_interface, ONLY :  IO_getScalar
  use Particles_data, ONLY : pt_randSeed
  use IO_data, ONLY: io_nextCheckpointTime, io_checkpointFileIntervalTime,\
                     io_nextPlotFileTime, io_plotFileIntervalTime
  use IOParticles_data, ONLY: io_nextParticleFileTime, io_particleFileIntervalTime

  !use Grid_data, ONLY : gr_smallrho

  implicit none
#include "constants.h"
#include "Flash.h"

  integer :: nozzle=1, clock, iexist
  real    :: maxPrecession
  logical :: resetOutputTime

  call Driver_getMype(MESH_COMM, sim_meshMe)

  ! Initialize the random number generator
  call system_clock(count=clock)
  pt_randSeed = (/clock, sim_meshMe/)

  call RuntimeParameters_get('smlrho', sim_smlrho)
  call RuntimeParameters_get('smallp', sim_smallp)
  call RuntimeParameters_get('smalle', sim_smalle)
  call RuntimeParameters_get('smallx', sim_smallx)

  !call RuntimeParameters_get('sim_pAmbient', sim_pAmbient)
  call RuntimeParameters_get('sim_Tcore', sim_Tcore)
  call RuntimeParameters_get('sim_Tout', sim_Tout)
  call RuntimeParameters_get('sim_rhoCore', sim_rhoCore)
  call RuntimeParameters_get('sim_rhoFloor', sim_rhoFloor)
  call RuntimeParameters_get('sim_mu', sim_mu)
  call RuntimeParameters_get('sim_windVel', sim_windVel)
  call RuntimeParameters_get('sim_gammaICM', sim_gamma)
  call RuntimeParameters_get('sim_bzAmbient', sim_bzAmbient)
  call RuntimeParameters_get('sim_densityProfile', sim_densityProfile)
  call RuntimeParameters_get('sim_rCore', sim_rCore)
  call RuntimeParameters_get('sim_rCoreT', sim_rCoreT)
  call RuntimeParameters_get('sim_densityBeta', sim_densityBeta)
  sim_rCut = sim_rCore*sqrt((sim_rhoFloor/sim_rhoCore)**(-2./3./sim_densityBeta)-1.0 )
  call RuntimeParameters_get('sim_powerJet', sim(nozzle)%power)
  !call RuntimeParameters_get('sim_rhoJet', sim(nozzle)%density)
  call RuntimeParameters_get('sim_velJet', sim(nozzle)%velocity)
  call RuntimeParameters_get('sim_velJet', sim(nozzle)%velJet)
  call RuntimeParameters_get('sim_machJet', sim(nozzle)%mach)
  call RuntimeParameters_get('sim_initMachJet', sim(nozzle)%initMach)
  call RuntimeParameters_get('sim_outflowRatio', sim(nozzle)%outflowR)
  call RuntimeParameters_get('sim_gammaJet', sim(nozzle)%gamma)
  call RuntimeParameters_get('sim_betaJet', sim(nozzle)%beta)
  call RuntimeParameters_get('sim_helicityJet', sim(nozzle)%helicity)
  call RuntimeParameters_get('sim_timeMHDon', sim(nozzle)%timeMHDon)
  call RuntimeParameters_get('sim_tOn', sim(nozzle)%tOn)
  call RuntimeParameters_get('sim_duration', sim(nozzle)%duration)
  call RuntimeParameters_get('nozzleRadius', sim(nozzle)%radius)
  call RuntimeParameters_get('nozzleHalfL', sim(nozzle)%length)
  call RuntimeParameters_get('zTorInj', sim(nozzle)%zTorInj)
  call RuntimeParameters_get('rFeatherIn', sim(nozzle)%rFeatherIn)
  call RuntimeParameters_get('rFeatherOut', sim(nozzle)%rFeatherOut)
  sim(nozzle)%rFeatherMix = gr_minCellSize*NGUARD
  call RuntimeParameters_get('zFeather', sim(nozzle)%zFeather)
  sim(nozzle)%zFeatherMix = gr_minCellSize
  call RuntimeParameters_get('initGeometry', sim(nozzle)%initGeometry)
  !call RuntimeParameters_get('derefine_z1', sim(nozzle)%derefine_z1)
  !call RuntimeParameters_get('derefine_z2', sim(nozzle)%derefine_z2)
  call RuntimeParameters_get('refine_jetR', sim(nozzle)%refine_jetR)
  call RuntimeParameters_get('derefine_jetR', sim(nozzle)%derefine_jetR)
  call RuntimeParameters_get('lrefine_0', sim(nozzle)%lrefine_0)
  call RuntimeParameters_get('sim_ptInitNum', sim_ptInitNum)
  call RuntimeParameters_get('sim_ptAddPeriod', sim_ptAddPeriod)
  call RuntimeParameters_get('sim_ptAddArea', sim_ptAddArea)
  call RuntimeParameters_get('sim_ptSmljet', sim_ptSmljet)
  call RuntimeParameters_get('sim_ptMaxRadius', sim_ptMaxRadius)
  call RuntimeParameters_get('sim_ptRemoveDigit', sim_ptRemoveDigit)
  call Runtimeparameters_get('nozzlePrecession', sim(nozzle)%precession)
  maxPrecession = 0.1*sim(nozzle)%velocity/(1.0/sim(nozzle)%duration &
                     *(sim(nozzle)%length+sim(nozzle)%zFeatherMix))

  !if (dr_globalMe==MASTER_PE) then
  !   write(*,'(A24,f8.3)') 'maxPrecession is', maxPrecession
  !endif

  if (sim(nozzle)%precession.gt.maxPrecession) then
     sim(nozzle)%precession = maxPrecession
     if (sim_meshMe==MASTER_PE) then
        print*, '!!!!!!!!'
        print*, 'Warning! nozzlePrecession is too large.'
        write(*,'(A24,f8.3)') 'nozzlePrecession is now ', maxPrecession
     endif
  endif

  call RuntimeParameters_get('nozzleNutation', sim(nozzle)%nutation)
  call RuntimeParameters_get('nozzlePrecAngle', sim(nozzle)%precangle)
  call RuntimeParameters_get('useTableJiggle', sim_useTableJiggle)

  if ((sim_meshMe == MASTER_PE) .and. (sim_useTableJiggle == .true.)) then
     call RuntimeParameters_get('nozzleVecInput', sim_nozVecInput)
     inquire(file=sim_nozVecInput, exist=iexist)
     if (.not.iexist) then
        call Driver_abortFlash('[Simulation_init] ERROR: opening nozzle vector file')
     endif
  endif

  call RuntimeParameters_get('lowerRefHalf', sim_lowerRefHalf)


  if (dr_restart) then
     call IO_getScalar('coneVecX', sim(nozzle)%coneVec(1))
     call IO_getScalar('coneVecY', sim(nozzle)%coneVec(2))
     call IO_getScalar('coneVecZ', sim(nozzle)%coneVec(3))
     call IO_getScalar('nozzlePosX', sim(nozzle)%pos(1))
     call IO_getScalar('nozzlePosY', sim(nozzle)%pos(2))
     call IO_getScalar('nozzlePosZ', sim(nozzle)%pos(3))
     call IO_getScalar('nozzleVecX', sim(nozzle)%jetvec(1))
     call IO_getScalar('nozzleVecY', sim(nozzle)%jetvec(2))
     call IO_getScalar('nozzleVecZ', sim(nozzle)%jetvec(3))
     call IO_getScalar('nozzleVecOldX', sim(nozzle)%jetvecOld(1))
     call IO_getScalar('nozzleVecOldY', sim(nozzle)%jetvecOld(2))
     call IO_getScalar('nozzleVecOldZ', sim(nozzle)%jetvecOld(3))
     call IO_getScalar('nozzleAngVelX', sim(nozzle)%angVel(1))
     call IO_getScalar('nozzleAngVelY', sim(nozzle)%angVel(2))
     call IO_getScalar('nozzleAngVelZ', sim(nozzle)%angVel(3))
     call IO_getScalar('nozzleLinVelX', sim(nozzle)%linVel(1))
     call IO_getScalar('nozzleLinVelY', sim(nozzle)%linVel(2))
     call IO_getScalar('nozzleLinVelZ', sim(nozzle)%linVel(3))
     call IO_getScalar('nozzlePressure', sim(nozzle)%pressure)
     call IO_getScalar('nozzleDensity', sim(nozzle)%density)
     call IO_getScalar('nozzleBz', sim(nozzle)%bz)
     call IO_getScalar('nozzleBphi', sim(nozzle)%bphi)
     call IO_getScalar('randomSeed', sim(nozzle)%randSeed(1))

     call RuntimeParameters_get('resetOutputTime', resetOutputTime)

     if (resetOutputTime) then
        io_nextCheckpointTime = dr_simTime + io_checkpointFileIntervalTime
        io_nextPlotFileTime = dr_simTime + io_plotFileIntervalTime
        io_nextParticleFileTime = dr_simTime + io_particleFileIntervalTime
     endif


  else

     call RuntimeParameters_get('coneVecX', sim(nozzle)%coneVec(1))
     call RuntimeParameters_get('coneVecY', sim(nozzle)%coneVec(2))
     call RuntimeParameters_get('coneVecZ', sim(nozzle)%coneVec(3))
     sim(nozzle)%coneVec = sim(nozzle)%coneVec/ sqrt(sum(sim(nozzle)%coneVec*sim(nozzle)%coneVec))
     call RuntimeParameters_get('nozzlePosX', sim(nozzle)%pos(1))
     call RuntimeParameters_get('nozzlePosY', sim(nozzle)%pos(2))
     call RuntimeParameters_get('nozzlePosZ', sim(nozzle)%pos(3))
     call RuntimeParameters_get('nozzleVecX', sim(nozzle)%jetvec(1))
     call RuntimeParameters_get('nozzleVecY', sim(nozzle)%jetvec(2))
     call RuntimeParameters_get('nozzleVecZ', sim(nozzle)%jetvec(3))
     sim(nozzle)%jetvec = sim(nozzle)%jetvec/ sqrt(sum(sim(nozzle)%jetvec*sim(nozzle)%jetvec))
     sim(nozzle)%jetvecOld = sim(nozzle)%jetvec
     !call RuntimeParameters_get('nozzleAngVelX', sim(nozzle)%angVel(1))
     !call RuntimeParameters_get('nozzleAngVelY', sim(nozzle)%angVel(2))
     !call RuntimeParameters_get('nozzleAngVelZ', sim(nozzle)%angVel(3))
     sim(nozzle)%angVel(:) = (/0.0, 0.0, 0.0/)
     call RuntimeParameters_get('nozzleLinVelX', sim(nozzle)%linVel(1))
     call RuntimeParameters_get('nozzleLinVelY', sim(nozzle)%linVel(2))
     call RuntimeParameters_get('nozzleLinVelZ', sim(nozzle)%linVel(3))
     sim(nozzle)%density = -1.0
     ! Use the jet on time for initialization
     call sim_jetNozzleUpdate(nozzle, sim(nozzle)%tOn, dr_dt)
     call RuntimeParameters_get('randomSeed', sim(nozzle)%randSeed(1))

     if (sim(nozzle)%zTorInj < sim(nozzle)%length+1.5*sim(nozzle)%zFeather .and.&
        sim_meshMe==MASTER_PE) then
        print*, '!!!!!!!!'
        print*, 'Warning! zTorInj is too small that it overlaps with the nozzle.'
        print*, 'Toroidal field will be smaller than it should be.'
        print*, '!!!!!!!!'
     endif
  endif

  if (sim_meshMe==MASTER_PE) then
     write(*,'(a, 2es11.3, f7.2)') '(p, rho, M)=', &
     sim(nozzle)%pressure, sim(nozzle)%density, &
     sim(nozzle)%velocity/sqrt(sim(nozzle)%gamma*sim(nozzle)%pressure/sim(nozzle)%density)

     !write(*,'(a, 2es11.3)') '(bz, bphi)=', sim(nozzle)%bz, sim(nozzle)%bphi
     write(*,'(a, es11.3)') ' rFeatherOut:' , sim(nozzle)%rFeatherOut
     write(*,'(a, es11.3)') ' rFeatherMix:' , sim(nozzle)%rFeatherMix
     write(*,'(a, es11.3)') ' zFeatherMix:' , sim(nozzle)%zFeatherMix
     write(*,'(a, es11.3)') ' sim_rCut:' , sim_rCut
  endif


end subroutine Simulation_init
