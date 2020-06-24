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
!!  Initializes all the data specified in Simulation_data.
!!  It calls RuntimeParameters_get rotuine for initialization.
!!  Initializes initial conditions for Sedov Spherical Explosion 
!!  problem.
!!Initializes all the parameters needed for the random
!!field generator
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!
!!***

subroutine Simulation_init()

  use Hydro_data, ONLY : hy_units
  use Simulation_data
  use Simulation_jetNozzleUpdate, ONLY : sim_jetNozzleUpdate
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_data, ONLY : dr_simTime, dr_dt, dr_restart
  use Driver_interface, ONLY : Driver_getMype,  Driver_abortFlash
  use Grid_data, ONLY : gr_minCellSize
  use IO_interface, ONLY :  IO_getScalar
  use Particles_data, ONLY : pt_randSeed
  use IO_data, ONLY: io_nextCheckpointTime, io_checkpointFileIntervalTime,\
                     io_nextPlotFileTime, io_plotFileIntervalTime
  use IOParticles_data, ONLY: io_nextParticleFileTime, io_particleFileIntervalTime
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get

  
  !use Grid_data, ONLY : gr_smallrho

  implicit none
#include "constants.h"
#include "Flash.h"
#include "/usr/fftw-3.3.4/include/fftw3.f"
#include "Simulation.h"

  integer :: nozzle=1, clock, iexist
  real    :: maxPrecession
  logical :: resetOutputTime
  !cluster
  integer, intent(in) :: myPE
  integer :: i, ii
  real :: tmp
  !random
  real :: rn1, rn2, bp, seedin, delr, del, rad, norm
  complex, dimension(NXAVEC,NYAVEC,NZAVEC) :: arr
  integer :: plan, i, j, k, ii, jj, kk, istat
  integer, dimension(1) :: seed
  real :: index1,index2,index3,k1,k2,knull,delx,dely,delz
  real :: xmin, xmax, ymin, ymax, zmin, zmax
  real, dimension(4) :: indices
  real, dimension(4,4) :: one_matrix

  knull=DBLE(NXAVEC*NYAVEC*NZAVEC)**(1.d0/3.d0)/2.0d0
  seedin=1021345499345.


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

     !from cluster
     

  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('xmin',sim_xMin)
  call RuntimeParameters_get('xmax',sim_xMax)
  call RuntimeParameters_get('ymin',sim_yMin)
  call RuntimeParameters_get('ymax',sim_yMax)
  call RuntimeParameters_get('zmin',sim_zMin)
  call RuntimeParameters_get('zmax',sim_zMax)
  call RuntimeParameters_get('eos_singleSpeciesA', sim_eos_singleSpeciesA)
  call RuntimeParameters_get('killdivb',sim_killdivb)
  call RuntimeParameters_get('Bfield0', sim_Bfield0)
  call RuntimeParameters_get('saturated', sim_saturatedConduction)
!  call RuntimeParameters_get('T_in',T_in)
!  call RuntimeParameters_get('T_out',T_out)
!  call RuntimeParameters_get('rho_in',rho_in)
!  call RuntimeParameters_get('Mvir',Mvir)
  call RuntimeParameters_get('M0',M0)
  call RuntimeParameters_get('rs',rs)
  call RuntimeParameters_get('rc',rc)

  call PhysicalConstants_get('Boltzmann', kb)
  call PhysicalConstants_get('Newton', newton)
  call PhysicalConstants_get('proton mass', mp)

!  h = 0.7e0
!  Mvir = Mvir/h
!  c  = 6.0* ((Mvir/(1.0e14/h))**(-0.2e0))
!  mc = log(1.0+c) - c/(1.0+c)
!  M0 = (Mvir*1.9889225e33)/(2.0*mc)
!  rvir = 4.95e24/h *(Mvir*h/5.0e14)**(1.0/3.0)
!  rs = rvir/c


  Msun = 1.9889225E33
  Mpc  = 3.0856775807E24
  km   = 1.0e5 
  kpc  = 3.0856775807E21

  rs = rs * kpc
  rc = rc * kpc
  M0 = M0 * Msun


  DeltaT = T_out-T_in

  !random
  
  call RuntimeParameters_get('high_cut_index', index1)
  call RuntimeParameters_get('low_cut_index', index2)
  call RuntimeParameters_get('psd_index', index3)
  call RuntimeParameters_get('high_cut_k', k1)
  call RuntimeParameters_get('low_cut_k', k2)
!  call RuntimeParameters_get('psd_seed', seedin)
  call RuntimeParameters_get('sim_bcentral', sim_bcentral)
  
  seed(1)=int(seedin)

  call RANDOM_SEED
  call RANDOM_SEED(SIZE=K)
  call RANDOM_SEED(PUT=seed(1:K))

  if (dr_restart) then

  else
  
     allocate(sim_Ax(NXAVEC,NYAVEC,NZAVEC),stat=istat)  
     allocate(sim_Ay(NXAVEC,NYAVEC,NZAVEC),stat=istat)  
     allocate(sim_Az(NXAVEC,NYAVEC,NZAVEC),stat=istat)  

     delx=(sim_xmax-sim_xmin)/DBLE(NXAVEC-1)
     dely=(sim_ymax-sim_ymin)/DBLE(NYAVEC-1)
     delz=(sim_zmax-sim_zmin)/DBLE(NZAVEC-1)
     del=(delx*dely*delz)**(1.d0/3.d0)
     sim_ax_coords=(/ (i, i=0, NXAVEC-1) /)*delx + sim_xmin
     sim_ay_coords=(/ (i, i=0, NYAVEC-1) /)*dely + sim_ymin
     sim_az_coords=(/ (i, i=0, NZAVEC-1) /)*delz + sim_zmin
     
     ! define our FFT: Inverse, size is NXAVEC,NYAVEC,NZAVEC
     call dfftw_plan_dft_3d(plan, NXAVEC,NYAVEC,NZAVEC, arr, arr, &
          FFTW_BACKWARD, FFTW_ESTIMATE)

     ! X-component of vector potential
     do i=1,NXAVEC
        do j=1,NYAVEC
           do k=1,NZAVEC

              call RANDOM_NUMBER(rn1)
              call RANDOM_NUMBER(rn2)
              rn1=2.*rn1 - 1.
              rn2=2.*rn2 - 1.
              arr(i,j,k)=CMPLX(rn1,rn2)

              rad=sqrt((i-1.)**2 + (j-1.)**2 + (k-1.)**2)
              if (rad.eq.0.) then
                  rad=1.0
              endif

              ! Normalize to power spectrum
              norm=sqrt(arr(i,j,k)*conjg(arr(i,j,k)))
              arr(i,j,k)=arr(i,j,k)/norm * rad**(-index3) * & 
                   exp(-(rad/k1)**index1)*exp(-(k2/rad)**index2) * &
                   (SIGN(0.5d0,1.d0-rad/knull)+0.5d0)
              
           enddo
        enddo
     enddo
     ! No constant term, so average should be zero
     arr(1,1,1)=CMPLX(0.,0.)

     call dfftw_execute_dft(plan, arr, arr)
     
     ! real part of inverse transform is Ax.
     sim_Ax=REAL(arr)*del

     ! Y-component of vector potential
     do i=1,NXAVEC
        do j=1,NYAVEC
           do k=1,NZAVEC

              call RANDOM_NUMBER(rn1)
              call RANDOM_NUMBER(rn2)
              rn1=2.*rn1 - 1.
              rn2=2.*rn2 - 1.
              arr(i,j,k)=CMPLX(rn1,rn2)

              rad=sqrt((i-1.)**2 + (j-1.)**2 + (k-1.)**2)
              if (rad.eq.0.) then
                  rad=1.0
              endif

              norm=sqrt(arr(i,j,k)*conjg(arr(i,j,k)))
              arr(i,j,k)=arr(i,j,k)/norm*rad**(-index3)* & 
                   exp(-(rad/k1)**index1)*exp(-(k2/rad)**index2) * &
                   (SIGN(0.5d0,1.d0-rad/knull)+0.5d0)

           enddo
        enddo
     enddo
     arr(1,1,1)=CMPLX(0,0)

     call dfftw_execute_dft(plan, arr, arr)
     sim_Ay=REAL(arr)*del

     ! Z-component of vector potential
     do i=1,NXAVEC
        do j=1,NYAVEC
           do k=1,NZAVEC
              call RANDOM_NUMBER(rn1)
              call RANDOM_NUMBER(rn2)
              rn1=2.*rn1 - 1.
              rn2=2.*rn2 - 1.
              arr(i,j,k)=CMPLX(rn1,rn2)

              rad=sqrt((i-1.)**2 + (j-1.)**2 + (k-1.)**2)
              if (rad.eq.0.) then
                  rad=1.0
              endif
              
              norm=sqrt(arr(i,j,k)*conjg(arr(i,j,k)))
              arr(i,j,k)=arr(i,j,k)/norm*rad**(-index3)* & 
                   exp(-(rad/k1)**index1)*exp(-(k2/rad)**index2) * &
                   (SIGN(0.5d0,1.d0-rad/knull)+0.5d0)

           enddo
        enddo
     enddo
     arr(1,1,1)=CMPLX(0,0)

     call dfftw_execute_dft(plan, arr, arr)
     sim_Az=REAL(arr)*del

     ! Normalize to desired plasma beta...
     bp=0.
     do i=1,NXAVEC
        ii=MODULO(i,NXAVEC) + 1
        do j=1,NYAVEC
           jj=MODULO(j,NXAVEC) + 1
           do k=1,NZAVEC
              kk=MODULO(k,NXAVEC) + 1
              bp=bp + 0.25 * &
                   (((sim_Az(ii,jj,k)-sim_Az(ii,j,k))/dely-(sim_Ay(ii,j,kk)-sim_Ay(ii,j,k))/delz + &
                     (sim_Az(i,jj,k) -sim_Az(i,j,k)) /dely-(sim_Ay(i,j,kk) -sim_Ay(i,j,k)) /delz)**2 + &
                    ((sim_Ax(i,jj,kk)-sim_Ax(i,jj,k))/delz-(sim_Az(ii,jj,k)-sim_Az(i,jj,k))/delx + &
                     (sim_Ax(i,j,kk) -sim_Ax(i,j,k)) /delz-(sim_Az(ii,j,k) -sim_Az(i,j,k)) /delx)**2 + &
                    ((sim_Ay(ii,j,kk)-sim_Ay(i,j,kk))/delx-(sim_Ax(i,jj,kk)-sim_Ax(i,j,kk))/dely + &
                     (sim_Ay(ii,j,k) -sim_Ay(i,j,k)) /delx-(sim_Ax(i,jj,k) -sim_Ax(i,j,k)) /dely)**2)
           enddo
        enddo
     enddo

     bp=sqrt(bp/NXAVEC/NYAVEC/NZAVEC)
     
     sim_Ax=sim_Ax/bp*sim_bcentral
     sim_Ay=sim_Ay/bp*sim_bcentral
     sim_Az=sim_Az/bp*sim_bcentral

     call dfftw_destroy_plan(plan)

     ! Set up matrix for cubic interpolation

     one_matrix(1,:) = (/  1.d0/2.d0  ,  -3.d0/2.d0  ,   3.d0/2.d0  , -1.d0/2.d0  /)
     one_matrix(2,:) = (/  1.d0/4.d0  ,  -1.d0/4.d0  ,  -1.d0/4.d0  ,  1.d0/4.d0  /)
     one_matrix(3,:) = (/ -1.d0/8.d0  ,  11.d0/8.d0  , -11.d0/8.d0  ,  1.d0/8.d0  /)
     one_matrix(4,:) = (/ -1.d0/16.d0 ,   9.d0/16.d0 ,   9.d0/16.d0 , -1.d0/16.d0 /)
     
     do i = 1,4
        do j = 1,4
           do k = 1,4
              do ii = 1,4
                 do jj = 1,4
                    do kk = 1,4
                       sim_cubic_matrix(16*(ii-1) + 4*(jj-1) + kk, 16*(i-1) + 4*(j-1) + k) = &
                            & one_matrix(ii,i)*one_matrix(jj,j)*one_matrix(kk,k)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
 
  end if



  
!jet
     
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

  !from cluster
  
!  open(unit=1,file='bx.dat',status='old',form='unformatted')
  open(unit=1,file='bx.dat',status='old')
  read(1,*)ntotal, tbx
  allocate(bx(ntotal,ntotal,ntotal))
  read(1,*) bx
  close (1)

!  open(unit=1,file='by.dat',status='old',form='unformatted')
  open(unit=1,file='by.dat',status='old')
  read(1,*)ntotal, tby
  allocate(by(ntotal,ntotal,ntotal))
  read(1,*) by
  close (1)

!  open(unit=1,file='bz.dat',status='old',form='unformatted')
  open(unit=1,file='bz.dat',status='old')
  read(1,*)ntotal, tbz
  allocate(bz(ntotal,ntotal,ntotal))
  read(1,*) bz
  close (1)

! write(*,*) bx(3,1,1), bx(8,2,120), bx(4,50,4)
! stop

  !reading in the ICs for the gas
  open(unit=1,file='A2199_gas_profile_linear.dat',status='old')
  read(1,*) ngas, dr
  allocate(rhogas(ngas))
  allocate(tempgas(ngas))
  do i = 1, ngas
     read(1,*) tmp, rhogas(i), tempgas(i) 
  enddo
  close (1)

end subroutine Simulation_init


  

end subroutine Simulation_init
