!!****if* source/Simulation/SimulationMain/Sedov/Simulation_data
!!
!! NAME
!!
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data 
!!
!!  DESCRIPTION
!!
!!  Stores the local data for Simulation setup: Sedov
!!  
!! PARAMETERS
!!

!!
!!***

module Simulation_data

  implicit none
#include "constants.h"

  !! *** Runtime Parameters *** !!

  !! *** Variables pertaining to this Simulation *** !!

  real, save    :: sim_gamma, sim_zinit, sim_hubble
  real, save    :: sim_xMin, sim_xMax, sim_yMin, sim_yMax, sim_zMin, sim_zMax, sim_eos_singleSpeciesA
  real, save    :: sim_Bfield0, T_in, T_out, rho_in 
  real, save    :: kb, mp, newton, Msun, Mpc, kpc, km
  logical, save :: sim_killdivb
  real, save    :: h, Mvir, c, mc, M0, rvir, rs, DeltaT, tbx, tby, tbz, dr, rc


  ! header variables  

  integer, save :: ntotal, ntot, ngas, ntot_tmp

  real(kind=4), allocatable, save    :: Positions_Ascii(:,:), Velocities_Ascii(:,:)
  real(kind=8), allocatable, save    :: Masses_Ascii(:)
  real(kind=4), allocatable, save    :: Positions_Ascii_tmp(:,:), Velocities_Ascii_tmp(:,:)
  real(kind=8), allocatable, save    :: Masses_Ascii_tmp(:)


  real, allocatable, save    :: bx(:,:,:), by(:,:,:), bz(:,:,:), rhogas(:), tempgas(:)

  logical :: sim_saturatedConduction


end module Simulation_data



