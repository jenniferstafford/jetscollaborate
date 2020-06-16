!!****if* source/Simulation/SimulationMain/magnetoHD/Random_field_test/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for the random field generator
!!
!! ARGUMENTS
!!
!!
!!***

module Simulation_data

  implicit none
#include "Simulation.h" 
  !! *** Runtime Parameters *** !!

  real, save, allocatable, dimension(:,:,:) :: sim_Ax, sim_Ay, sim_Az
  real, save, dimension(NXAVEC) :: sim_ax_coords
  real, save, dimension(NYAVEC) :: sim_ay_coords
  real, save, dimension(NZAVEC) :: sim_az_coords
  real, save :: sim_xmin, sim_xmax, sim_ymin, sim_ymax, sim_zmin, sim_zmax
  real, save, dimension(64,64) :: sim_cubic_matrix
  real, save :: sim_bcentral
  integer, save :: psd_seed

end module Simulation_data



