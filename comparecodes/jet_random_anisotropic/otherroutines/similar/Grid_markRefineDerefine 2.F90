!!****if* source/Grid/GridMain/paramesh/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  Grid_markRefineDerefine(integer(IN) :: myPE)
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!! ARGUMENTS
!!  myPE : my processor id.
!! 
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_myPE or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine(MyPE)

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var,gr_refineOnParticleCount
  use tree, ONLY : newchild, refine, derefine, stay, lrefine_min, lrefine_max
  use Grid_interface   !, ONLY : Grid_fillGuardCells

  use Simulation_data

  implicit none

#include "constants.h"
#include "Flash.h"

  integer,intent(IN) :: MyPE
  
  real :: ref_cut,deref_cut,ref_filter

  real :: ictr, jctr, kctr, RefRad

  integer       :: l,i,iref
  logical :: doEos=.true.
  integer,parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  logical,dimension(maskSize) :: gcMask

  ! that are implemented in this file need values in guardcells

  gcMask=.false.
  do i = 1,gr_numRefineVars
     iref = gr_refine_var(i)
     if (iref > 0) gcMask(iref) = .TRUE.
  end do

  call Grid_fillGuardCells(MyPE,CENTER_FACES,ALLDIR,doEos=.true.,&
       maskSize=maskSize, mask=gcMask, makeMaskConsistent=.true.)

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  do l = 1,gr_numRefineVars
     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)
     ref_filter = gr_refine_filter(l)
     call gr_markRefineDerefine(MyPE,iref,ref_cut,deref_cut,ref_filter)
  end do

#ifdef FLASH_GRID_PARAMESH2
  ! Make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(MyPE,-1, 0.0, 0.0, 0.0)
  end if
#endif

  if(gr_refineOnParticleCount)call gr_ptMarkRefineDerefine()


  ictr = 0.0  !0.5*(sim_xMax-sim_xMin)
  jctr = 0.0  !0.5*(sim_yMax-sim_yMin)
  kctr = 0.0  !0.5*(sim_zMax-sim_zMin)
!  RefRad = 0.75*ictr
!  RefRad = 0.50*ictr
!  RefRad = 0.3333333*ictr
  RefRad = 0.2*sim_xMax  ! 100kpc


!  call gr_markInRadius(ictr, jctr, kctr, RefRad, lrefine_min)
  call gr_markInRadius(ictr, jctr, kctr, RefRad, lrefine_min, lrefine_max)

  
  return
end subroutine Grid_markRefineDerefine

