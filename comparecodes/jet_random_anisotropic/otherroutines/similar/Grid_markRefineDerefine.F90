!!****if* source/Grid/GridMain/paramesh/Grid_markRefineDerefine
!!
!! NAME
!!  Grid_markRefineDerefine
!!
!! SYNOPSIS
!!
!!  call Grid_markRefineDerefine()
!!  
!! DESCRIPTION 
!!  Mark blocks for refinement or derefinement
!!  This routine is used with AMR only where individual 
!!  blocks are marked for refinement or derefinement based upon
!!  some refinement criterion. The Uniform Grid does not need
!!  this routine, and uses the stub.
!!
!!  This routine is normally called by the implementation of
!!  Grid_updateRefinement. It may also get called repeatedly
!!  during the initial construction of the Grid from
!!  Grid_initDomain.
!!
!! ARGUMENTS
!!
!!  none
!!
!! SEE ALSO
!!
!!  Grid_updateRefinement
!!  Grid_initDomain
!!  gr_expandDomain
!!
!! NOTES
!!
!! Every unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. For Grid unit these variables begin with "gr_"
!! like, gr_meshMe or gr_eosMode, and are stored in fortran
!! module Grid_data (in file Grid_data.F90). The other variables
!! are local to the specific routines and do not have the prefix "gr_"
!!
!!
!!***

subroutine Grid_markRefineDerefine()

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var,gr_refineOnParticleCount,&
                        gr_enforceMaxRefinement, gr_maxRefine,&
                        gr_lrefineMaxByTime,&
                        gr_lrefineMaxRedDoByTime,&
                        gr_lrefineMaxRedDoByLogR,&
                        gr_lrefineCenterI,gr_lrefineCenterJ,gr_lrefineCenterK,&
                        gr_eosModeNow
  !use tree, ONLY : newchild, refine, derefine, stay, nodetype
!!$  use physicaldata, ONLY : force_consistency
  use Logfile_interface, ONLY : Logfile_stampVarMask
  !use Grid_interface, ONLY : Grid_fillGuardCells
  use Particles_interface, only: Particles_sinkMarkRefineDerefine
!<-- ychen
  use tree, ONLY : newchild, refine, derefine, stay, nodetype, lrefine_max
  use Grid_interface, ONLY : Grid_fillGuardCells, Grid_markRefineSpecialized
  use Simulation_data
!--> ychen
  implicit none

#include "constants.h"
#include "Flash.h"
#include "Simulation.h"

  
  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref
  logical,save :: gcMaskArgsLogged = .FALSE.
  integer,save :: eosModeLast = 0
  logical :: doEos=.true.
  integer,parameter :: maskSize = NUNK_VARS+NDIM*NFACE_VARS
  logical,dimension(maskSize) :: gcMask
!<-- ychen
  integer :: nozzle
  real :: r_ref
!--> ychen

  if(gr_lrefineMaxRedDoByTime) then
     call gr_markDerefineByTime()
  end if
  
  if(gr_lrefineMaxByTime) then
     call gr_setMaxRefineByTime()
  end if

  if (gr_eosModeNow .NE. eosModeLast) then
     gcMaskArgsLogged = .FALSE.
     eosModeLast = gr_eosModeNow
  end if

  ! that are implemented in this file need values in guardcells

  gcMask=.false.
  do i = 1,gr_numRefineVars
     iref = gr_refine_var(i)
     if (iref > 0) gcMask(iref) = .TRUE.
  end do

  gcMask(NUNK_VARS+1:min(maskSize,NUNK_VARS+NDIM*NFACE_VARS)) = .TRUE.
!!$  gcMask(NUNK_VARS+1:maskSize) = .TRUE.


  if (.NOT.gcMaskArgsLogged) then
     call Logfile_stampVarMask(gcMask, .true., '[Grid_markRefineDerefine]', 'gcArgs')
  end if

!!$  force_consistency = .FALSE.
  call Grid_fillGuardCells(CENTER_FACES,ALLDIR,doEos=.true.,&
       maskSize=maskSize, mask=gcMask, makeMaskConsistent=.true.,doLogMask=.NOT.gcMaskArgsLogged,&
       selectBlockType=ACTIVE_BLKS)
     gcMaskArgsLogged = .TRUE.
!!$  force_consistency = .TRUE.

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  do l = 1,gr_numRefineVars
     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)
     ref_filter = gr_refine_filter(l)
     call gr_markRefineDerefine(iref,ref_cut,deref_cut,ref_filter)
  end do
  
!<-- ychen
  do nozzle=1, NOZZLES 
     ! Mark the jet using momentum criterion
     call gr_markJet(nozzle)

     ! Mark the nozzle region (a sphere) for maximum refinement
     r_ref = 5.0*max(sim(nozzle)%radius, sim(nozzle)%length)
     call Grid_markRefineSpecialized(INRADIUS, 4, (/ sim(nozzle)%pos(1), &
     sim(nozzle)%pos(2), sim(nozzle)%pos(3), r_ref/), lrefine_max )

     !call gr_markCylinder(sim(nozzle)%pos, sim(nozzle)%jetvec,&
     !                     sim(nozzle)%radius, sim(nozzle)%length, lrefine_max)

     !rvec = cross(sim(nozzle)%jetvec, (/0.0, gr_smallx, 1.0/))
     !rvec = rvec/sqrt(sum(rvec**2))
     !rb = sim(nozzle)%pos + sim(nozzle)%length*sim(nozzle)%jetvec &
     !                     + sim(nozzle)%radius*rvec &
     !                     + sim(nozzle)%radius*cross(sim(nozzle)%jetvec, rvec)
     !lb = sim(nozzle)%pos - sim(nozzle)%length*sim(nozzle)%jetvec &
     !                     - sim(nozzle)%radius*rvec &
     !                     - sim(nozzle)%radius*cross(sim(nozzle)%jetvec, rvec)

     !call gr_markInRectangle(min(lb(1), rb(1)), max(lb(1), rb(1)),&
     !                        min(lb(2), rb(2)), max(lb(2), rb(2)),&
     !                        min(lb(3), rb(3)), max(lb(3), rb(3)), lrefine_max, 0)

  end do
!--> ychen

#ifdef FLASH_GRID_PARAMESH2
  ! For PARAMESH2, call gr_markRefineDerefine here if it hasn't been called above.
  ! This is necessary to make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(-1, 0.0, 0.0, 0.0)
  end if
#endif

  if(gr_refineOnParticleCount)call gr_ptMarkRefineDerefine()

  if(gr_enforceMaxRefinement) call gr_enforceMaxRefine(gr_maxRefine)

  if(gr_lrefineMaxRedDoByLogR) &
       call gr_unmarkRefineByLogRadius(gr_lrefineCenterI,&
       gr_lrefineCenterJ,gr_lrefineCenterK)
  
  call Particles_sinkMarkRefineDerefine()

  ! When the flag arrays are passed to Paramesh for processing, only leaf
  ! blocks should be marked. - KW
  where (nodetype(:) .NE. LEAF)
     refine(:)   = .false.
     derefine(:) = .false.
  end where
  
  return
end subroutine Grid_markRefineDerefine

