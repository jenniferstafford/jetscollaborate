!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/Spitzer/Conductivity_init
!!
!! NAME
!!
!!  Conductivity_init
!!
!! SYNOPSIS
!!
!!  Conductivity_init(integer(IN) :: myPE)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   myPE -- local processor number
!!
!!
!!
!!***

subroutine Conductivity_init(myPE)

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get

  use Conductivity_data, ONLY: cond_myPE, cond_numProcs, SpitzerFraction, isotropic
  use Grid_interface, ONLY: Grid_getNumProcs
  implicit none

  integer, intent(IN) :: myPE

  ! Everybody should know this
  cond_myPE = myPE
  call Grid_getNumProcs(cond_numProcs)

  call RuntimeParameters_get("SpitzerFraction", SpitzerFraction)
  call RuntimeParameters_get("isotropic", isotropic)

end subroutine Conductivity_init
