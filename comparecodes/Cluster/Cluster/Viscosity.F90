!!****f* source/physics/materialProperties/Viscosity/Viscosity
!!
!!  NAME    
!!   Viscosity
!!
!!  SYNOPSIS
!!   Viscosity(real, intent(IN)  :: xtemp, 
!!             real, intent(IN)  :: xden, 
!!             real, inteng(IN)  :: massfrac(NSPECIES),
!!             real, intent(OUT) :: visc)
!!
!!
!! DESCRIPTION
!!   A generic viscosity routine.  Just returns a constant viscosity
!!   from the run-time parameter, 'diff_visc_nu'
!!
!! ARGUMENTS
!!
!!  INPUTS
!!   xtemp   :   REAL    temperature (in K)
!!   xden    :   REAL    density (in g/cm**3)
!!   massfrac:   REAL    mass fractions of the composition
!!
!!  OUTPUTS
!!   visc    :   REAL    viscosity
!!
!!***

subroutine Viscosity(xtemp,xden,massfrac,visc_dyn,visc)
!! True Stub
  implicit none

  real,INTENT(in)    :: xtemp
  real,INTENT(in)    :: xden
  real,INTENT(in)    :: massfrac
  real,INTENT(out)   :: visc, visc_dyn

  !dummy value assigned in stub
  visc_dyn = 0.
  visc = 0.

  return 
end subroutine Viscosity
