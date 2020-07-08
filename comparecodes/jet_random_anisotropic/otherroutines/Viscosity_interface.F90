!!****h* source/physics/materialProperties/Viscosity/Viscosity_interface
!!
!! NAME    
!!   Viscosity_interface
!!
!! SYNOPSIS
!!
!!  use Viscosity_interface, ONLY : Viscosity
!!
!! DESCRIPTION
!!
!! This is the header file for the Viscosity Unit that defines its
!! public interfaces.
!!***

#include "Flash.h"

module Viscosity_interface
  interface

     subroutine Viscosity_init(myPE)
       integer, intent(IN) :: myPE 
     end subroutine Viscosity_init
  end interface

  interface
     subroutine Viscosity_finalize()
     end subroutine Viscosity_finalize
  end interface
  
  interface
     subroutine Viscosity(xtemp,xden,massfrac,visc_dyn,visc)
       real,INTENT(out)   :: visc, visc_dyn
       real,INTENT(in)    :: xtemp
       real,INTENT(in)    :: xden
       real,INTENT(in)    :: massfrac(NSPECIES)
     end subroutine Viscosity
  end interface

end module Viscosity_interface
