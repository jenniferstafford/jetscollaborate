!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/Spitzer/Conductivity
!!
!! NAME
!!
!!  Conductivity
!!
!! SYNOPSIS
!!
!!  Conductivity(real, intent(IN)  :: xtemp,
!!               real, intent(IN)  :: xden,
!!               real, dimension(NSPECIES), intent(IN)  :: massfrac,
!!               real, intent(OUT)  :: cond,
!!               real, intent(OUT)  :: diff_coeff)
!!
!! DESCRIPTION
!!
!!   Thermal conductivity and diffusion coefficient given by
!!   Spitzer (1962)
!!
!! ARGUMENTS
!!
!!   xtemp          temperature (in K)
!!   xden           density (in g/cm**3)
!!   massfrac       mass fractions of the composition
!!   cond           conductivity
!!   diff_coeff     diffusion coefficient ( = cond/(rho*cv))
!!
!!  NOTES
!!   See: Spitzer L. 1962, In: `Physics of fully ionized gases', 
!!   (New York: Wiley Interscience)
!!
!!***

subroutine Conductivity(xtemp,xden,massfrac,cond,diff_coeff)

  use Simulation_data, ONLY : sim_eos_singleSpeciesA
  use Conductivity_data, ONLY : SpitzerFraction
  use Hydro_data,        ONLY : hy_kref, hy_qref, hy_Rconst

  use Eos_interface, ONLY : Eos

  implicit none
  
#include "constants.h"  
#include "Flash.h"
#include "Eos.h"
  
  real, intent(IN) :: xtemp, xden
  real, intent(OUT) ::  diff_coeff, cond
  real, dimension(NSPECIES), intent(IN) :: massfrac
  real :: xtemp_tmp, n_e  
  logical :: cluster

  real, dimension(EOS_NUM+NSPECIES) :: eos_arr
!  logical, dimension(EOS_NUM) :: mask
  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask
  integer :: mode, vecLen
  
!  real, parameter :: ck = 9.2e-7
  real, parameter :: ck = 1.84e-5/40.0e0
!  real, parameter :: cexp = 2.5, K_max = 1.0e32
  real, parameter :: cexp = 2.5, K_max = 1.0e33  !3.0e33

  vecLen = 1
  mode = MODE_DENS_PRES
  eos_arr(EOS_TEMP) = xtemp
  eos_arr(EOS_DENS) = xden
  mask = .false.
!  mask(EOS_TEMP) = .true.
!  mask(EOS_DENS) = .true.
  mask(EOS_CV) = .true.
  mask(EOS_CP) = .true.
  mask(EOS_DET) = .true.

  call Eos(mode,vecLen,eos_arr,massfrac,mask)

  n_e = xden/1.67e-24  ! assuming pure ionized hydrogen

  cond = ck*xtemp**cexp
!  if (xtemp .gt. 1.0e8) cond = ck * (1.0e8**cexp)   ! temporary "saturation"

  diff_coeff = cond/(xden*eos_arr(EOS_CV))

! setting limit to K *unless* we are inside the cluster
!  cluster = (n_e .gt. 1.0e-5) .or. (xtemp .gt. 1.0e4)
!  if ( .not.cluster .and. ( diff_coeff .gt. K_max) ) cond = xden*eos_arr(EOS_CV)*K_max

!  if ( diff_coeff .gt. K_max ) cond = xden*eos_arr(EOS_CV)*K_max

  if ( diff_coeff .gt. K_max ) cond = xden*eos_arr(EOS_CV)*K_max

!  if (n_e .le. 1.0e-7)  cond = 0.0
!  if (xtemp .le. 1.0e3) cond = 0.0

  cond       = cond*SpitzerFraction  
  diff_coeff = cond/(xden*eos_arr(EOS_CV))

!  if (n_e .le. 1.0e-7)          diff_coeff = 1.0e-20
!  if (xtemp .le. 1.0e3)          diff_coeff = 1.0e-20

  return 
end subroutine Conductivity
