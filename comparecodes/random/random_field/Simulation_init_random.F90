!!****if* source/Simulation/SimulationMain/magnetoHD/Random_field_test/Simulation_init
!!NAME Simulation_init SYNOPSIS Simulation_init()
!!DESCRIPTION Initializes all the parameters needed for the random
!!field generator ARGUMENTS ***

subroutine Simulation_init()

  use Hydro_data, ONLY : hy_units
  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Driver_data, ONLY : dr_restart
  use Driver_interface,  ONLY : Driver_getMype ! added 
  implicit none

#include "constants.h"
#include "Flash.h"
#include "/usr/fftw-3.3.4/include/fftw3.f"
#include "Simulation.h"

  real :: rn1, rn2, bp, seedin, delr, del, rad, norm
  complex, dimension(NXAVEC,NYAVEC,NZAVEC) :: arr
  integer :: plan,i, j, k, ii, jj, kk, istat
  integer, dimension(1) :: seed
  real :: index1,index2,index3,k1,k2,knull,delx,dely,delz
  real :: xmin, xmax, ymin, ymax, zmin, zmax

  real, dimension(4) :: indices
  real, dimension(4,4) :: one_matrix

  knull=DBLE(NXAVEC*NYAVEC*NZAVEC)**(1.d0/3.d0)/2.0d0
  seedin=1021345499345.
!  seedin=6453451021345.

  call RuntimeParameters_get('xmin', xmin)
  call RuntimeParameters_get('xmax', xmax)
  call RuntimeParameters_get('ymin', ymin)
  call RuntimeParameters_get('ymax', ymax)
  call RuntimeParameters_get('zmin', zmin)
  call RuntimeParameters_get('zmax', zmax)
  call RuntimeParameters_get('sim_xmin', sim_xmin)
  call RuntimeParameters_get('sim_xmax', sim_xmax)
  call RuntimeParameters_get('sim_ymin', sim_ymin)
  call RuntimeParameters_get('sim_ymax', sim_ymax)
  call RuntimeParameters_get('sim_zmin', sim_zmin)
  call RuntimeParameters_get('sim_zmax', sim_zmax)

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

end subroutine Simulation_init
