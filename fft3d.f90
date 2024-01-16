module fft3d

use stafft

implicit none

contains

  subroutine ptospc(n,var,trig,fac)

    integer, parameter :: dp = kind(1.d0) 
    real(dp), dimension(n,n,n) :: var
    real(dp), dimension(n) :: a
    real(dP), dimension(2*n) :: trig
    integer, dimension(5) :: fac
    integer :: ix,iy,iz
    integer :: n

!   FFT over x   
    do iz=1,n
       do iy=1,n
          !$omp parallel do private(ix) shared(var,a)
          do ix=1,n
             a(ix) = var(ix,iy,iz)
          enddo 
          !somp end parallel do
          call forfft(1,n,a,trig,fac)
          !$omp parallel do private(ix) shared(var,a)         
          do ix=1,n
             var(ix,iy,iz) = a(ix)             
          enddo
          !s omp end parallel do
       enddo
    enddo  

!    FFT over y
    do iz=1,n
       do ix=1,n
          !somp parallel do private(iy) shared(var,a)
          do iy=1,n
             a(iy) = var(ix,iy,iz)
          enddo   
          !somp end parallel do
          call forfft(1,n,a,trig,fac)
          !somp parallel do private(iy) shared(var,a)
          do iy=1,n
             var(ix,iy,iz) = a(iy)
          enddo
         !somp end parallel do
       enddo   
    enddo
   
!    FFT over z
    do ix=1,n
       do iy=1,n
          !somp parallel do private(iz) shared(var,a)
          do iz=1,n
             a(iz) = var(ix,iy,iz)
          enddo   
          !somp end parallel do
          call forfft(1,n,a,trig,fac)
          !somp parallel do private(iz) shared(var,a)
          do iz=1,n
             var(ix,iy,iz) = a(iz)
          enddo
         !somp end parallel do
       enddo   
    enddo

  return
    
  end subroutine ptospc  

  subroutine spctop(n,var,trig,fac)

    integer, parameter :: dp = kind(1.d0) 
    real(dp), dimension(n,n,n) :: var
    real(dp), dimension(n) :: a
    real(dP), dimension(2*n) :: trig
    integer, dimension(5) :: fac
    integer :: ix,iy,iz
    integer :: n

!   FFT over x   
    do iz=1,n
       do iy=1,n
          !$omp parallel do private(ix) shared(var,a)
          do ix=1,n
             a(ix) = var(ix,iy,iz)
          enddo 
          !somp end parallel do
          call revfft(1,n,a,trig,fac)
          !$omp parallel do private(iy) shared(var,a)         
          do ix=1,n
             var(ix,iy,iz) = a(ix)             
          enddo
          !somp end parallel do
       enddo
    enddo  

!    FFT over y
    do iz=1,n
       do ix=1,n
          !somp parallel do private(iy) shared(var,a)
          do iy=1,n
             a(iy) = var(ix,iy,iz)
          enddo   
          !somp end parallel 
          call revfft(1,n,a,trig,fac)
          !somp parallel do private(iy) shared(var)
          do iy=1,n
             var(ix,iy,iz) = a(iy)
          enddo
         !somp end parallel do
       enddo   
    enddo
   
!    FFT over z
    do ix=1,n
       do iy=1,n
          !somp parallel do private(iy) shared(var,a)
          do iz=1,n
             a(iz) = var(ix,iy,iz)
          enddo   
          !somp end parallel do
          call revfft(1,n,a,trig,fac)
          !somp parallel do private(iz) shared(var,a)
          do iz=1,n
             var(ix,iy,iz) = a(iz)
          enddo
         !somp end parallel do
       enddo   
    enddo

   return    
  end subroutine spctop  

end module fft3d
