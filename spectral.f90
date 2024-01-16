module spectral
implicit none

! ========================================

contains

! ========================================

  subroutine initarr(n,rk,wk,rks,rksi,wkf)

    integer, parameter :: dp = kind(1.d0)    
    integer :: n
    integer :: ix,iy,iz
    real(dp), dimension(n,n,n) :: rks,rksi,wkf
    real(dp), dimension(n) :: rk,wk
    real(dp) :: fnw
    real(dp), parameter :: fac = 36.d0
   
    fnw = 2.d0/dble(n)

    rk(1) = 0.d0
!$omp parallel do private(ix) shared(rk)
    do ix=2,n/2
       rk(ix) = dble(ix-1)
       rk(n+2-ix) = rk(ix)
    enddo 
!$omp end parallel  do
    rk(n/2+1) = dble(n/2)

!$omp parallel do private(ix) shared(wk,rk)
    do ix=1,n
       wk(ix) = exp(-fac*(rk(ix)*fnw)**fac)
    enddo 
!$omp end parallel  do


!$omp parallel do private(ix,iy,iz) shared(rks,wkf,rk)
    do iz=1,n
       do iy=1,n
          do ix=1,n
             rks(ix,iy,iz) = rk(ix)**2 + rk(iy)**2 + rk(iz)**2
             wkf(ix,iy,iz) = wk(ix)*wk(iy)*wk(iz)
          enddo
        enddo
    enddo
!$omp end parallel do

!$omp parallel do private(ix,iy,iz) shared(rks,rksi)
    do iz=1,n
       do iy=1,n
          do ix=1,n
             if ( (ix==1) .and. (iy==1) .and. (iz==1) ) then 
                rksi(ix,iy,iz) = 0.d0
             else
                rksi(ix,iy,iz) = 1.d0/rks(ix,iy,iz)
             endif
          enddo
        enddo
     enddo
!$omp end parallel do
    return

  end subroutine initarr

! ========= X - DERIV ============ !

  subroutine derivx(n,var,der,rk)

    integer, parameter :: dp = kind(1.d0) 
    integer :: n
    integer :: nw,nw1,np2,ic 
    integer :: ix,iy,iz
    real(dp), dimension(n,n,n) :: var,der
    real(dp), dimension(n) :: rk 
 
    nw = n/2
    nw1  = nw + 1
    np2 = n+2

!$omp parallel do  private(ix,iy,iz) shared(der)  
    do iz = 1,n
       do iy = 1,n
          der(1,iy,iz) = 0.d0
          der(nw1,iy,iz) = 0.d0
       enddo   
    enddo
!$omp end parallel do

!$omp parallel do private(ix,iy,iz) shared(der,var,rk)   
    do iz = 1,n
       do iy = 1,n
          do ix = 2,n-nw
             ic = np2 - ix
             der(ix,iy,iz) = -rk(ix)*var(ic,iy,iz)
             der(ic,iy,iz) =  rk(ix)*var(ix,iy,iz)
          enddo 
       enddo
     enddo   
!$omp end parallel do    

   return
   end subroutine derivx

! ========= Y - DERIV ============ !

   subroutine derivy(n,var,der,rk)

    integer, parameter :: dp = kind(1.d0) 
    integer :: n
    integer :: nw,nw1,np2,ic 
    integer :: ix,iy,iz
    real(dp), dimension(n,n,n) :: var,der
    real(dp), dimension(n) :: rk 
 
    nw = n/2
    nw1 = nw + 1
    np2 = n+2
    
!$omp parallel do private(ix,iy,iz) shared(der,var,rk)   
    do iz = 1,n
       do ix = 1,n
          der(ix,1,iz) = 0.d0
          der(ix,nw1,iz) = 0.d0
       enddo   
    enddo
!$omp end parallel do

!$omp parallel do
    do iz = 1,n
       do ix = 1,n
          do iy = 2,n-nw
             ic = np2 - iy
             der(ix,iy,iz) = -rk(iy)*var(ix,ic,iz)
             der(ix,ic,iz) =  rk(iy)*var(ix,iy,iz)
          enddo 
       enddo
     enddo   
!$omp end parallel do    
  return
  end subroutine derivy

! ========= Z - DERIV ============ !

   subroutine derivz(n,var,der,rk)

    integer, parameter :: dp = kind(1.d0) 
    integer :: n
    integer :: nw,nw1,np2,ic 
    integer :: ix,iy,iz
    real(dp), dimension(n,n,n) :: var,der
    real(dp), dimension(n) :: rk 
 
    nw = n/2
    nw1 = nw + 1
    np2 = n+2
    
!$omp parallel do private(ix,iy,iz) shared(der,var,rk)   
    do iy = 1,n
       do ix = 1,n
          der(ix,iy,1) = 0.d0
          der(ix,iy,nw1) = 0.d0
       enddo   
    enddo
!$omp end parallel do

!$omp parallel do private(ix,iy,iz) shared(der,var,rk)   
    do iy = 1,n
       do ix = 1,n
          do iz = 2,n-nw
             ic = np2 - iz
             der(ix,iy,iz) = -rk(iz)*var(ix,iy,ic)
             der(ix,iy,ic) =  rk(iz)*var(ix,iy,iz)
          enddo 
       enddo
     enddo   
!$omp end parallel do    
  return

  end subroutine derivz

end module spectral
