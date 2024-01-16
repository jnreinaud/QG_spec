module initq
  implicit none

  integer, parameter :: ddp=kind(1.d0)

contains

  real(ddp) function sinc(x)
    real(ddp), parameter :: pi = 4.d0*atan(1.d0)
    real(ddp), parameter :: small = 1.d-12
    real(ddp) :: x 
    if (x < small) then
       sinc = 1.d0 
    else
       sinc = sin(x)/x
    endif
    return
  end function sinc

  subroutine init(n,q)
    real(ddp), parameter :: pi = 4.d0*atan(1.d0)
    real(ddp), parameter :: g1 = 5.763459196895054d0
    real(ddp), parameter :: eps = 1.02d0
    integer :: n
    integer :: ix,iy,iz

    real(ddp) :: x,y,z,r,rp,th
    real(ddp) :: dg    

    real(ddp), dimension(n,n,n) :: q

    dg = 2.d0*pi/dble(n)

    do iz=1,n
       z = - pi + dble(iz)*dg
       r = z**2
       do iy=1,n
          y = - pi + dble(iy)*dg
          r = y**2 + r
          do ix=1,n
             x = - pi + dble(ix)*dg
             r = x**2 + r
             th = atan2(y,x)
             rp  = (x/eps)**2+y**2 + z**2
             if ( rp <= 1.d0) then
                q(ix,iy,iz) = -(sinc(g1*sqrt(rp))-sinc(g1))/sinc(g1)
             else
                q(ix,iy,iz) = 0.d0
             endif
       enddo
    enddo
  enddo
  
  return

  end subroutine init

end module initq
