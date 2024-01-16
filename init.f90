module initq
  implicit none

  contains

  subroutine init(n,q)

  integer, parameter :: dp = kind(1.d0) 
  real(dp), parameter :: pi = 4.d0*atan(1.d0)
  real(dp), parameter :: fac = 8.d0
  real(dp), parameter :: ref = 0.25d0, cef = 0.25d0*ref 
  
  integer :: n
  integer :: ix,iy,iz

  real(dp) :: x,y,z
  real(dp) :: r1,r2,r3,r4,r5,r6
  real(dp) :: f1,f2,f3,f4,f5,f6

  real(dp) :: dg    
   
  real(dp), dimension(n,n,n) :: q
  real(dp) :: dx,dz
  real(dp), parameter :: qs = 0.3d0, qc = -1.d0

  write(*,*) ' Enter dx/r and dz/r '
  read(*,*) dx,dz
  dx = 0.5d0*ref*dx
  dz = 0.5d0*ref*dz
    
  dg = 2.d0*pi/dble(n)

  do iz=1,n
     do iy=1,n
        do ix=1,n
           q(ix,iy,iz) = 0.d0
        enddo
     enddo   
  enddo
  
  do iz=1,n
    z = - pi + dble(iz)*dg
    do iy=1,n
       y = - pi + dble(iy)*dg
       do ix=1,n
          x = - pi + dble(ix)*dg
          r1 = ((x-dx)**2 + y**2)/ref**2 + ((z-ref-cef-dz)**2)/cef**2
          r2 = ((x-dx)**2 + y**2)/ref**2 + ((z-dz)**2)/ref**2
          r3 = ((x-dx)**2 + y**2)/ref**2 + ((z+ref+cef-dz)**2)/cef**2
          r4 = ((x+dx)**2 + y**2)/ref**2 + ((z-ref-cef+dz)**2)/cef**2
          r5 = ((x+dx)**2 + y**2)/ref**2 + ((z+dz)**2)/ref**2
          r6 = ((x+dx)**2 + y**2)/ref**2 + ((z+ref+cef+dz)**2)/cef**2
          f1 = qs*exp(-fac*r1**fac)
          f2 = qc*exp(-fac*r2**fac)
          f3 = qs*exp(-fac*r3**fac)
          f4 = qs*exp(-fac*r4**fac)
          f5 = qc*exp(-fac*r5**fac)
          f6 = qs*exp(-fac*r6**fac)
          q(ix,iy,iz) = q(ix,iy,iz) + f1+f2+f3+f4+f5+f6
       enddo
    enddo

   enddo
  
  call dumpinit(n,q)

  return

  end subroutine init

  subroutine dumpinit(n,q)

  integer, parameter :: dp = kind(1.d0) 
  real(dp), dimension(n,n,n) :: q
  integer :: n
  integer :: ix,iy,iz
  open(unit=99,file='init_cross.dat',status='replace',action='write')
  do ix=1,n
     do iz=1,n
        write(99,*) q(ix,n/2,iz)
     enddo
  enddo
  close(99)
  return
  end subroutine dumpinit
  
end module initq
