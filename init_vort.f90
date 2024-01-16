module initq
  implicit none

  contains

  subroutine init(n,q)

  integer, parameter :: dp = kind(1.d0) 
  real(dp), parameter :: pi = 4.d0*atan(1.d0)
  
  integer :: n
  integer :: ix,iy,iz

  integer :: nv,iv

  real(dp) :: x,y,z
  real(dp) :: xv,yv,zv,rv,qv
  real(dp) :: fac

  real(dp) :: r

  real(dp) :: dg    
   
  real(dp), dimension(n,n,n) :: q
  real(dp) :: dx,dz
  
  write(*,*) ' Enter the number vortices'
  read(*,*) nv 
  write(*,*) ' Enter fac q ~ exp(-r**fac) '
  read(*,*) fac
      
  dg = 2.d0*pi/dble(n)

  do iz=1,n
     do iy=1,n
        do ix=1,n
           q(ix,iy,iz) = 0.d0
        enddo
     enddo   
  enddo
  
  do iv=1,nv
     write(*,*) '  - Enter the center x,y,z, the radius r, and the pv of vortex ',iv
     read(*,*) xv,yv,zv,rv,qv
     do iz=1,n
        z = - pi + dble(iz)*dg
        do iy=1,n
           y = - pi + dble(iy)*dg
           do ix=1,n
              x = - pi + dble(ix)*dg
              r = sqrt((x-xv)**2 + (y-yv)**2 + (z-zv)**2)/rv
              q(ix,iy,iz) = q(ix,iy,iz) + qv*exp(-r**fac)
           enddo
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
