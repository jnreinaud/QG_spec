program qgcs

  use commons
  use initq
  use fft3d
  use stafft
  use spectral

  implicit none

  integer, parameter :: dp = kind(1.d0)
  
  real(dp), parameter :: pi = 4.d0*atan(1.d0)
  real(dp), parameter :: cfl = 0.6d0

  real(dp), dimension(n,n,n) :: q,psi
  real(dp), dimension(n,n,n) :: u,v,qx,qy
  real(dp), dimension(n,n,n) :: qi,qf
  real(dp), dimension(n,n,n) :: rks,rksi,wkf
  real(dp), dimension(n) :: rk,wk

  real(dp), dimension(2*n) :: trig
  integer, dimension(5) :: fac

  real(dp) :: dt,dfac,nu,dg,umax

  real(dp) :: t,tmax,tstep,tsav

  integer :: ix,iy,iz

  integer :: nbytes,ite,nloop

  dg = 2.d0*pi/dble(n) 

  call readin
  call initfft(n,fac,trig)
  call initarr(n,rk,wk,rks,rksi,wkf)
  call init(n,q)

  call ptospc(n,q,trig,fac)

  ite = 0

  call dumpq
  ite = ite + 1

  t = 0.d0
  nloop = 1
  tsav = tstep
  
  do while ( t < tmax )
     call stepping
     write(*,*) t,dt
     t = t + dt
     if (t > tsav) then
        call dumpq
        call dumpz(n/2)
        call specene
        ite = ite + 1
        tsav = tsav + tstep
      endif
  enddo

  stop

  contains

! ==================================

! -----  time stepping -------------
  
    subroutine stepping

      real(dp) :: dt2
      real(dp) :: tmp

      dt2 = 2.d0*dt
      tmp = dt2*nu
!$omp parallel do private(ix,iy,iz) shared(qi,q)
      do iz=1,n
         do iy=1,n
            do ix=1,n
               qf(ix,iy,iz) = q(ix,iy,iz)
            enddo
         enddo   
      enddo
!$omp end parallel do
      call tend(0)

      do iz=1,n
         do iy=1,n
            do ix=1,n
               q(ix,iy,iz) = (qi(ix,iy,iz)+dt2*psi(ix,iy,iz))/(1.d0+tmp*rks(ix,iy,iz))
            enddo
         enddo   
      enddo   
!$omp end parallel do
      do iz=1,n
         do iy=1,n
            do ix=1,n
               qi(ix,iy,iz) = qf(ix,iy,iz)
            enddo
         enddo   
      enddo

      return

    end subroutine stepping

    subroutine euler

      real(dp) :: tmp

      tmp = dt*nu
!$omp parallel do private(ix,iy,iz) shared(qi,q)
      do iz=1,n
         do iy=1,n
            do ix=1,n
               qi(ix,iy,iz) = q(ix,iy,iz)
            enddo
         enddo   
      enddo
!$omp end parallel do
      call tend(0)

      do iz=1,n
         do iy=1,n
            do ix=1,n
               q(ix,iy,iz) = (qi(ix,iy,iz)+dt*psi(ix,iy,iz))/(1.d0+tmp*rks(ix,iy,iz))
            enddo
         enddo   
      enddo   
!$omp end parallel do
      do iz=1,n
         do iy=1,n
            do ix=1,n
               qi(ix,iy,iz) = qf(ix,iy,iz)
            enddo
         enddo   
      enddo

      return

    end subroutine euler


    
! ------ tendencies -------------------

    subroutine tend(iopt)

      integer :: iopt
      real(dp) :: tmp1,tmp2

!$omp parallel do private(ix,iy,iz) shared(q,rksi,psi)
      do iz=1,n
         do iy=1,n
            do ix=1,n
               psi(ix,iy,iz) = -rksi(ix,iy,iz)*q(ix,iy,iz)
            enddo
         enddo
      enddo
!$omp end parallel do

      call derivx(n,q,qx,rk)
      call derivy(n,q,qy,rk)
      call derivx(n,psi,v,rk)
      call derivy(n,psi,u,rk)

!$omp parallel do private(ix,iy,iz) shared(wkf,qx,qy,u,v)
      do iz=1,n
         do iy=1,n
            do ix=1,n
               qx(ix,iy,iz) =  qx(ix,iy,iz)*wkf(ix,iy,iz)
               qy(ix,iy,iz) =  qy(ix,iy,iz)*wkf(ix,iy,iz)
               u(ix,iy,iz)  =  -u(ix,iy,iz)*wkf(ix,iy,iz)
               v(ix,iy,iz)  =   v(ix,iy,iz)*wkf(ix,iy,iz)
            enddo
         enddo
      enddo
!$omp end parallel do

      call spctop(n,qx,trig,fac)
      call spctop(n,u,trig,fac)
      call spctop(n,qy,trig,fac)
      call spctop(n,v,trig,fac)
       
!$omp parallel do private(ix,iy,iz,tmp1,tmp2) shared(psi,qx,qy,u,v)
      do iz=1,n
         do iy=1,n
            do ix=1,n
               tmp1 = u(ix,iy,iz)*qx(ix,iy,iz)
               tmp2 = v(ix,iy,iz)*qy(ix,iy,iz)
               psi(ix,iy,iz) = - (tmp1+tmp2)
            enddo
         enddo
      enddo
!$omp end parallel do

      call ptospc(n,psi,trig,fac)

      if (iopt == 0) then
         umax = 0.d0
         do iz=1,n
            do iy=1,n
               do ix=1,n
                  tmp1 = u(ix,iy,iz)**2 + v(ix,iy,iz)**2
                  if (tmp1 > umax ) then
                     umax = tmp1
                  endif
               enddo
            enddo
          enddo  
          umax = sqrt(umax)
      endif

      return

    end subroutine tend

! ------  init --------------------

    subroutine readin
      open(unit=10,file='param.dat',status='old')
      read(10,*) nu,tmax,tstep
      close(10)
      return
    end subroutine readin

! ------ dunp ---------------------

    subroutine dumpq

      character(len=20) :: fname 
      fname = 'q_t0000.dat'
      write(fname(4:7),'(i4.4)') ite

      open(99,file=fname,status='replace',action='write',form='unformatted')

!$omp parallel do private(ix,iy,iz) shared(psi,q)
      do iz=1,n
      do iy=1,n
      do ix=1,n
         psi(ix,iy,iz) = q(ix,iy,iz)
      enddo
      enddo
      enddo
!$omp end parallel do

      call spctop(n,psi,trig,fac)

      write(99) sngl(t),sngl(psi)
      close(99)
      return
    end subroutine dumpq

! ------ dumpx - crossection x=constant '
    subroutine dumpx(nx)
      integer :: nx
      character(len=20) :: fname 
      fname = 'x_0000_t0000.dat'
      write(fname(3:6),'(i4.4)') nx
      write(fname(9:12),'(i4.4)') ite
      open(99,file=fname,status='replace',action='write') 
      do iy=1,n
         do iz=1,n
            write(99,*) psi(nx,iy,iz)
         enddo
       enddo
       close(99)
       return
     end subroutine dumpx


! ------ dumpz - crossection z=constant '
    subroutine dumpz(nz)
      integer :: nz
      character(len=20) :: fname 
      fname = 'z_0000_t0000.dat'
      write(fname(3:6),'(i4.4)') nz
      write(fname(9:12),'(i4.4)') ite
      open(99,file=fname,status='replace',action='write') 
      do ix=1,n
         do iy=1,n
            write(99,*) psi(ix,iy,nz)
         enddo
       enddo
       close(99)
       return
     end subroutine dumpz

! ------ dumpy - crossection y=constant '
    subroutine dumpy(ny)
      integer :: ny
      character(len=20) :: fname 
      fname = 'y_0000_t0000.dat'
      write(fname(3:6),'(i4.4)') ny
      write(fname(9:12),'(i4.4)') ite
      open(99,file=fname,status='replace',action='write') 
      do ix=1,n
         do iz=1,n
            write(99,*) psi(ix,ny,iz)
         enddo
       enddo
       close(99)
       return
     end subroutine dumpy

! ---------- kinetic energy ---- !
! === CAREFUL recompute psi  === !

     subroutine specene
       character(len=20) :: fname 
       real(dp), dimension(n) :: sene
       integer :: kene
       fname = 'spectrum_0000.dat'
       write(fname(10:13),'(i4.4)') ite
       open(99,file=fname,status='replace',action='write')

      do iz=1,n
         do iy=1,n
            do ix=1,n
               psi(ix,iy,iz) = -rksi(ix,iy,iz)*q(ix,iy,iz)
            enddo
         enddo
      enddo

      call derivx(n,psi,v,rk)
      call derivy(n,psi,u,rk)

       do kene=1,n
          sene(kene) = 0.d0
       enddo 
       do ix=1,n
          do iy=1,n
             do iz=1,n
                kene = sqrt(rk(ix)**2+rk(iy)**2+rk(iz)**2)+1
                sene(kene) = sene(kene) + u(ix,iy,iz)**2 + v(ix,iy,iz)**2
             enddo
          enddo
       enddo
       do ix=1,n/2+1
          write(99,*) sene(ix)
       enddo 
       close(99)
       return
     end subroutine specene  
   
end program qgcs

