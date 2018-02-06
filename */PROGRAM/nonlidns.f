
       subroutine nonlidns(n1,n2,n3,n,nh,vx,vy,vz,fx,fy,fz,dt,ktrunk,
     &      k2x,k2y,k2z,kx,ky,kz,cfl,it,
     &      wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)


c*****************************************************************
c Computation of non-linear term :
c
c           F(k,t) = Pij [V x W](k,t)
c
c*****************************************************************

      implicitnone

      integer n1,n2,n3,n,nh,ktrunk,it
      real dt,cfl

      real kx(0:nh),ky(0:nh),kz(0:nh)
      real k2x(0:nh),k2y(0:nh),k2z(0:nh)

      real vx(0:n1+1,0:n2,0:n3)
      real vy(0:n1+1,0:n2,0:n3)
      real vz(0:n1+1,0:n2,0:n3)

      real fx(0:n1+1,0:n2,0:n3)
      real fy(0:n1+1,0:n2,0:n3)
      real fz(0:n1+1,0:n2,0:n3)

!------------------- temperton ----------------------------
      real wk(4*(n+2)**2,0:1)
      complex ex1(0:n),ex2(0:n),ex3(0:n)
      integer ifax1(100),ifax2(100),ifax3(100)
!------------------- ESSL ---------------------------------
!      real   wk(62000+3*n1+(2*n3+257)*((n1/2+1)*n2+5))
!      integer ex1,ex2,ex3,ifax1,ifax2,ifax3
!---------------------------------------------------------

      real vmax(3),vvmax(0:n3-1,3),vlimit(3),vcfl(3)
      integer iz

c --> f = w(k)

      call rotv(n1,n2,n3,nh,vx,vy,vz,fx,fy,fz,kx,ky,kz)

c --> v =  u(x)


      call ifft(vx,n1,n2,n3,n,wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)
      call ifft(vy,n1,n2,n3,n,wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)
      call ifft(vz,n1,n2,n3,n,wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)

!--------------  test -----------------
!      call energy(vx,vy,vz,n1,n2,n3)
!--------------------------------------
      
c --> check CFL

      if (mod(it,10) .eq. 0) then


!$OMP  PARALLEL SHARED(vvmax)
!$OMP DO SCHEDULE(RUNTIME)

        do iz=0,n3-1

          vvmax(iz,1) = maxval(abs(vx(:,:,iz:iz)))
          vvmax(iz,2) = maxval(abs(vy(:,:,iz:iz)))
          vvmax(iz,3) = maxval(abs(vz(:,:,iz:iz)))

        enddo

!$OMP END DO
!$OMP END PARALLEL 

        vmax(1) = maxval(vvmax(:,1))
        vmax(2) = maxval(vvmax(:,2))
        vmax(3) = maxval(vvmax(:,3))

        vlimit(1) = 2./(float(n1/2)*dt)
        vlimit(2) = 2./(float(n2/2)*dt)
        vlimit(3) = 2./(float(n3/2)*dt)

        vcfl(:) = vmax(:)/vlimit(:)

        cfl = maxval(vcfl)
 
        if (cfl .gt. 1.0) then
          write(*,*)' '
          write(*,*)'=================================='
          write(*,*)'  CFL advection too large :',cfl
          write(*,*)'=================================='
          write(*,*)' '
          stop
        endif
  
      endif



c --> f =  w(x)

      call ifft(fx,n1,n2,n3,n,wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)
      call ifft(fy,n1,n2,n3,n,wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)
      call ifft(fz,n1,n2,n3,n,wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)

c --> f =  (u x w)(x)

      call vecpro(n1,n2,n3,vx,vy,vz,fx,fy,fz)


c --> f =  (u x w)(k)

      call fft(fx,n1,n2,n3,n,wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)
      call fft(fy,n1,n2,n3,n,wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)
      call fft(fz,n1,n2,n3,n,wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)

c --> truncation in Fourier space for k > ktrunc

      call trunck(n1,n2,n3,nh,fx,fy,fz,k2x,k2y,k2z,ktrunk)

c --> projection to remove pressure


      call project(n1,n2,n3,nh,fx,fy,fz,
     &             kx,ky,kz,k2x,k2y,k2z)



102   format(A9,e11.5)

      return
      end
