

      subroutine rdmfld(n1,n2,n3,nh,vx,vy,vz,
     &   kx,ky,kz,k2x,k2y,k2z)

***************************************************************************
c
c  computes a random isotropic of vortical 3d field starting with
c  given spectra
c     
***************************************************************************

      implicitnone

      integer n1,n2,n3,nh

      complex vx(0:n1/2,0:n2,0:n3)
      complex vy(0:n1/2,0:n2,0:n3)
      complex vz(0:n1/2,0:n2,0:n3)

      real   kx(0:nh),ky(0:nh),kz(0:nh)
      real   k2x(0:nh),k2y(0:nh),k2z(0:nh)

      real    v11,v12,v21,v22,v13,v23,ro1,ro2,ro3 
      real    theta1,theta2,theta3
      real    pid,epsi,epsilon,ks,ko

      complex s

      integer i,j,l,k

      real   vect1(0:nh),vect2(0:nh)

c initial constants

      pid     = 2.*acos(-1.0)
      epsi    = 1.e-06
      epsilon = 1.e-22
      ko      = 2.0

c tableaux a zero.

      vx(:,:,:) = 0.
      vy(:,:,:) = 0.
      vz(:,:,:) = 0.

      vect1(:) = 0.
      vect2(:) = 0.

c--------------------------------------------------------------------
c       choose the energy spectra of the random velocity field
c--------------------------------------------------------------------

      vect1(0) = 0.
      do k=1,nh
         vect1(k)=0.06*exp(-(float(k)-ko)**2)
      enddo

c--------------------------------------------------------------------
c       compute the number of wavenumbers per spectrum shell
c--------------------------------------------------------------------

c for x=0

      do l=0,n3-1
      do j=0,n2-1
        k=int(sqrt(k2z(l)+k2y(j))+0.5+epsilon)
        if (k.le.nh) vect2(k)=vect2(k)+1.
      enddo
      enddo

c for x=[1,kmax]

      do l=0,n3-1
      do j=0,n2-1
      do i=1,n1/2
        k=int(sqrt(k2x(i)+k2y(j)+k2z(l))+0.5+epsilon)
        if (k.le.nh) vect2(k)=vect2(k)+2.
      enddo
      enddo
      enddo

      do k=1,nh
        if (vect2(k).ne.0.) vect2(k)=sqrt(vect1(k)/vect2(k))
      enddo

c--------------------------------------------------------------------
c       compute random velocity field
c--------------------------------------------------------------------



      do l=0,n3-1
      do j=0,n2-1
      do i=0,n1/2

        k=int(sqrt(k2x(i)+k2y(j)+k2z(l)+epsilon)+0.5)

        if ((k.lt.nh) .and. (vect2(k) .gt. epsi))   then

          call  randa(v11)
          call  randa(v21)
          call  randa(v12)
          call  randa(v22)

          call  randa(v13)
          call  randa(v23)
          theta1 = v11*pid
          theta2 = v12*pid

          theta3 = v13*pid
          v21 = min(v21,1.-epsi)
          v21 = max(v21,epsilon)
          v22 = min(v22,1.-epsi)
          v22 = max(v22,epsilon)
          v23 = min(v23,1.-epsi)
          v23 = max(v23,epsilon)

          ro1    = sqrt(-alog(1.-v21))
          ro2    = sqrt(-alog(1.-v22))
          ro3    = sqrt(-alog(1.-v23))

          v11    = ro1*cos(theta1)
          v21    = ro1*sin(theta1)
          v12    = ro2*cos(theta2)
          v22    = ro2*sin(theta2)
          v13    = ro3*cos(theta3)
          v23    = ro3*sin(theta3)

          vx(i,j,l)=vect2(k)*cmplx(v11,v21)
          vy(i,j,l)=vect2(k)*cmplx(v12,v22)
          vz(i,j,l)=vect2(k)*cmplx(v13,v23)

       else
         vx(i,j,l)= cmplx(0.,0.)
         vx(i,j,l)= cmplx(0.,0.)
         vy(i,j,l)= cmplx(0.,0.)
       endif


      enddo
      enddo
      enddo

 10   format(3(I3,1x),6(E10.3,2x))


c--------------------------------------------------------------------
c imposing the conjugate hermitian condition et solinoidal 
c--------------------------------------------------------------------


      do i=1,n2/2-1
      do j=1,n3/2-1
        vx(0,i,j)=cmplx(real(vx(0,n2-i,n3-j)),-imag(vx(0,n2-i,n3-j)))
        vy(0,i,j)=cmplx(real(vy(0,n2-i,n3-j)),-imag(vy(0,n2-i,n3-j)))
        vz(0,i,j)=cmplx(real(vz(0,n2-i,n3-j)),-imag(vz(0,n2-i,n3-j)))
      enddo
      enddo

      do j=n3/2+1,n3-1
      do i=1,n2/2-1
          vx(0,i,j)=cmplx(real(vx(0,n2-i,n3-j)),-imag(vx(0,n2-i,n3-j)))
          vy(0,i,j)=cmplx(real(vy(0,n2-i,n3-j)),-imag(vy(0,n2-i,n3-j)))
          vz(0,i,j)=cmplx(real(vz(0,n2-i,n3-j)),-imag(vz(0,n2-i,n3-j)))
      enddo
      enddo

      do j=1,n3/2-1
          vx(0,0,j)=cmplx(real(vx(0,0,n3-j)),-imag(vx(0,0,n3-j)))
          vy(0,0,j)=cmplx(real(vy(0,0,n3-j)),-imag(vy(0,0,n3-j)))
          vz(0,0,j)=cmplx(real(vz(0,0,n3-j)),-imag(vz(0,0,n3-j)))
      enddo

      do i=1,n2/2-1
          vx(0,i,0)=cmplx(real(vx(0,n2-i,0)),-imag(vx(0,n2-i,0)))
          vy(0,i,0)=cmplx(real(vy(0,n2-i,0)),-imag(vy(0,n2-i,0)))
          vz(0,i,0)=cmplx(real(vz(0,n2-i,0)),-imag(vz(0,n2-i,0)))
      enddo

      vx(0,0,0)=cmplx(0.,0.)
      vy(0,0,0)=cmplx(0.,0.)
      vz(0,0,0)=cmplx(0.,0.)

c avoiding side effects.

      call zerkam(vx,n1,n2,n3)
      call zerkam(vy,n1,n2,n3)
      call zerkam(vz,n1,n2,n3)

c renders V solenoidal



      do l=0,n3-1
      do j=0,n2-1
      do i=0,n1/2-1

      ks=k2x(i)+k2y(j)+k2z(l)+epsilon

      s=   (  kx(i)*vx(i,j,l)
     &      + ky(j)*vy(i,j,l) 
     &      + kz(l)*vz(i,j,l) )/ks

      vx(i,j,l)=vx(i,j,l)-s*kx(i)
      vy(i,j,l)=vy(i,j,l)-s*ky(j)
      vz(i,j,l)=vz(i,j,l)-s*kz(l)

      enddo
      enddo
      enddo

      return
      end
