

       subroutine project(n1,n2,n3,nh,fx,fy,fz,
     &                   kx,ky,kz,k2x,k2y,k2z)

***********************************************************************
c   calcul de la projection : Pij( V x W)(k)   dans l'espace spectral.
c**********************************************************************

       implicitnone


       integer n1,n2,n3,nh

       complex fx(0:n1/2,0:n2,0:n3)
       complex fy(0:n1/2,0:n2,0:n3)
       complex fz(0:n1/2,0:n2,0:n3)

       real k2x(0:nh),k2y(0:nh),k2z(0:nh)
       real kx(0:nh),ky(0:nh),kz(0:nh)

       real    xik2,epsilon
       complex txc,tyc,tzc
       integer ix,iy,iz

      epsilon = 1.E-22

!$OMP  PARALLEL PRIVATE(xik2,txc,tyc,tzc)
!$OMP DO SCHEDULE(RUNTIME)


      do iz=0,n3-1
      do iy=0,n2-1
      do ix=0,n1/2-1

       xik2  = 1./(k2x(ix)+k2y(iy)+k2z(iz)+epsilon)
       
       txc           = (1.-k2x(ix)*xik2)*fx(ix,iy,iz)
     &               - ( kx(ix)*ky(iy)*fy(ix,iy,iz)
     &               +   kx(ix)*kz(iz)*fz(ix,iy,iz) )*xik2

       tyc           = (1.-k2y(iy)*xik2)*fy(ix,iy,iz)
     &               - ( ky(iy)*kx(ix)*fx(ix,iy,iz)
     &               +   ky(iy)*kz(iz)*fz(ix,iy,iz) )*xik2

       tzc           = (1.-k2z(iz)*xik2)*fz(ix,iy,iz)
     &               - ( kz(iz)*kx(ix)*fx(ix,iy,iz)
     &               +   kz(iz)*ky(iy)*fy(ix,iy,iz) )*xik2

       fx(ix,iy,iz)=txc
       fy(ix,iy,iz)=tyc
       fz(ix,iy,iz)=tzc

       enddo
       enddo
       enddo

!$OMP END DO
!$OMP END PARALLEL 

       return
       end
