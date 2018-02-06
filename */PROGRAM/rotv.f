
       subroutine rotv(n1,n2,n3,nh,vx,vy,vz,wx,wy,wz,kx,ky,kz)

************************************************
c calcul du rotationnel du champ v=(vx,vy,vz).
************************************************

       implicitnone


       integer n1,n2,n3,nh

       complex vx(0:n1/2,0:n2,0:n3)
       complex vy(0:n1/2,0:n2,0:n3)
       complex vz(0:n1/2,0:n2,0:n3)

       complex wx(0:n1/2,0:n2,0:n3)
       complex wy(0:n1/2,0:n2,0:n3)
       complex wz(0:n1/2,0:n2,0:n3)

       real   kx(0:nh),ky(0:nh),kz(0:nh)

       integer ix,iy,iz
       complex wxx,wyy,wzz

!$OMP  PARALLEL PRIVATE(wxx,wyy,wzz)
!$OMP DO SCHEDULE(RUNTIME)

        do iz=0,n3-1
        do iy=0,n2-1
        do ix=0,n1/2-1

       wxx=
     &   cmplx(-ky(iy)*imag(vz(ix,iy,iz))+kz(iz)*imag(vy(ix,iy,iz)),
     &         -kz(iz)*real(vy(ix,iy,iz))+ky(iy)*real(vz(ix,iy,iz)))

       wyy =
     &   cmplx(-kz(iz)*imag(vx(ix,iy,iz))+kx(ix)*imag(vz(ix,iy,iz)),
     &         -kx(ix)*real(vz(ix,iy,iz))+kz(iz)*real(vx(ix,iy,iz)))

       wzz =
     &   cmplx(-kx(ix)*imag(vy(ix,iy,iz))+ky(iy)*imag(vx(ix,iy,iz)),
     &         -ky(iy)*real(vx(ix,iy,iz))+kx(ix)*real(vy(ix,iy,iz)))

       wx(ix,iy,iz)= wxx
       wy(ix,iy,iz)= wyy
       wz(ix,iy,iz)= wzz

       enddo
       enddo
       enddo

!$OMP END DO
!$OMP END PARALLEL 

       return
       end
