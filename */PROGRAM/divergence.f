      subroutine divergence(n1,n2,n3,nh,vox,voy,voz,vx,vy,vz,
     &                 kx,ky,kz,k2x,k2y,k2z)

      integer n1,n2,n3


      complex vx(0:n1/2,0:n2,0:n3)
      complex vy(0:n1/2,0:n2,0:n3)
      complex vz(0:n1/2,0:n2,0:n3)

      complex vox(0:n1/2,0:n2,0:n3)
      complex voy(0:n1/2,0:n2,0:n3)
      complex voz(0:n1/2,0:n2,0:n3)

      real   kx(0:nh),ky(0:nh),kz(0:nh)
      real   k2x(0:nh),k2y(0:nh),k2z(0:nh)

      real div,div1,ks,epsilon
      integer ix,iy,iz

      complex s
 
      div1 = 0.
      epsilon = 1.E-22


      do iz=0,n3-1
      do iy=0,n2-1
      do ix=0,n1/2-1

        div = abs( kx(ix)*vx(ix,iy,iz)
     &            +ky(iy)*vy(ix,iy,iz)
     &            +kz(iz)*vz(ix,iy,iz) )

        div1=amax1(div,div1)

      enddo
      enddo
      enddo

      write(*,2) div1

      if (div1.ge.1.E-09) then

!$OMP PARALLEL DO SCHEDULE(AUTO)
 
        do l=0,n3-1
        do j=0,n2-1
        do i=0,n1/2-1

          ks=k2x(i)+k2y(j)+k2z(l)+epsilon

          s=   (  kx(i)*vx(i,j,l)
     &          + ky(j)*vy(i,j,l) 
     &          + kz(l)*vz(i,j,l) )/ks

          vx(i,j,l)=vx(i,j,l)-s*kx(i)
          vy(i,j,l)=vy(i,j,l)-s*ky(j)
          vz(i,j,l)=vz(i,j,l)-s*kz(l)

          s=   (  kx(i)*vox(i,j,l)
     &          + ky(j)*voy(i,j,l) 
     &          + kz(l)*voz(i,j,l) )/ks

          vox(i,j,l)=vox(i,j,l)-s*kx(i)
          voy(i,j,l)=voy(i,j,l)-s*ky(j)
          voz(i,j,l)=voz(i,j,l)-s*kz(l)

        enddo
        enddo
        enddo

!$OMP END PARALLEL DO


        write(*,*)' '
        write(*,*)'======================================'
        write(*,*)'      DIVERGENCE RESET TO ZERO        '
        write(*,*)'======================================'
        write(*,*)' '


      endif

2     format(/' Maximum for div( V ) is :',e18.11,/)

       return
       end
