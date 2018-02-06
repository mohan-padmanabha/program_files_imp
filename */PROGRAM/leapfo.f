

       subroutine leapfo(n1,n2,n3,nh,vx,vy,vz,
     &                  fx,fy,fz,eev1,eev2,eev3,dt)

**************************************************************
c
c leap-frog centered at t+dt/2.
c
c (v(n+1) - v(n-1))/ dt = f( v(n+1/2) )
c
**************************************************************


        implicitnone

        integer n1,n2,n3,nh

        complex vx(0:n1/2 ,0:n2,0:n3)
        complex vy(0:n1/2 ,0:n2,0:n3)
        complex vz(0:n1/2 ,0:n2,0:n3)

        complex fx(0:n1/2 ,0:n2,0:n3)
        complex fy(0:n1/2 ,0:n2,0:n3)
        complex fz(0:n1/2 ,0:n2,0:n3)
 
        real  eev1(0:nh),eev2(0:nh),eev3(0:nh)
  
        real dt

        real av,expv,ev2

        integer ix,iy,iz


!$OMP  PARALLEL PRIVATE(expv,ev2,av)
!$OMP DO SCHEDULE(RUNTIME)

        do iz=0,n3-1
        do iy=0,n2-1
        do ix=0,n1/2-1

          expv  = eev1(ix)*eev2(iy)*eev3(iz)


          ev2 = expv*expv
          av  = expv*dt

          vx(ix,iy,iz)=vx(ix,iy,iz)*ev2+fx(ix,iy,iz)*av
          vy(ix,iy,iz)=vy(ix,iy,iz)*ev2+fy(ix,iy,iz)*av
          vz(ix,iy,iz)=vz(ix,iy,iz)*ev2+fz(ix,iy,iz)*av

        enddo
        enddo
        enddo

!$OMP END DO
!$OMP END PARALLEL 

        return
        end
