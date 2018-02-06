

       subroutine eulero(n1,n2,n3,nh,vx,vy,vz,
     &                  fx,fy,fz,eev1,eev2,eev3,dt)

*******************************************
c
c  Euler at t with dt/2
c
c   ( v(n+1/2)-v(n) ) / (dt/2)  =  f( v(n) )
c
*******************************************

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

        integer ix,iy,iz 
        real    expv,dth

        dth = dt/2.

!$OMP  PARALLEL PRIVATE(expv)
!$OMP DO SCHEDULE(RUNTIME)

        do iz=0,n3-1
        do iy=0,n2-1
        do ix=0,n1/2-1

        
        expv  = eev1(ix)*eev2(iy)*eev3(iz)


        vx(ix,iy,iz)=(vx(ix,iy,iz)+fx(ix,iy,iz)*dth)*expv
        vy(ix,iy,iz)=(vy(ix,iy,iz)+fy(ix,iy,iz)*dth)*expv
        vz(ix,iy,iz)=(vz(ix,iy,iz)+fz(ix,iy,iz)*dth)*expv


        enddo
        enddo
        enddo

!$OMP END DO
!$OMP END PARALLEL 


        return
        end
