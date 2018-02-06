

      subroutine force_velocity(n1,n2,n3,nh,vx,vy,vz,
     &               dt,cz,k2x,k2y,k2z,inj_energy)

***************************************************************************
c
c  Compute a stirring force F(k,t) and
c  add this force to the velocity field.
c
c input:   v  = V(k,t)             output :  v  = V(k,t)  + Fo(k,t)
c          f  = ...                          f  = ...
c
***************************************************************************

      implicitnone

      integer n1,n2,n3,nh
      real inj_energy,dt

      complex vx(0:n1/2,0:n2,0:n3)
      complex vy(0:n1/2,0:n2,0:n3)
      complex vz(0:n1/2,0:n2,0:n3)

      real  cz(0:nh)
      real  k2x(0:nh),k2y(0:nh),k2z(0:nh)


      integer ix,iy,iz
      real    coef,fact,inj,ko,ks,val
      real    alpha,beta,gamma
      complex fx,fy,fz

*-------------------------------------------------------
*           parameters of the forcing
*-------------------------------------------------------

      inj = 1.0*dt
      ko  = 2.

*-------------------------------------------------------

      alpha  = 0.
      beta   = 0.
      gamma  = 0.

!$OMP PARALLEL PRIVATE(ks,fact,val)
!$OMP DO REDUCTION(+:alpha,beta,gamma)

      do iz=0,n3-1
      do iy=0,n2-1
      do ix=0,n1/2-1

       ks=k2x(ix)+k2y(iy)+k2z(iz)

       fact = exp(-abs(ks-ko**2))

       val = cz(ix)*(
     &      real(vx(ix,iy,iz))*real(vx(ix,iy,iz)) +
     &      imag(vx(ix,iy,iz))*imag(vx(ix,iy,iz)) +
     &      real(vy(ix,iy,iz))*real(vy(ix,iy,iz)) +
     &      imag(vy(ix,iy,iz))*imag(vy(ix,iy,iz)) +
     &      real(vz(ix,iy,iz))*real(vz(ix,iy,iz)) +
     &      imag(vz(ix,iy,iz))*imag(vz(ix,iy,iz)) )

       alpha  = alpha + val
       beta   = beta  + val*fact
       gamma  = gamma + val*fact**2

      enddo
      enddo
      enddo

!$OMP END PARALLEL 

      coef = (-beta + sqrt(beta**2 + gamma*inj))/gamma

      inj_energy = 0.

!$OMP PARALLEL PRIVATE(ks,fact,fx,fy,fz)
!$OMP DO REDUCTION(+:inj_energy)

      do iz=0,n3-1
      do iy=0,n2-1
      do ix=0,n1/2-1

       ks=k2x(ix)+k2y(iy)+k2z(iz)

       fact = exp(-abs(ks-ko**2))

       fx=vx(ix,iy,iz)*fact*coef
       fy=vy(ix,iy,iz)*fact*coef
       fz=vz(ix,iy,iz)*fact*coef

         inj_energy = inj_energy +  cz(ix)*(
     &      real(fx)*real(fx) +
     &      imag(fx)*imag(fx) +
     &      real(fy)*real(fy) +
     &      imag(fy)*imag(fy) +
     &      real(fz)*real(fz) +
     &      imag(fz)*imag(fz) +
     &  2.*(
     &      real(fx)*real(vx(ix,iy,iz)) +
     &      imag(fx)*imag(vx(ix,iy,iz)) +
     &      real(fy)*real(vy(ix,iy,iz)) +
     &      imag(fy)*imag(vy(ix,iy,iz)) +
     &      real(fz)*real(vz(ix,iy,iz)) +
     &      imag(fz)*imag(vz(ix,iy,iz))  ))

c----> add force to velocity

         vx(ix,iy,iz) = vx(ix,iy,iz) + fx
         vy(ix,iy,iz) = vy(ix,iy,iz) + fy
         vz(ix,iy,iz) = vz(ix,iy,iz) + fz

      enddo
      enddo
      enddo

!$OMP END PARALLEL 


      inj_energy = inj_energy/dt

      return
      end
