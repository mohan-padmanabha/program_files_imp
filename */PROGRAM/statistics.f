
      subroutine statistics(n1,n2,n3,nh,vx,vy,vz,fx,fy,fz,
     &               k2x,k2y,k2z,cz,nu,ktrunk,
     &               eners,enert,enerf,
     &               reynolds,reynolds_taylor,
     &               integral_scale,taylor_scale,kol_scale,
     &               vrms,energy,inj_energy,dissipation,
     &               turnover_time,time)

      integer n1,n2,n3,nh,ktrunk
      real inj_energy

      complex vx(0:n1/2,0:n2,0:n3)
      complex vy(0:n1/2,0:n2,0:n3)
      complex vz(0:n1/2,0:n2,0:n3)

      complex fx(0:n1/2,0:n2,0:n3)
      complex fy(0:n1/2,0:n2,0:n3)
      complex fz(0:n1/2,0:n2,0:n3)

      real   eners(0:nh),ensts(0:nh)
      real   enert(0:nh),enstt(0:nh)
      real   enerf(0:nh),enstf(0:nh)

      real reynolds,reynolds_taylor
      real integral_scale,taylor_scale,kol_scale,rh
      real vrms,energy,enstrophy,dissipation
      real turnover_time

      real  k2x(0:nh),k2y(0:nh),k2z(0:nh),cz(0:nh)
      real  nu

      integer ix,iy,iz,k,l
      real    xk,xk2,xkk2,u1n2,epsilon,pi

      epsilon = 1.E-22
      pi = acos(-1.)

      eners(:)=0.
      ensts(:)=0.
      enert(:)=0.
      enstt(:)=0.
      enerf(:)=0.
      enstf(:)=0.

      energy         = 0.
      enstrophy      = 0.
      integral_scale = 0.
      dissipation    = 0.

!$OMP PARALLEL PRIVATE(xk2,xk,k,u1n2,flux)
!$OMP DO REDUCTION(+:enert,enstt,eners,ensts)
!$OMP& REDUCTION(+:energy,enstrophy,integral_scale,dissipation)

      do iz=0,n3-1
      do iy=0,n2-1
      do ix=0,n1/2-1

        xk2  = k2x(ix)+k2y(iy)+k2z(iz)+epsilon
        xk   = sqrt(xk2)
        k    = int(xk+.5)

        u1n2 = cz(ix)*(
     &      real(vx(ix,iy,iz))*real(vx(ix,iy,iz)) +
     &      imag(vx(ix,iy,iz))*imag(vx(ix,iy,iz)) +
     &      real(vy(ix,iy,iz))*real(vy(ix,iy,iz)) +
     &      imag(vy(ix,iy,iz))*imag(vy(ix,iy,iz)) +
     &      real(vz(ix,iy,iz))*real(vz(ix,iy,iz)) +
     &      imag(vz(ix,iy,iz))*imag(vz(ix,iy,iz))  )


        flux =cz(ix)*(
     &      real(fx(ix,iy,iz))*real(vx(ix,iy,iz)) +
     &      imag(fx(ix,iy,iz))*imag(vx(ix,iy,iz)) +
     &      real(fy(ix,iy,iz))*real(vy(ix,iy,iz)) +
     &      imag(fy(ix,iy,iz))*imag(vy(ix,iy,iz)) +
     &      real(fz(ix,iy,iz))*real(vz(ix,iy,iz)) +
     &      imag(fz(ix,iy,iz))*imag(vz(ix,iy,iz))  )

        enert(k)  = enert(k) + flux
        enstt(k)  = enstt(k) + flux*xk2

        eners(k)  = eners(k)  + u1n2
        ensts(k)  = ensts(k)  + u1n2*xk2

        energy    = energy    + u1n2
        enstrophy = enstrophy + u1n2*xk2

        if (k.ne.0) integral_scale = integral_scale + u1n2/xk

        dissipation  = dissipation +2.*u1n2*nu*xk2

      enddo
      enddo
      enddo

!$OMP END PARALLEL 

      do k=0,nh-1
      do l=0,k
         enerf(k) = enert(l) + enerf(k)
         enstf(k) = enstt(l) + enstf(k)
      enddo
      enddo


c-------------------------------------------------------
c integral scale : 3/8 x L x [Sum_k E(k)/k]/[Sum_k E(k)] 
c-------------------------------------------------------

       integral_scale   = 0.75*pi*integral_scale/energy

c-------------------------------------------------------
c Taylor micro scale : sqrt(5) x sqrt[ energy / enstrophy ]
c-------------------------------------------------------

       taylor_scale=sqrt(5.*energy/enstrophy)

c-------------------------------------------------------
c Vrms:  sqrt[ energy / enstrophy ]
c-------------------------------------------------------

       vrms=sqrt(2.*energy/3.)

c-------------------------------------------------------
c Turnover time:
c-------------------------------------------------------

       turnover_time=integral_scale/vrms

c-------------------------------------------------------
c Reynolds number (based on integral scale)
c-------------------------------------------------------

       reynolds=vrms*integral_scale/nu

c-------------------------------------------------------
c Reynolds number (based on Taylor scale)
c-------------------------------------------------------

       reynolds_taylor=vrms*taylor_scale/nu

c-------------------------------------------------------
c Kolmogorov scale: (viscosity**3/dissipation)**1/4 
c-------------------------------------------------------

      kol_scale=((nu**3)/dissipation)**(.25)

      rh=float(ktrunk)/(pi/kol_scale)


c -------------------------------------------------------------

       print *,' '
       write(*,13)'Viscosity....................:',nu
       write(*,13)'Reynolds Number..............:',reynolds
       write(*,13)'Taylor Reynolds Number.......:',reynolds_taylor
       write(*,13)'integral scale...............:',integral_scale
       write(*,13)'Taylor micro-scale...........:',taylor_scale
       write(*,13)'Kolmogorov scale (1/kd)......:',kol_scale
       write(*,13)'ktrunc/kd....................:',rh
       write(*,13)'rms velocity.................:',vrms
       write(*,13)'total kenetic energy.........:',energy
       write(*,13)'total enstrophy..............:',enstrophy
       write(*,13)'energy dissipation ..........:',dissipation
       write(*,13)'energy injection.............:',inj_energy
       write(*,13)'turnover time................:',turnover_time
       print *,' '

13    format(a30,1x,E15.9)

      return
      end
