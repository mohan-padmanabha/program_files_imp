

c=======================================================================
c
c                    Program to simulate
c
c                  Navier Stokes equations
c
c
c                   *** December 2010 ***
c
c=======================================================================
c
c   d U_i / dt - Pij (U x W })_j =    nu*d_jj U_i
c
c Where:
c         U                                 is the velocity
c         Pij = \delta_ij - k_i k_j / k^2   is the projector 
c
************************************************************************

      program navisto

      implicit none

      include "netcdf.inc"

*-------------------------------------------------------
*        resolution of the DNS
*-------------------------------------------------------

      integer n1,n2,n3,n,nh,nw

      parameter( n1  = 128 )    ! Resolution in X direction
      parameter( n2  = 128 )    ! Resolution in Y direction
      parameter( n3  = 128 )    ! Resolution in Z direction
      parameter( n   = 128 )    ! max(n1,n2,n3)
      parameter( nh  = 128 )    ! Size of spectra to be computed
                               ! nh >= max[n1,n2,n3] all the lines has to be same
      integer ktrunk
      parameter(ktrunk = 42)    ! 2/3*(n/2) for alayasing remove truncation error

*-------------------------------------------------------

* --- time parameters and  statistics

      integer itmax,it,nstat,nfield 
      real    time,dt,cfl

      real   eners(0:nh) ! store enery spectra
      real   enert(0:nh) !
      real   enerf(0:nh)

      real reynolds,reynolds_taylor ! values which are computed
      real integral_scale,taylor_scale,kol_scale
      real vrms,energy,enstrophy,dissipation
      real turnover_time,inj_energy

* --- NetCDF ID

      integer nsid,iddims,idtims,idtime,idreyn,idreyn_tay,idinteg_sc ! values for netcdf
      integer idtayl_sc,idkolm_sc,idrms_vel,iddissip,idinj,idturn,idener
      integer idener_sp,idener_tr,idener_fl,status


* --- dissipation parameters

      real    nu,nuf

* --- wavenumbers

      real   kx(0:nh),ky(0:nh),kz(0:nh) ! wave number
      real   k2x(0:nh),k2y(0:nh),k2z(0:nh),cz(0:nh)
      real   eev1(0:nh),eev2(0:nh),eev3(0:nh)
 
* --- FFT arrays

!------------------- temperton ----------------------------
      real wk(4*(n+2)**2,0:1) ! for fft convertion !! 
      complex ex1(0:n),ex2(0:n),ex3(0:n) ! to save computation form fft
      integer ifax1(100),ifax2(100),ifax3(100)
!------------------- ESSL ---------------------------------
!      real   wk(62000+3*n1+(2*n3+257)*((n1/2+1)*n2+5))
!      integer ex1,ex2,ex3,ifax1,ifax2,ifax3
!---------------------------------------------------------

* --- 3D arrays

      real vx(0:n1+1,0:n2,0:n3) ! 3d variables and size is big because fft use until n-1
      real vy(0:n1+1,0:n2,0:n3)
      real vz(0:n1+1,0:n2,0:n3)

      real vox(0:n1+1,0:n2,0:n3) 
      real voy(0:n1+1,0:n2,0:n3)
      real voz(0:n1+1,0:n2,0:n3)

      real fx(0:n1+1,0:n2,0:n3)
      real fy(0:n1+1,0:n2,0:n3)
      real fz(0:n1+1,0:n2,0:n3)

      integer  i,nb_procs,OMP_GET_NUM_THREADS
      real t1,t2

      INTEGER nb_periodes_initial
      INTEGER nb_periodes_final
      INTEGER nb_periodes_max
      INTEGER nb_periodes_sec
      INTEGER nb_periodes     
      REAL  temps_elapsed  


*=================================================================
*                        Initialisations.
*=================================================================

      nb_procs = 1

!$OMP PARALLEL
!$      nb_procs = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL


      write(*,*)'---------------------------'
      write(*,*)'Number of processors :',nb_procs
      write(*,*)'---------------------------'

*-------------------------------------------------------
*        parameters of the DNS
*-------------------------------------------------------

      nuf   = 1.E-29           ! viscosity
      dt    = 1.000E-03        ! time step  (adapt to cfl)      
      itmax = 100              ! number of time step   

*-------------------------------------------------------
*       initialization of all constant and parameters
*-------------------------------------------------------

      time       = 0.
      it         = 0
      nstat      = 0
      nfield     = 0
      inj_energy = 0.

      nw = (n1+2)*(n2+1)*(n3+1)
      nu = nuf


*-------------------------------------------------------
*          initialisation of wavenumbers
*-------------------------------------------------------

      call defk(n1,n2,n3,nh,kx,ky,kz,k2x,k2y,k2z,cz)
 
      call precal(n1,n2,n3,nh,nu,dt,k2x,k2y,k2z,eev1,eev2,eev3)

*-------------------------------------------------------
*          initialization of random number generator
*-------------------------------------------------------

      call randini

*-------------------------------------------------------
*                initialization of  FFTs
*-------------------------------------------------------


      call initfft(n1,n2,n3,n,
     &  wk,ex1,ex2,ex3,ifax1,ifax2,ifax3 ) ! init of fft

*-------------------------------------------------------
*              initialization of arrays
*-------------------------------------------------------

       call initial(vox,nw,0.)
       call initial(voy,nw,0.)
       call initial(voz,nw,0.)

*-------------------------------------------------------
*   get velocity initial conditions 
*-------------------------------------------------------

      call rdmfld(n1,n2,n3,nh,vox,voy,voz,
     &   kx,ky,kz,k2x,k2y,k2z) !generate randam field in a ginven spectra

      call divergence(n1,n2,n3,nh,vx,vy,vz,vox,voy,voz,
     &                kx,ky,kz,k2x,k2y,k2z)

      call trunck(n1,n2,n3,nh,vox,voy,voz,k2x,k2y,k2z,ktrunk)


      call dcopy(nw,vox,1,vx,1)
      call dcopy(nw,voy,1,vy,1)
      call dcopy(nw,voz,1,vz,1)

*-------------------------------------------------------
*            computation of non-linear term
*-------------------------------------------------------


      call nonlidns(n1,n2,n3,n,nh,vox,voy,voz,fx,fy,fz,dt,ktrunk,
     &      k2x,k2y,k2z,kx,ky,kz,cfl,it,
     &      wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)

      call dcopy(nw,vx,1,vox,1)
      call dcopy(nw,vy,1,voy,1)
      call dcopy(nw,vz,1,voz,1)

*=================================================================
*                         time stepping
*=================================================================

 
      call CPU_TIME(t1)
      CALL SYSTEM_CLOCK(COUNT_RATE=nb_periodes_sec, 
     &                   COUNT_MAX=nb_periodes_max)
      CALL SYSTEM_CLOCK(COUNT=nb_periodes_initial)

1     continue


      write(*,102)'time =',time,'it =',it,
     &            'CFL adv =',cfl

*-------------------------------------------------------
*                     check divergence
*-------------------------------------------------------

      if (mod(it,20) .eq. 0) then ! operation is perfomed every 20 steps
        call divergence(n1,n2,n3,nh,vox,voy,voz,vx,vy,vz,
     &                  kx,ky,kz,k2x,k2y,k2z)
      endif

*-------------------------------------------------------
*                     save field (in NetCDF)
*-------------------------------------------------------

c$$$      if (mod(it,100) .eq. 0) then
c$$$        call ifft(vx,n1,n2,n3,n,wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)
c$$$        call ifft(vy,n1,n2,n3,n,wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)
c$$$        call ifft(vz,n1,n2,n3,n,wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)
c$$$        call write_field(n1,n2,n3,nfield,time,vx,vy,vz)
c$$$        call dcopy(nw,vox,1,vx,1)
c$$$        call dcopy(nw,voy,1,vy,1)
c$$$        call dcopy(nw,voz,1,vz,1)
c$$$        nfield = nfield + 1
c$$$      endif

*-------------------------------------------------------
*                     compute statistics
*-------------------------------------------------------

      if (mod(it,20) .eq. 0) then

        call statistics(n1,n2,n3,nh,vx,vy,vz,fx,fy,fz,
     &               k2x,k2y,k2z,cz,nu,ktrunk,
     &               eners,enert,enerf,
     &               reynolds,reynolds_taylor,
     &               integral_scale,taylor_scale,kol_scale,
     &               vrms,energy,inj_energy,dissipation,
     &               turnover_time,time) !computed for calculation of statistics

*-------------------------------------------------------
*                     save statistics (in NetCDF)
*-------------------------------------------------------

c$$$        call write_stat(nh,eners,enert,enerf,
c$$$     &   reynolds,reynolds_taylor,
c$$$     &   integral_scale,taylor_scale,kol_scale,
c$$$     &   vrms,energy,inj_energy,dissipation,
c$$$     &   turnover_time,time,nstat,
c$$$     &   nsid,iddims,idtims,idtime,idreyn,idreyn_tay,idinteg_sc,
c$$$     &   idtayl_sc,idkolm_sc,idrms_vel,iddissip,idinj,idturn,idener,
c$$$     &   idener_sp,idener_tr,idener_fl)

      endif

*-------------------------------------------------------
*                     Runge-Kutta integration
*-------------------------------------------------------

      call runge(n1,n2,n3,n,nh,nw,vox,voy,voz,vx,vy,vz,
     &           fx,fy,fz,dt,ktrunk,cfl,it,inj_energy,
     &           kx,ky,kz,k2x,k2y,k2z,cz,eev1,eev2,eev3,
     &           wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)

*-------------------------------------------------------
*                     time incrementation
*-------------------------------------------------------

      it    = it + 1
      time  = dt + time

*-------------------------------------------------------
*                     adapting time step
*-------------------------------------------------------

!      dt = dt*0.5/cfl

c$$$      if (time .le. 0.0) then
c$$$          nu = nuf*(1 + 10.D0*exp(-70.D0*time))
c$$$      else
c$$$          nu = nuf
c$$$      endif

!      call precal(n1,n2,n3,nh,nu,dt,k2x,k2y,k2z,eev1,eev2,eev3)

      if (it.le.itmax) goto 1

*=================================================================
*                    end  time stepping
*=================================================================

!      status = nf_close(nsid)

 2    continue

      call CPU_TIME(t2)
      CALL SYSTEM_CLOCK(COUNT=nb_periodes_final)
      nb_periodes = nb_periodes_final - nb_periodes_initial
      IF (nb_periodes_final < nb_periodes_initial) THEN
        nb_periodes = nb_periodes + nb_periodes_max
      ENDIF
      temps_elapsed   = REAL(nb_periodes) / nb_periodes_sec
      write(*,*)'ELAPSED TIME :',temps_elapsed
      write(*,*)'CPU TIME :',t2-t1

      write(*,*)' '
      write(*,*)'END PROGRAM'
      write(*,*)' '

      stop

102   format(A6,e13.5,3x,A4,i7,3x,A9,e11.5)

      end
