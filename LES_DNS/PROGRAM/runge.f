
      subroutine runge(n1,n2,n3,n,nh,nw,vox,voy,voz,vx,vy,vz,
     &                 fx,fy,fz,dt,ktrunk,cfl,it,inj_energy,
     &                 kx,ky,kz,k2x,k2y,k2z,cz,eev1,eev2,eev3,
     &                 wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)

      implicitnone

      integer n1,n2,n3,n,nh,nw,ktrunk,it
      real dt,cfl,inj_energy

      real  k2x(0:nh),k2y(0:nh),k2z(0:nh)
      real kx(0:nh),ky(0:nh),kz(0:nh), cz(0:nh)

      real  eev1(0:nh),eev2(0:nh),eev3(0:nh)

      real vx(0:n1+1,0:n2,0:n3)
      real vy(0:n1+1,0:n2,0:n3)
      real vz(0:n1+1,0:n2,0:n3)

      real vox(0:n1+1,0:n2,0:n3)
      real voy(0:n1+1,0:n2,0:n3)
      real voz(0:n1+1,0:n2,0:n3)

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


c----------------------------------------------------------------
c
c   Step 1a (Euler) :  (V(t+dt/2)-V(t))/(dt/2)  =  F(t)
c
c input:   vo = V(k,t)          output : vo = V(k,t+dt/2) 
c          f  = F(k,t)                   f  = F(k,t)   
c         [v  = V(k,t)]                 [v  = V(k,t)]   
c
c----------------------------------------------------------------

      call eulero(n1,n2,n3,nh,vox,voy,voz,fx,fy,fz,eev1,eev2,eev3,dt)

*-----------------------------------------------------------------
c
c   Step 1b : computation of  F(t+dt/2)
c
c input:   vo = V(k,t+dt/2)       output : vo = dummy      
c          f  = F(k,t)                     f  = F(k,t+dt/2) 
c         [v  = V(k,t)]                   [v  = V(k,t)]   
c
*-----------------------------------------------------------------

      call nonlidns(n1,n2,n3,n,nh,vox,voy,voz,fx,fy,fz,dt,ktrunk,
     &    k2x,k2y,k2z,kx,ky,kz,cfl,it,
     &    wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)


*------------------------------------------------------------------------
c
c   Step 2a (leap-frog) :  ( V(t+dt)-V(t) )/(dt)  =  F(t+dt/2) 
c
c  input: v  = V(k,t)              output :  v  = V(k,t+dt) 
c        [vo = dummy]                       [vo = dummy]
c        f  = F(k,t+dt/2)                    f  = F(k,t+dt/2)
c
*------------------------------------------------------------------------

       call leapfo(n1,n2,n3,nh,vx,vy,vz,fx,fy,fz,eev1,eev2,eev3,dt)

*------------------------------------------------------------------------
c
c   Step 2b : put  V(k,t+dt) into  vo
c
c input:  v  = V(k,t+dt)            output : v  = V(k,t+dt) 
c        [f  = F(k,t+dt/2)]                 [f  = F(k,t+dt/2)]
c        [vo = dummy]                       [vo = V(k,t+dt)]
*------------------------------------------------------------------------

      call dcopy(nw,vx,1,vox,1)
      call dcopy(nw,vy,1,voy,1)
      call dcopy(nw,vz,1,voz,1)

*-------------------------------------------------------------------------
c
c   Step 2c  :  computation of  F(t+dt) for next time step
c
c input: vo = V(k,t+dt)                 output : vo = dummy
c        f  = F(k,t+dt/2)                        f  = F(k,t+dt) 
c       [v  = V(k,k+dt)]                        [v  = V(k,t+dt)]           
c
*-------------------------------------------------------------------------

       call nonlidns(n1,n2,n3,n,nh,vox,voy,voz,fx,fy,fz,dt,ktrunk,
     &    k2x,k2y,k2z,kx,ky,kz,cfl,it,
     &    wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)

*-------------------------------------------------------------------------
c
c   Step 3a  :  force  V(k,t+dt)  
c
c input: v  = V(k,t+dt)                 Output : v  = V(k,t+dt) + forcing        
c       [f  = F(k,t+dt)]                        [f  = F(k,t+dt)]
c       [vo = dummy]                            [vo = dummy]
c
*-------------------------------------------------------------------------

       call force_velocity(n1,n2,n3,nh,vx,vy,vz,
     &               dt,cz,k2x,k2y,k2z,inj_energy)

*-------------------------------------------------------------------------
c
c   Step 3b  :  put  V(k,t+dt)+forcing  into v
c
c input: vo = dummy                     output : vo = V(k,t+dt) + forcing
c       [v  = V(k,t+dt) + forcing]              [v  = V(k,t+dt) + forcing]        
c       [f  = F(k,t+dt)]                        [f  = F(k,t+dt)] 
c
*-------------------------------------------------------------------------

      call dcopy(nw,vx,1,vox,1)
      call dcopy(nw,vy,1,voy,1)
      call dcopy(nw,vz,1,voz,1)

       return
       end
