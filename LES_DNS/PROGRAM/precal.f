

       subroutine precal(n1,n2,n3,nh,nu,dt,k2x,k2y,k2z,eev1,eev2,eev3)

       implicitnone

       integer n1,n2,n3,nh
       real    dt,nu

       real   k2x(0:nh),k2y(0:nh),k2z(0:nh),cz(0:nh)
       real   eev1(0:nh),eev2(0:nh),eev3(0:nh)

       integer i

 
       do i=0,n1
           eev1(i)=exp(-nu*k2x(i)*0.5*dt)
       enddo
       do i=0,n2
           eev2(i)=exp(-nu*k2y(i)*0.5*dt)
       enddo
       do i=0,n3
           eev3(i)=exp(-nu*k2z(i)*0.5*dt)
       enddo

       eev1(0)  =1.
       eev1(n1/2)=1.
       eev1(n1) =1.
       eev2(0)  =1.
       eev2(n2/2)=1.
       eev2(n2) =1.
       eev3(0)  =1.
       eev3(n3/2)=1.
       eev3(n3) =1.

       return 
       end
