       subroutine energy(vx,vy,vz,n1,n2,n3)

*********************************************
c compute energy cinetique in physical space
*********************************************

       implicitnone

       integer n1,n2,n3

       real vx(0:n1+1,0:n2,0:n3)
       real vy(0:n1+1,0:n2,0:n3)
       real vz(0:n1+1,0:n2,0:n3)

       real nrj
       integer ix,iy,iz

       nrj = 0.
       
        do iz=0,n3-1
        do iy=0,n2-1
        do ix=0,n1-1
           nrj=nrj+vx(ix,iy,iz)*vx(ix,iy,iz)+
     &             vy(ix,iy,iz)*vy(ix,iy,iz)+
     &             vz(ix,iy,iz)*vz(ix,iy,iz)      
       enddo
       enddo
       enddo

       print *,'nrj=',0.5*nrj/float(n1*n2*n3)

       return
       end
