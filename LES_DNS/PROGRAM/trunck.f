

       subroutine trunck(n1,n2,n3,nh,vx,vy,vz,k2x,k2y,k2z,ktrunk)

       implicitnone

       integer n1,n2,n3,nh,ktrunk

       complex vx(0:n1/2,0:n2,0:n3)
       complex vy(0:n1/2,0:n2,0:n3)
       complex vz(0:n1/2,0:n2,0:n3)

       real   k2x(0:nh),k2y(0:nh),k2z(0:nh)

       real ktrunk2,kk
       integer ix,iy,iz

       ktrunk2 = float(ktrunk)*float(ktrunk)

 
!$OMP  PARALLEL PRIVATE(kk)
!$OMP DO SCHEDULE(RUNTIME)

       do iz=0,n3-1
       do iy=0,n2-1
       do ix=0,n1/2
         kk = k2x(ix)+k2y(iy)+k2z(iz)
         if (kk .gt. ktrunk2) then
           vx(ix,iy,iz)=cmplx(0.,0.)
           vy(ix,iy,iz)=cmplx(0.,0.)
           vz(ix,iy,iz)=cmplx(0.,0.)
         endif
       enddo
       enddo
       enddo

!$OMP END DO
!$OMP END PARALLEL 

      return
      end
