

       subroutine vecpro(n1,n2,n3,ax,ay,az,bx,by,bz)

*********************************************
c produit vectoriel dans l'espace physique.
*********************************************

       implicitnone

       integer n1,n2,n3

       real ax(0:n1+1,0:n2,0:n3)
       real ay(0:n1+1,0:n2,0:n3)
       real az(0:n1+1,0:n2,0:n3)

       real bx(0:n1+1,0:n2,0:n3)
       real by(0:n1+1,0:n2,0:n3)
       real bz(0:n1+1,0:n2,0:n3)

       integer ix,iy,iz
       real    alp,bet,gam

!$OMP  PARALLEL PRIVATE(alp,bet,gam)
!$OMP DO SCHEDULE(RUNTIME)

        do iz=0,n3-1
        do iy=0,n2-1
        do ix=0,n1-1
         alp=ay(ix,iy,iz)*bz(ix,iy,iz)-az(ix,iy,iz)*by(ix,iy,iz)
         bet=az(ix,iy,iz)*bx(ix,iy,iz)-ax(ix,iy,iz)*bz(ix,iy,iz)
         gam=ax(ix,iy,iz)*by(ix,iy,iz)-ay(ix,iy,iz)*bx(ix,iy,iz)
         bx(ix,iy,iz)=alp
         by(ix,iy,iz)=bet
         bz(ix,iy,iz)=gam 
       enddo
       enddo
       enddo


!$OMP END DO
!$OMP END PARALLEL 

       return
       end
