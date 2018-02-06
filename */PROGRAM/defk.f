

       subroutine defk(n1,n2,n3,nh,kx,ky,kz,k2x,k2y,k2z,cz)

       implicitnone

       integer n1,n2,n3,nh

       real   kx(0:nh),ky(0:nh),kz(0:nh)
       real   k2x(0:nh),k2y(0:nh),k2z(0:nh),cz(0:nh)
 
       integer i
       real    lgx,lgy,lgz

       lgx = 1.
       lgy = 1.
       lgz = 1.

       do i=0,n1/2-1
         kx(i)=float(i)/lgx
       enddo
       do i=0,n2/2-1
         ky(i)=float(i)/lgy
       enddo
       do i=0,n3/2-1
         kz(i)=float(i)/lgz
       enddo

       do  i=n1/2,n1
         kx(i)=float(i-n1)/lgx
       enddo
       do  i=n2/2,n2
         ky(i)=float(i-n2)/lgy
       enddo
       do  i=n3/2,n3
         kz(i)=float(i-n3)/lgz
       enddo

       do i=0,n1
         k2x(i)=kx(i)*kx(i)
       enddo

       do i=0,n2
         k2y(i)=ky(i)*ky(i)
       enddo

       do i=0,n3
         k2z(i)=kz(i)*kz(i)
       enddo


       do  i=1,nh
         cz(i)=1.
       enddo

       cz(0)=0.5

       return
       end
