
      subroutine fftz(a,n1,n2,n3,n,wk,kx1,kx2,isign,ex3,ifax3)


*********************************
* fft complex --> complex along z
*********************************

      implicitnone


      integer n1,n2,n3,n
      real wk(4*(n+2)**2,0:1)
      complex ex1(0:n),ex2(0:n),ex3(0:n)
      integer ifax1(100),ifax2(100),ifax3(100)

      real    a(2,0:n1/2,0:n2,0:n3)

      integer nft,inc,jump,ikx,kx1,kx2,isign

      nft  = n2
      inc  = (n1+2)*(n2+1)
      jump = n1+2

!$OMP  PARALLEL PRIVATE(wk)
!$OMP DO SCHEDULE(RUNTIME)

      do ikx=kx1,kx2

        call cfftmlt(a(1,ikx,0,0),a(2,ikx,0,0),wk,ex3,ifax3,
     &              inc,jump,n3,nft,isign)

      enddo

!$OMP END DO
!$OMP END PARALLEL 

      return
      end
