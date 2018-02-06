

      subroutine ffty(a,n1,n2,n3,n,wk,kx1,kx2,isign,ex2,ifax2)


**********************************
* fft complex-->complex along y
**********************************

      implicitnone


      integer n1,n2,n3,n
      real wk(4*(n+2)**2,0:1)
      complex ex1(0:n),ex2(0:n),ex3(0:n)
      integer ifax1(100),ifax2(100),ifax3(100)

      real    a(2,0:n1/2,0:n2,0:n3)
 
      integer nft,inc,jump,ikx,kx1,kx2,isign

      nft  = n3
      inc  = n1+2
      jump = (n1+2)*(n2+1)

!$OMP  PARALLEL PRIVATE(wk)
!$OMP DO SCHEDULE(RUNTIME)

      do ikx=kx1,kx2
        call cfftmlt(a(1,ikx,0,0),a(2,ikx,0,0),wk,ex2,ifax2,
     &              inc,jump,n2,nft,isign)
      enddo

!$OMP END DO
!$OMP END PARALLEL 




      return
      end
