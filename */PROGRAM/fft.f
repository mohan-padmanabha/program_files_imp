

      subroutine fft(a,n1,n2,n3,n,wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)


***********************************************
* direct fft 3D 
************************************************

      implicitnone


      integer n1,n2,n3,n

!------------------- temperton ----------------------------
      real    wk(4*(n+2)**2,0:1)
      complex ex1(0:n),ex2(0:n),ex3(0:n)
      integer ifax1(100),ifax2(100),ifax3(100)
!------------------- ESSL ---------------------------------
!      real   wk(62000+3*n1+(2*n3+257)*((n1/2+1)*n2+5))
!      integer naux,inc2x,inc3x,inc2y,inc3y
!      integer ex1,ex2,ex3,ifax1,ifax2,ifax3
!---------------------------------------------------------

      real a(0:n1+1,0:n2,0:n3)


      integer isign,ierr
      real    fnorm

      call zero(a,n1,n2,n3)

!------------------- temperton ----------------------------
      isign = -1
      fnorm = 1./(float(n2)*float(n3))
      call fftx(a,n1,n2,n3,n,wk,0,n2-1,  isign,ex1,ifax1)
      call ffty(a,n1,n2,n3,n,wk,0,n1/2-1,isign,ex2,ifax2)
      call fftz(a,n1,n2,n3,n,wk,0,n1/2-1,isign,ex3,ifax3)
      a(:,:,:) = fnorm*a(:,:,:)
!------------------- ESSL ---------------------------------
!      isign = 1
!      fnorm = 1./(float(n1)*float(n2)*float(n3))
!      naux  = 62000+3*n1+(2*n3+257)*((n1/2+1)*n2+5)
!      inc2x = n1+2
!      inc3x = inc2x*(n2+1)
!      inc2y = inc2x/2
!      inc3y = inc3x/2
!      call drcft3(a,inc2x,inc3x,a,inc2y,inc3y,
!     &            n1,n2,n3,isign,fnorm,wk,naux)
!---------------------------------------------------------

      call zerkam(a,n1,n2,n3)

      return
      end
