

      subroutine fftx(a,n1,n2,n3,n,wk,iy1,iy2,isign,ex1,ifax1)


**********************************
* fft real --> complex along x
**********************************


      implicitnone


      integer n1,n2,n3,n

      real wk(4*(n+2)**2,0:1)
      complex ex1(0:n),ex2(0:n),ex3(0:n)
      integer ifax1(100),ifax2(100),ifax3(100)

      real    a(0:n1+1,0:n2,0:n3)

      integer nft,inc,jump,iy,iy1,iy2,isign

      nft  = n3
      inc  = 1
      jump = (n1+2)*(n2+1)

!==================================================================
! a : tableau a transformer
! wk: tableau de travail
! ex1: vecteur initialisée
! ifax1 :  vecteur initialisée
! inc : distance entre deux element à transformer
! jump : distance entre les elements des vecteurs à tranformer
! l1 : taille de la fft
! nft : nombre de fft simultanée
! isign : signe de la FFT
!==================================================================


!$OMP  PARALLEL PRIVATE(wk)
!$OMP DO SCHEDULE(RUNTIME)

      do  iy=iy1,iy2
       call rfftmlt(a(0,iy,0),wk,ex1,ifax1,inc,jump,n1,nft,isign)
      enddo

!$OMP END DO
!$OMP END PARALLEL 


      return
      end
