
      subroutine initfft(l1,l2,l3,lm,
     &     wk,ex1,ex2,ex3,ifax1,ifax2,ifax3)



      integer l1,l2,l3,lm
      real wk(4*(lm+2)**2,0:1)
      complex ex1(0:lm),ex2(0:lm),ex3(0:lm)
      integer ifax1(100),ifax2(100),ifax3(100)

      call fftfax(l1,ifax1,ex1)
      call cftfax(l2,ifax2,ex2)
      call cftfax(l3,ifax3,ex3)

      return
      end
