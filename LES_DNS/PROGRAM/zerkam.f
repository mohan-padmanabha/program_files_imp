
      subroutine zerkam(z,l1,l2,l3)

      implicitnone

      integer l1,l2,l3
      complex z(0:l1/2,0:l2,0:l3)
      integer i,j,k

C$OMP PARALLEL DO SCHEDULE(RUNTIME)
      do j = 0,l2
      do i = 0,l1/2 
        z(i,j,l3) =cmplx(0.,0.)
        z(i,j,l3/2)=cmplx(0.,0.)
      enddo
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO SCHEDULE(RUNTIME)
      do k=0,l3
      do i=0,l1/2 
        z(i,l2, k)=cmplx(0.,0.)
        z(i,l2/2,k)=cmplx(0.,0.)
      enddo
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO SCHEDULE(RUNTIME)
      do k=0,l3
      do j=0,l2
        z(l1/2,j,k)=cmplx(0.,0.)
      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
