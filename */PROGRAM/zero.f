
      subroutine zero(a,l1,l2,l3)

      implicitnone

      integer l1,l2,l3
      real a(0:l1+1,0:l2,0:l3)
      integer i,j,k

C$OMP PARALLEL DO SCHEDULE(RUNTIME)
      do j=0,l2
      do k=0,l3
        a(l1,j,k)=0.
        a(l1+1,j,k)=0.
      enddo
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO SCHEDULE(RUNTIME)
      do i=0,l1 +1
      do k=0,l3
        a(i,l2,k)=0.
      enddo
      enddo
C$OMP END PARALLEL DO

C$OMP PARALLEL DO SCHEDULE(RUNTIME)
      do i=0,l1 +1
      do j=0,l2
        a(i,j,l3)=0.
      enddo
      enddo
C$OMP END PARALLEL DO

      return
      end
