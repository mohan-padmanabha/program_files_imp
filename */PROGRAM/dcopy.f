       subroutine dcopy(n,a,incx,b,incy)

       real    a(n),b(n)
       integer n,incx,incy

!$OMP PARALLEL DO SCHEDULE(RUNTIME)
       do i=1,n
         b(i)=a(i)
       enddo
!$OMP END PARALLEL DO

       return
       end
