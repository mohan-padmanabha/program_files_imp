       subroutine initial(a,n,val)

       real    a(n)
       integer n

!$OMP PARALLEL DO SCHEDULE(RUNTIME)
 

       do i=1,n
         a(i)=val
       enddo

!$OMP END PARALLEL DO

       return
       end
