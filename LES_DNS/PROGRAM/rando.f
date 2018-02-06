

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine randini
c
c  randini is used to initialize the random number generators.  iseed is a
c  positive integer less than 1.0e9.  The basic random number generator is
c  taken from Press et al., Numerical Recipes, p. 199 and is based on
c  Knuth's suggestion for a portable random number generator.  mseed is
c  any large number less than m=1.0e9.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        parameter (m=1000000000,mseed=161803398,rm=1.0/m)

c       The dimension 55 is special and should not be modified; see Knuth.

        integer ma(55)
        real randar(500),xnorm
        common /randnos/ randar,irptr,nrand,ma,inext,inextp
        common /randnor/ xnorm,inorm
        save /randnos/,/randnor/
c
        nrand=157
C!!
C!        write(*,*) 'Enter random number seed < 1.0e9'
C!        read(*,*) iseed
        iseed=217
C!!
        iseed=mod(iseed,m)

c       Initialize ma(55).

        mj=mseed-iseed
        mj=mod(mj,m)
        if (mj.lt.0) mj=mj+m
        ma(55)=mj
        mk=1

c       Now initialize the rest of the table, in a slightly random order,
c       with numbers that are not especially random.

        do 10 i=1,54
           ii=mod(21*i,55)
           ma(ii)=mk
           mk=mj-mk
           if (mk.lt.0) mk=mk+m
           mj=ma(ii)
 10     continue

c       Randomize them by "warming up the generator."

        do 30 k=1,4
           do 20 i=1,55
              ma(i)=ma(i)-ma(1+mod(i+30,55))
              if (ma(i).lt.0) ma(i)=ma(i)+m
 20        continue
 30     continue
        inext=0
        inextp=31

c       Exercise generator before storing in shuffling table.

        do 40 i=1,nrand
           inext=inext+1
           if (inext.eq.56) inext=1
           inextp=inextp+1
           if (inextp.eq.56) inextp=1
           mj=ma(inext)-ma(inextp)
           if (mj.lt.0) mj=mj+m
           ma(inext)=mj
 40     continue

c       Now fill shuffling table.

        do 50 i=1,nrand
           inext=inext+1
           if (inext.eq.56) inext=1
           inextp=inextp+1
           if (inextp.eq.56) inextp=1
           mj=ma(inext)-ma(inextp)
           if (mj.lt.0) mj=mj+m
           ma(inext)=mj
           randar(i)=mj*rm
 50     continue
        inext=inext+1
        if (inext.eq.56) inext=1
        inextp=inextp+1
        if (inextp.eq.56) inextp=1
        mj=ma(inext)-ma(inextp)
        if (mj.lt.0) mj=mj+m
        ma(inext)=mj
        irptr=int(mj*rm*nrand)+1
        xnorm=0.0
        inorm=0
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine randa(x)
c
c  randa generates uniform random numbers in the interval [0,1).
c  It must be initialized with a call to randini.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        parameter (m=1000000000,mseed=161803398,rm=1.0/m)

c       The dimension 55 is special and should not be modified; see Knuth.

        integer ma(55)
        real randar(500),x
        common /randnos/ randar,irptr,nrand,ma,inext,inextp
        save /randnos/
c
c       Extract random number from shuffling table.

        x=randar(irptr)
        irptr=int(nrand*x)+1

c       Generate a new random number.

        inext=inext+1
        if (inext.eq.56) inext=1
        inextp=inextp+1
        if (inextp.eq.56) inextp=1
        mj=ma(inext)-ma(inextp)
        if (mj.lt.0) mj=mj+m
        ma(inext)=mj
        randar(irptr)=mj*rm
        return
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine randg(x)
c
c  randg generates standard normal deviates.  It uses randa, which must be
c  initialized by a call to randini.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        parameter(twopi=6.28318531)
        real x,xr,xnorm
        common /randnor/ xnorm,inorm
        save /randnor/
c
        if (inorm.eq.1) then
           x=xnorm
           inorm=0
        else
 10        call randa(x)
c          -------------
           if (x.le.0.0) go to 10
           xr=sqrt(-2.0*log(x))
           call randa(x)
c          -------------
           xnorm=xr*sin(twopi*x)
           x=xr*cos(twopi*x)
           inorm=1
        end if
        return
        end
