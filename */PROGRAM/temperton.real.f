



















 




















 
 










































































































































































































 













































































































































 






































































************************************************************************
      subroutine cfftmlt(a,b,w,ex,ifax,inc,jump,n,nft,isign)
      real a(*),b(*),w(*),ex(*)
      integer ifax(*),inc,jump,n,nft,isign
      call gpfa(a,b,ex,inc,jump,n,ifax(2),ifax(3),ifax(5),nft,isign)
      return
      end
************************************************************************
      subroutine cftfax(n,ifax,ex)
      real ex(*)
      integer n,ifax(*)
      call primfac(n,ifax)
      call trig37(ex,n,ifax(2),ifax(3),ifax(5))
      return
      end
************************************************************************
      subroutine fax(ifax,n,mode)
c                                                                      *
c c06-summatiom
c c06-summation of series                                      b6.1/3  *
c                                                                      *
c                                                              fftrig  *
c                                                              fax     *
c                                                                      *
c                                                                      *
c sup
c subprogram                    fftrig                                 *
c                               fax                                    *
c                                                                      *
c purpose          setup routines for fft packages                     *
c                                                                      *
c                                                                      *
c version          cyber                         cray-1                *
c                                                                      *
c                  jan 1979 original             jan 1979 original     *
c                                                                      *
c usage                                                                *
c                  call fftrig(trigs,n,3)                              *
c                  call fax   (ifax ,n,3)                              *
c                                                                      *
c arguments        1.dimension                                         *
c                       trigs(dimension 3*n/2 - add 1 if n/2 is odd)   *
c                       ifax(10)                                       *
c                                                                      *
c                  2.input                                             *
c                      n - the lenght of the transforms to be performed*
c                          n must be even.                             *
c                          the number of words of ifax used increases  *
c                          logarithmically with n.                     *
c                          ifax(10) suffices for practical purposes.   *
c                          (transforms of lenght at least 10000)       *
c                                                                      *
c                  3.output                                            *
c                      trigs - fftrig returns an array of trigonometric*
c                              function values subsequently used by    *
c                              fft routines.                           *
c                      ifax  - fax factorizes n/2 into a product of    *
c                              4's and 2's and higher prime numbers.   *
c                              ifax(1) contains the number of factors. *
c                              and the factors themselves are stored   *
c                              in ascending order in ifax(2),ifax(3).. *
c                              if fax is called with n odd ,ifax(1)    *
c                              is set to -99(error condition) and no   *
c                              factorization is done.                  *
c                                                                      *
c write up         none                                                *
c                                                                      *
c entry points           fftrig,  fax                                  *
c                                                                      *
c common blocks    none                                                *
c i/o              none                                                *
c precision        single                                              *
c other routines   none                                                *
c       required                                                       *
c 7/80                     fftrig-1                                    *
c                                                                      *
c***********************************************************************
c                                                                      *
c co6-summation of series                                       b6.1/3 *
c                                                                      *
c                                                              fftrig  *
c                                                              fax     *
c                                                                      *
c acsses (object)  cyber:                                              *
c                           attach,eclib.                              *
c                           ldset(lib=eclib)                           *
c                  cray 1:                                             *
c                           ldr(lib=eclib...)                          *
c                                                                      *
c access (source)           attach,oldpl,eclibpl                       *
c                                                                      *
c                  cyber :         %define cyber                       *
c                  cray:           %define cray                        *
c                                  %c   fftrig,   fax                  *
c                                                                      *
c language         fortran                                             *
c                                                                      *
c specialist       clive temperton                                     *
c                                                                      *
c history          written by c.temperton      jan     1979            *
c                                                                      *
c algorithm                                                            *
c references                                                           *
c                                                                      *
c object size               fftrig  fax  (octal words)                 *
c                  cyber:     145   127                                *
c                  cray :     221   157                                *
c                                                                      *
c                                                                      *
c accuracy                                                             *
c                                                                      *
c timing                                                               *
c                                                                      *
c portability      standard fortran                                    *
c                                                                      *
c system routines  none                                                *
c        required                                                      *
c
c 7/80                      fftrig-2                                   *
c
c***********************************************************************
c     end
      dimension ifax(10)
      nn=n
      if (iabs(mode).eq.1) go to 10
      if (iabs(mode).eq.8) go to 10
      nn=n/2
      if ((nn+nn).eq.n)  go to 10
      ifax(1)=-99
      return
   10 k=1
c     test for factors of 4
   20 if (mod(nn,4).ne.0) go to 30
      k=k+1
      ifax(k)=4
      nn=nn/4
      if (nn.eq.1) go to 80
      go to 20
c     test for extra factor of 2
   30 if (mod(nn,2).ne.0) go to 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      if (nn.eq.1) go to 80
c     test for factors of 3
   40 if (mod(nn,3).ne.0) go to 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      if (nn.eq.1) go to 80
      go to 40
 
c     now find remaining factors
   50 l=5
      inc=2
c     inc alternatively takes on values 2 and 4
   60 if (mod(nn,l).ne.0) go to 70
      k=k+1
      ifax(k)=l
      nn=nn/l
      if (nn.eq.1) go to 80
      go to 60
   70 l=l+inc
      inc=6-inc
      go to 60
   80 ifax(1)=k-1
c     ifax(1) contains number of factors
      nfax=ifax(1)
c     sort factors into ascending order
      if (nfax.eq.1) go to 110
      do 100 ii=2,nfax
      istop=nfax+2-ii
      do 90 i=2,istop
      if (ifax(i+1).ge.ifax(i)) go to 90
      item=ifax(i)
      ifax(i)=ifax(i+1)
      ifax(i+1)=item
   90 continue
  100 continue
  110 continue
      return
      end
************************************************************************
      subroutine fft99(a,work,trigs,ifax,inc,jump,n,lot,isign)
*                                                                      *
* c06-summation of series                                       b6.1/3 *
*                                                                      *
*                                                               fft99  *
*                                                               fft991 *
*                                                                      *
*                                                                      *
* subprogram                     fft99                                 *
*                                fft991                                *
*                                                                      *
* purpose          perform multiple fast fourier transforms            *
*                                                                      *
*                                                                      *
* version          cyber                         cray-1                *
*                                                                      *
*                  jan 1979 original             jan 1979 original     *
*                                                                      *
* usage                                                                *
*                  call fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)   *
*                  call fft991(a,work,trigs,ifax,inc,jump,n,m,isign)   *
*                                                                      *
* arguments        1.dimension                                         *
*                       a(idim),work((n+1)*m),trigs(3*n/2),ifax(10)    *
*                       work is a work array                           *
*                                                                      *
*                  2.input                                             *
*                      a - an array containing the input data or       *
*                          coefficient vectors.                        *
*                         this array is overwritten by the results.    *
*                      trigs and ifax - arrays set up by fftrig and fax*
*                                     - see writeup of fftrig and fax  *
*                      inc - the word increment between successive     *
*                           elements of each data or coefficient vector*
*                           e.g. inc=1 for consecutively stored data.  *
*                      jump - the word increment between the first     *
*                            elements of successive data or coefficient*
*                            vectors.                                  *
*                      n - the length of each transform. (see note x)  *
*                      m - the number of transforms to be done         *
*                          simultaneously.                             *
*                      isign - +1 for a transform from fourier         *
*                              coefficients to data values.            *
*                              -1 for a transform from data values     *
*                              to fourier coefficients.                *
*                                                                      *
*                  3.output                                            *
*                      a - contains either the coefficients or the     *
*                          data values,depending on isign.             *
*                          in each case n independent quantities       *
*                          occupy n+2 words.   the coefficients are    *
*                          stored as successive pairs of real and      *
*                          imaginary parts -                           *
*                          a(k),b(k) , k=0,1,...n/2                    *
*                          b(0) and b(n/2) are stored although they    *
*                          must be 0.                                  *
*                      for fft99 the data is stored with explicit      *
*                          periodicity -                               *
*                          x(n-1),x(0),x(1),....x(n-1),x(0)            *
*                      for fft991 the data appears as -                *
*                          x(0),x(1),x(2),......x(n-1),0,0             *
*                                                                      *
* notes            1. on cray-1, arrange data so that jump is not a    *
*                     multiple of 8 (to avoid memory bank conflicts)   *
*                                                                      *
* write up         computer bulletin b6.6/1                            *
*                                                                      *
* entry points        fft99,fft991                                     *
*                                                                      *
* common blocks    none                                                *
*                                                                      *
* i/o              none                                                *
*                                                                      *
* precision        single                                              *
*                                                                      *
* other routines   fft99a,fft99b,vpassm          (cy)                  *
*       required   cal99,cpass                   (cr)                  *
*                                                                      *
*                                                                      *
* 7/80                      fft99-1                                    *
*                                                                      *
************************************************************************
*                                                                      *
* c06-summation of series                                       b6.1/3 *
*                                                                      *
*                                                               fft99  *
*                                                               fft991 *
*                                                                      *
* access (object)  cyber:                                              *
*                           attach,eclib.                              *
*                           ldset(lib=eclib)                           *
*                  cray 1:                                             *
*                           ldr(lib=eclib...)                          *
*                                                                      *
* access (source)           attach,oldpl,eclibpl                       *
*                                                                      *
*                  cyber :         %define cyber                       *
*                  cray:           %define cray                        *
*                                  %c    fft99,fft991                  *
*                                                                      *
* language         fortran                                             *
*                  but cray implementation of pass is in cal           *
*                                                                      *
* specialist       clive temperton                                     *
*                                                                      *
* history          written by c.temperton      jan     1979            *
*                                                                      *
* algorithm        the algorithm is the self-sorting (temperton)       *
*                  version of the fast fourier transform               *
*                                                                      *
* references       ecmwf technical report no.3                         *
*                  ecmwf internal report no.21 -   c.temperton         *
*                                                                      *
* object size               fft991  fft99  (octal words)               *
*                  cyber:    2665    2676                              *
*                  cray :    1250    1260                              *
*                                                                      *
*                                                                      *
* accuracy                                                             *
*                                                                      *
* timing           vectorization is on vectors of length m.      (cr)  *
*                  hence timing is strongly dependent on m.            *
*                  time per transform on cray-1 (microseconds)         *
*                  n    m=4    m=16    m=64                            *
*                 64     46      17      10                            *
*                128     81      33      21                            *
*                180    150      58      37                            *
*                192    149      58      36                            *
*                240    192      76      49                            *
*                256    191      76      49                            *
*                288    219      89      58                            *
*                300    253     102      68                            *
*                320    248     101      66                            *
*                360    286     118      79                            *
*               1024    898     359     238                            *
*                                                                      *
* portability      standard fortran                                    *
*                  standard cal  (cr)                                  *
*                                                                      *
* system routines  none                                                *
*        required                                                      *
*                                                                      *
* 7/80                      fft99-1                                    *
*                                                                      *
************************************************************************
c
c     routine 'fft99' - multiple fast real periodic transform
c     corresponding to old scalar routine fft9
c     procedure used to convert to half-length complex transform
c     is given by cooley, lewis ' welch (j. sound vib., vol. 12
c     (1970), 315-337)
c
c     a is the array containing input ' output data
c     work is an area of size (n+1)*lot
c     trigs is a previously prepared list of trig function values
c     ifax is a previously prepared list of factors of n/2
c     inc is the increment within each data 'vector'
c         (e.g. inc=1 for consecutively stored data)
c     jump is the increment between the start of each data vector
c     n is the length of the data vectors
c     lot is the number of data vectors
c     isign = +1 for transform from spectral to gridpoint
c           = -1 for transform from gridpoint to spectral
c
c     ordering of coefficients:
c         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
c         where b(0)=b(n/2)=0; (n+2) locations required
c
c     ordering of data:
c         x(n-1),x(0),x(1),x(2),...,x(n),x(0)
c         i.e. explicit cyclic continuity; (n+2) locations required
c
c     vectorization is achieved on cray by doing the transforms in
c     parallel
c
c     *** n.b. n is assumed to be an even number
c
c     definition of transforms:
c     -------------------------
c
c     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
c         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
c
c     isign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
c               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
c
c
c     line following next is not routine header(only comment)
c
c     routine fft99(a,work,trigs,ifax,inc,jump,n,lot,isign)
c
c     end
      dimension a(n),work(n),trigs(n),ifax(1)
c
      nfax=ifax(1)
      if(nfax.le.0) go to 99
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 30
c
c     if necessary, transfer data to work area
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=inc+1
      jbase=1
      do 20 l=1,lot
      i=ibase
      j=jbase
cdir$ ivdep
      do 10 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   10 continue
      ibase=ibase+jump
      jbase=jbase+nx
   20 continue
c
      igo=60
      go to 40
c
c     preprocessing (isign=+1)
c     ------------------------
c
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
c
c     complex transform
c     -----------------
c
   40 continue
      ia=inc+1
      la=1
      do 80 k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
      call vpassm(a(ia),a(ia+inc),work(1),work(2),trigs,
     *   ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm(work(1),work(2),a(ia),a(ia+inc),trigs,
     *    2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
c
      if (isign.eq.-1) go to 130
c
c     if necessary, transfer data from work area
      if (mod(nfax,2).eq.1) go to 110
      ibase=1
      jbase=ia
      do 100 l=1,lot
      i=ibase
      j=jbase
cdir$ ivdep
      do 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 continue
      ibase=ibase+nx
      jbase=jbase+jump
  100 continue
c
c     fill in cyclic boundary points
  110 continue
      ia=1
      ib=n*inc+1
cdir$ ivdep
      do 120 l=1,lot
      a(ia)=a(ib)
      a(ib+inc)=a(ia+inc)
      ia=ia+jump
      ib=ib+jump
  120 continue
      go to 140
c
c     postprocessing (isign=-1):
c     --------------------------
c
  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)
c
  140 continue
      return
c  ** error exit   ifax(1) le 0 **
99    print *,' fft99 called but factors not supplied '
      call abort
      end
************************************************************************
      subroutine rfftmlt(a,work,trigs,ifax,inc,jump,n,lot,isign)
c
c     routine 'fft991' - multiple real/half-complex periodic
c     fast fourier transform
c
*                                                                      *
* c06-summation of series                                       b6.1/3 *
*                                                                      *
*                                                               fft99  *
*                                                               fft991 *
*                                                                      *
*                                                                      *
* subprogram                     fft99                                 *
*                                fft991                                *
*                                                                      *
* purpose          perform multiple fast fourier transforms            *
*                                                                      *
*                                                                      *
* version          cyber                         cray-1                *
*                                                                      *
*                  jan 1979 original             jan 1979 original     *
*                                                                      *
* usage                                                                *
*                  call fft99 (a,work,trigs,ifax,inc,jump,n,m,isign)   *
*                  call fft991(a,work,trigs,ifax,inc,jump,n,m,isign)   *
*                                                                      *
* arguments        1.dimension                                         *
*                       a(idim),work((n+1)*m),trigs(3*n/2),ifax(10)    *
*                       work is a work array                           *
*                                                                      *
*                  2.input                                             *
*                      a - an array containing the input data or       *
*                          coefficient vectors.                        *
*                         this array is overwritten by the results.    *
*                      trigs and ifax - arrays set up by fftrig and fax*
*                                     - see writeup of fftrig and fax  *
*                      inc - the word increment between successive     *
*                           elements of each data or coefficient vector*
*                           e.g. inc=1 for consecutively stored data.  *
*                      jump - the word increment between the first     *
*                            elements of successive data or coefficient*
*                            vectors.                                  *
*                      n - the length of each transform. (see note x)  *
*                      m - the number of transforms to be done         *
*                          simultaneously.                             *
*                      isign - +1 for a transform from fourier         *
*                              coefficients to data values.            *
*                              -1 for a transform from data values     *
*                              to fourier coefficients.                *
*                                                                      *
*                  3.output                                            *
*                      a - contains either the coefficients or the     *
*                          data values,depending on isign.             *
*                          in each case n independent quantities       *
*                          occupy n+2 words.   the coefficients are    *
*                          stored as successive pairs of real and      *
*                          imaginary parts -                           *
*                          a(k),b(k) , k=0,1,...n/2                    *
*                          b(0) and b(n/2) are stored although they    *
*                          must be 0.                                  *
*                      for fft99 the data is stored with explicit      *
*                          periodicity -                               *
*                          x(n-1),x(0),x(1),....x(n-1),x(0)            *
*                      for fft991 the data appears as -                *
*                          x(0),x(1),x(2),......x(n-1),0,0             *
*                                                                      *
* notes            1. on cray-1, arrange data so that jump is not a    *
*                     multiple of 8 (to avoid memory bank conflicts)   *
*                                                                      *
* write up         computer bulletin b6.6/1                            *
*                                                                      *
* entry points        fft99,fft991                                     *
*                                                                      *
* common blocks    none                                                *
*                                                                      *
* i/o              none                                                *
*                                                                      *
* precision        single                                              *
*                                                                      *
* other routines   fft99a,fft99b,vpassm          (cy)                  *
*       required   cal99,cpass                   (cr)                  *
*                                                                      *
*                                                                      *
* 7/80                      fft99-1                                    *
*                                                                      *
************************************************************************
*                                                                      *
* c06-summation of series                                       b6.1/3 *
*                                                                      *
*                                                               fft99  *
*                                                               fft991 *
*                                                                      *
* access (object)  cyber:                                              *
*                           attach,eclib.                              *
*                           ldset(lib=eclib)                           *
*                  cray 1:                                             *
*                           ldr(lib=eclib...)                          *
*                                                                      *
* access (source)           attach,oldpl,eclibpl                       *
*                                                                      *
*                  cyber :         %define cyber                       *
*                  cray:           %define cray                        *
*                                  %c    fft99,fft991                  *
*                                                                      *
* language         fortran                                             *
*                  but cray implementation of pass is in cal           *
*                                                                      *
* specialist       clive temperton                                     *
*                                                                      *
* history          written by c.temperton      jan     1979            *
*                                                                      *
* algorithm        the algorithm is the self-sorting (temperton)       *
*                  version of the fast fourier transform               *
*                                                                      *
* references       ecmwf technical report no.3                         *
*                  ecmwf internal report no.21 -   c.temperton         *
*                                                                      *
* object size               fft991  fft99  (octal words)               *
*                  cyber:    2665    2676                              *
*                  cray :    1250    1260                              *
*                                                                      *
*                                                                      *
* accuracy                                                             *
*                                                                      *
* timing           vectorization is on vectors of length m.      (cr)  *
*                  hence timing is strongly dependent on m.            *
*                  time per transform on cray-1 (microseconds)         *
*                  n    m=4    m=16    m=64                            *
*                 64     46      17      10                            *
*                128     81      33      21                            *
*                180    150      58      37                            *
*                192    149      58      36                            *
*                240    192      76      49                            *
*                256    191      76      49                            *
*                288    219      89      58                            *
*                300    253     102      68                            *
*                320    248     101      66                            *
*                360    286     118      79                            *
*               1024    898     359     238                            *
*                                                                      *
* portability      standard fortran                                    *
*                  standard cal  (cr)                                  *
*                                                                      *
* system routines  none                                                *
*        required                                                      *
*                                                                      *
* 7/80                      fft99-1                                    *
*                                                                      *
************************************************************************
c
c     same as fft99 except that ordering of data corresponds to
c     that in mrfft2
c
c     procedure used to convert to half-length complex transform
c     is given by cooley, lewis ' welch (j. sound vib., vol. 12
c     (1970), 315-337)
c
c     a is the array containing input ' output data
c     work is an area of size (n+1)*lot
c     trigs is a previously prepared list of trig function values
c     ifax is a previously prepared list of factors of n/2
c     inc is the increment within each data 'vector'
c         (e.g. inc=1 for consecutively stored data)
c     jump is the increment between the start of each data vector
c     n is the length of the data vectors
c     lot is the number of data vectors
c     isign = +1 for transform from spectral to gridpoint
c           = -1 for transform from gridpoint to spectral
c
c     ordering of coefficients:
c         a(0),b(0),a(1),b(1),a(2),b(2),...,a(n/2),b(n/2)
c         where b(0)=b(n/2)=0; (n+2) locations required
c
c     ordering of data:
c         x(0),x(1),x(2),...,x(n-1)
c
c     vectorization is achieved on cray by doing the transforms in
c     parallel
c
c     *** n.b. n is assumed to be an even number
c
c     definition of transforms:
c     -------------------------
c
c     isign=+1: x(j)=sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/n))
c         where c(k)=a(k)+i*b(k) and c(n-k)=a(k)-i*b(k)
c
c     isign=-1: a(k)=(1/n)*sum(j=0,...,n-1)(x(j)*cos(2*j*k*pi/n))
c               b(k)=-(1/n)*sum(j=0,...,n-1)(x(j)*sin(2*j*k*pi/n))
c
c     routine fft991(a,work,trigs,ifax,inc,jump,n,lot,isign)
c     end
      dimension a(n),work(n),trigs(n),ifax(1)
c
      nfax=ifax(1)
      if(nfax.le.0) go to 99
      nx=n+1
      nh=n/2
      ink=inc+inc
      if (isign.eq.+1) go to 30
c
c     if necessary, transfer data to work area
      igo=50
      if (mod(nfax,2).eq.1) goto 40
      ibase=1
      jbase=1
      do 20 l=1,lot
      i=ibase
      j=jbase
cdir$ ivdep
      do 10 m=1,n
      work(j)=a(i)
      i=i+inc
      j=j+1
   10 continue
      ibase=ibase+jump
      jbase=jbase+nx
   20 continue
c
      igo=60
      go to 40
c
c     preprocessing (isign=+1)
c     ------------------------
c
   30 continue
      call fft99a(a,work,trigs,inc,jump,n,lot)
      igo=60
c
c     complex transform
c     -----------------
c
   40 continue
      ia=1
      la=1
      do 80 k=1,nfax
      if (igo.eq.60) go to 60
   50 continue
      call vpassm(a(ia),a(ia+inc),work(1),work(2),trigs,
     *   ink,2,jump,nx,lot,nh,ifax(k+1),la)
      igo=60
      go to 70
   60 continue
      call vpassm(work(1),work(2),a(ia),a(ia+inc),trigs,
     *    2,ink,nx,jump,lot,nh,ifax(k+1),la)
      igo=50
   70 continue
      la=la*ifax(k+1)
   80 continue
c
      if (isign.eq.-1) go to 130
c
c     if necessary, transfer data from work area
      if (mod(nfax,2).eq.1) go to 110
      ibase=1
      jbase=1
      do 100 l=1,lot
      i=ibase
      j=jbase
cdir$ ivdep
      do 90 m=1,n
      a(j)=work(i)
      i=i+1
      j=j+inc
   90 continue
      ibase=ibase+nx
      jbase=jbase+jump
  100 continue
c
c     fill in zeros at end
  110 continue
      ib=n*inc+1
cdir$ ivdep
      do 120 l=1,lot
      a(ib)=0.0
      a(ib+inc)=0.0
      ib=ib+jump
  120 continue
      go to 140
c
c     postprocessing (isign=-1):
c     --------------------------
c
  130 continue
      call fft99b(work,a,trigs,inc,jump,n,lot)
c
  140 continue
      return
c   **  error     ifax(1) le 0  **
99    print *,' fft991 called but factors not supplied '
      call abort
      end
************************************************************************
      subroutine fft99a(a,work,trigs,inc,jump,n,lot)
c
c     routine fft99a - preprocessing step for fft99, isign=+1
c     (spectral to gridpoint transform)
c
c     line following next is not routine header(only comment)
c
c     routine fft99a(a,work,trigs,inc,jump,n,lot)
c     end
      dimension a(n),work(n),trigs(n)
      nh=n/2
      nx=n+1
      ink=inc+inc
c
c     a(0) ' a(n/2)
      ia=1
      ib=n*inc+1
      ja=1
      jb=2
cdir$ ivdep
      do 10 l=1,lot
      work(ja)=a(ia)+a(ib)
      work(jb)=a(ia)-a(ib)
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   10 continue
c
c     remaining wavenumbers
      iabase=2*inc+1
      ibbase=(n-2)*inc+1
      jabase=3
      jbbase=n-1
c
      do 30 k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
cdir$ ivdep
      do 20 l=1,lot
      work(ja)=(a(ia)+a(ib))-
     *    (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(jb)=(a(ia)+a(ib))+
     *    (s*(a(ia)-a(ib))+c*(a(ia+inc)+a(ib+inc)))
      work(ja+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))+
     *    (a(ia+inc)-a(ib+inc))
      work(jb+1)=(c*(a(ia)-a(ib))-s*(a(ia+inc)+a(ib+inc)))-
     *    (a(ia+inc)-a(ib+inc))
      ia=ia+jump
      ib=ib+jump
      ja=ja+nx
      jb=jb+nx
   20 continue
      iabase=iabase+ink
      ibbase=ibbase-ink
      jabase=jabase+2
      jbbase=jbbase-2
   30 continue
c
      if (iabase.ne.ibbase) go to 50
c     wavenumber n/4 (if it exists)
      ia=iabase
      ja=jabase
cdir$ ivdep
      do 40 l=1,lot
      work(ja)=2.0*a(ia)
      work(ja+1)=-2.0*a(ia+inc)
      ia=ia+jump
      ja=ja+nx
   40 continue
c
   50 continue
      return
      end
***********************************************************************
      subroutine fft99b(work,a,trigs,inc,jump,n,lot)
c
c     fft99b - postprocessing step for fft99, isign=-1
c     (gridpoint to spectral transform)
c
c
c     line follwing next is not routine header(only comment)
c
c     routine fft99b(work,a,trigs,inc,jump,n,lot)
c     end
      dimension work(n),a(n),trigs(n)
c
      nh=n/2
      nx=n+1
      ink=inc+inc
c
c     a(0) ' a(n/2)
      scale=1.0/float(n)
      ia=1
      ib=2
      ja=1
      jb=n*inc+1
cdir$ ivdep
      do 10 l=1,lot
      a(ja)=scale*(work(ia)+work(ib))
      a(jb)=scale*(work(ia)-work(ib))
      a(ja+inc)=0.0
      a(jb+inc)=0.0
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   10 continue
c
c     remaining wavenumbers
      scale=0.5*scale
      iabase=3
      ibbase=n-1
      jabase=2*inc+1
      jbbase=(n-2)*inc+1
c
      do 30 k=3,nh,2
      ia=iabase
      ib=ibbase
      ja=jabase
      jb=jbbase
      c=trigs(n+k)
      s=trigs(n+k+1)
cdir$ ivdep
      do 20 l=1,lot
      a(ja)=scale*((work(ia)+work(ib))
     *   +(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(jb)=scale*((work(ia)+work(ib))
     *   -(c*(work(ia+1)+work(ib+1))+s*(work(ia)-work(ib))))
      a(ja+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1)))
     *    +(work(ib+1)-work(ia+1)))
      a(jb+inc)=scale*((c*(work(ia)-work(ib))-s*(work(ia+1)+work(ib+1)))
     *    -(work(ib+1)-work(ia+1)))
      ia=ia+nx
      ib=ib+nx
      ja=ja+jump
      jb=jb+jump
   20 continue
      iabase=iabase+2
      ibbase=ibbase-2
      jabase=jabase+ink
      jbbase=jbbase-ink
   30 continue
c
      if (iabase.ne.ibbase) go to 50
c     wavenumber n/4 (if it exists)
      ia=iabase
      ja=jabase
      scale=2.0*scale
cdir$ ivdep
      do 40 l=1,lot
      a(ja)=scale*work(ia)
      a(ja+inc)=-scale*work(ia+1)
      ia=ia+nx
      ja=ja+jump
   40 continue
c
   50 continue
      return
      end
************************************************************************
      subroutine fftfax(n,ifax,ex)
      integer ifax(*),n
      real ex(*)
      call set99(ex,ifax,n)
      return
      end
*************************************************************************
      subroutine fftrig (trigs,n,mode)
      dimension trigs(n)
c    fftrig returns an array of trigonometric function values
c           subsequently used by    f f t    routines
c    see comments in routine    f a x
c     end
      pi=2.0*asin(1.d0)
      imode=iabs(mode)
      nn=n
      if (imode.gt.1.and.imode.lt.6) nn=n/2
      del=(pi+pi)/float(nn)
      l=nn+nn
      do 10 i=1,l,2
 
      angle=0.5*float(i-1)*del
      trigs(i)=cos(angle)
      trigs(i+1)=sin(angle)
   10 continue
      if (imode.eq.1) return
      if (imode.eq.8) return
      del=0.5*del
      nh=(nn+1)/2
      l=nh+nh
      la=nn+nn
      do 20 i=1,l,2
      angle=0.5*float(i-1)*del
      trigs(la+i)=cos(angle)
      trigs(la+i+1)=sin(angle)
   20 continue
      if (imode.le.3) return
      del=0.5*del
      la=la+nn
      if (mode.eq.5) go to 40
      do 30 i=2,nn
      angle=float(i-1)*del
      trigs(la+i)=2.0*sin(angle)
   30 continue
      return
   40 continue
      del=0.5*del
      do 50 i=2,n
      angle=float(i-1)*del
      trigs(la+i)=sin(angle)
   50 continue
      return
      end
************************************************************************
*        subroutine 'gpfa'
*        self-sorting in-place generalized prime factor (complex) fft
*
*        call gpfa(a,b,trigs,inc,jump,n,ip,iq,ir,lot,isign)
*
*        a is first real input/output vector
*        b is first imaginary input/output vector
*        trigs is a table of twiddle factors, precalculated
*              by calling subroutine 'trig37'
*        inc is the increment within each data vector
*        jump is the increment between data vectors
*        n is the length of the transforms:
*          -----------------------------------
*            n = (2**ip) * (3**iq) * (5**ir)
*          -----------------------------------
*        lot is the number of transforms
*        isign = +1 for forward transform
*              = -1 for inverse transform
*
*        written by clive temperton ecmwf 1990
*
*----------------------------------------------------------------------
*
*        definition of transform
*        -----------------------
*
*        x(j) = sum(k=0,...,n-1)(c(k)*exp(isign*2*i*j*k*pi/n))
*
*---------------------------------------------------------------------
*
*        for a mathematical development of the algorithm used,
*        see:
*
*        c temperton : "a generalized prime factor fft algorithm
*          for any n = (2**p)(3**q)(5**r)",
*        submitted to siam j. sci. stat. comp.
*
*----------------------------------------------------------------------
*
      subroutine gpfa(a,b,trigs,inc,jump,n,ip,iq,ir,lot,isign)
*
      dimension a(*), b(*), trigs(*)
*
      i = 1
      if (ip.gt.0) then
         call gpfa2f(a,b,trigs,inc,jump,n,ip,lot,isign)
         i = i + 2 * ( 2**ip)
      endif
      if (iq.gt.0) then
         call gpfa3f(a,b,trigs(i),inc,jump,n,iq,lot,isign)
         i = i + 2 * (3**iq)
      endif
      if (ir.gt.0) then
         call gpfa5f(a,b,trigs(i),inc,jump,n,ir,lot,isign)
      endif
*
      return
      end
*     fortran version of *gpfa2* -
*     radix-2 section of self-sorting, in-place, generalized pfa
*     central radix-2 and radix-8 passes included
*      so that transform length can be any power of 2
*
*-------------------------------------------------------------------
*
      subroutine gpfa2f(a,b,trigs,inc,jump,n,mm,lot,isign)
      dimension a(*), b(*), trigs(*)
*
      n2 = 2**mm
      inq = n/n2
      jstepx = (n2-n) * inc
      ninc = n * inc
      ink = inc * inq
*
      m2 = 0
      m8 = 0
      if (mod(mm,2).eq.0) then
         m = mm/2
      else if (mod(mm,4).eq.1) then
         m = (mm-1)/2
         m2 = 1
      else if (mod(mm,4).eq.3) then
         m = (mm-3)/2
         m8 = 1
      endif
      mh = (m+1)/2
*
      nblox = 1 + (lot-1)/64
      left = lot
      s = float(isign)
      istart = 1
*
*  loop on blocks of 64 transforms
*  -------------------------------
      do 500 nb = 1 , nblox
*
      if (left.le.64) then
         nvex = left
      else if (left.lt.128) then
         nvex = left/2
      else
         nvex = 64
      endif
      left = left - nvex
*
      la = 1
*
*  loop on type i radix-4 passes
*  -----------------------------
      mu = mod(inq,4)
      if (isign.eq.-1) mu = 4 - mu
      ss = 1.0
      if (mu.eq.3) ss = -1.0
*
      if (mh.eq.0) go to 200
*
      do 160 ipass = 1 , mh
      jstep = (n*inc) / (4*la)
      jstepl = jstep - ninc
*
*  k = 0 loop (no twiddle factors)
*  -------------------------------
      do 120 jjj = 0 , (n-1)*inc , 4*jstep
      ja = istart + jjj
*
*     "transverse" loop
*     -----------------
      do 115 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      j = 0
*
*  loop across transforms
*  ----------------------
cdir$ ivdep, shortloop
      do 110 l = 1 , nvex
      t0 = a(ja+j) + a(jc+j)
      t2 = a(ja+j) - a(jc+j)
      t1 = a(jb+j) + a(jd+j)
      t3 = ss * ( a(jb+j) - a(jd+j) )
      u0 = b(ja+j) + b(jc+j)
      u2 = b(ja+j) - b(jc+j)
      u1 = b(jb+j) + b(jd+j)
      u3 = ss * ( b(jb+j) - b(jd+j) )
      a(ja+j) = t0 + t1
      a(jc+j) = t0 - t1
      b(ja+j) = u0 + u1
      b(jc+j) = u0 - u1
      a(jb+j) = t2 - u3
      a(jd+j) = t2 + u3
      b(jb+j) = u2 + t3
      b(jd+j) = u2 - t3
      j = j + jump
  110 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  115 continue
  120 continue
*
*  finished if n2 = 4
*  ------------------
      if (n2.eq.4) go to 490
      kk = 2 * la
*
*  loop on nonzero k
*  -----------------
      do 150 k = ink , jstep-ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s*trigs(3*kk+2)
*
*  loop along transform
*  --------------------
      do 140 jjj = k , (n-1)*inc , 4*jstep
      ja = istart + jjj
*
*     "transverse" loop
*     -----------------
      do 135 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      j = 0
*
*  loop across transforms
*  ----------------------
cdir$ ivdep,shortloop
      do 130 l = 1 , nvex
      t0 = a(ja+j) + a(jc+j)
      t2 = a(ja+j) - a(jc+j)
      t1 = a(jb+j) + a(jd+j)
      t3 = ss * ( a(jb+j) - a(jd+j ) )
      u0 = b(ja+j) + b(jc+j)
      u2 = b(ja+j) - b(jc+j)
      u1 = b(jb+j) + b(jd+j)
      u3 = ss * ( b(jb+j) - b(jd+j) )
      a(ja+j) = t0 + t1
      b(ja+j) = u0 + u1
      a(jb+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jb+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jc+j) = co2*(t0-t1) - si2*(u0-u1)
      b(jc+j) = si2*(t0-t1) + co2*(u0-u1)
      a(jd+j) = co3*(t2+u3) - si3*(u2-t3)
      b(jd+j) = si3*(t2+u3) + co3*(u2-t3)
      j = j + jump
  130 continue
*-----( end of loop across transforms )
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  135 continue
  140 continue
*-----( end of loop along transforms )
      kk = kk + 2*la
  150 continue
*-----( end of loop on nonzero k )
      la = 4*la
  160 continue
*-----( end of loop on type i radix-4 passes)
*
*  central radix-2 pass
*  --------------------
  200 continue
      if (m2.eq.0) go to 300
*
      jstep = (n*inc) / (2*la)
      jstepl = jstep - ninc
*
*  k=0 loop (no twiddle factors)
*  -----------------------------
      do 220 jjj = 0 , (n-1)*inc , 2*jstep
      ja = istart + jjj
*
*     "transverse" loop
*     -----------------
      do 215 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      j = 0
*
*  loop across transforms
*  ----------------------
cdir$ ivdep, shortloop
      do 210 l = 1 , nvex
      t0 = a(ja+j) - a(jb+j)
      a(ja+j) = a(ja+j) + a(jb+j)
      a(jb+j) = t0
      u0 = b(ja+j) - b(jb+j)
      b(ja+j) = b(ja+j) + b(jb+j)
      b(jb+j) = u0
      j = j + jump
  210 continue
*-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  215 continue
  220 continue
*
*  finished if n2=2
*  ----------------
      if (n2.eq.2) go to 490
*
      kk = 2 * la
*
*  loop on nonzero k
*  -----------------
      do 260 k = ink , jstep - ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
*
*  loop along transforms
*  ---------------------
      do 250 jjj = k , (n-1)*inc , 2*jstep
      ja = istart + jjj
*
*     "transverse" loop
*     -----------------
      do 245 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      j = 0
*
*  loop across transforms
*  ----------------------
      if (kk.eq.n2/2) then
cdir$ ivdep, shortloop
      do 230 l = 1 , nvex
      t0 = ss * ( a(ja+j) - a(jb+j) )
      a(ja+j) = a(ja+j) + a(jb+j)
      a(jb+j) = ss * ( b(jb+j) - b(ja+j) )
      b(ja+j) = b(ja+j) + b(jb+j)
      b(jb+j) = t0
      j = j + jump
  230 continue
*
      else
*
cdir$ ivdep, shortloop
      do 240 l = 1 , nvex
      t0 = a(ja+j) - a(jb+j)
      a(ja+j) = a(ja+j) + a(jb+j)
      u0 = b(ja+j) - b(jb+j)
      b(ja+j) = b(ja+j) + b(jb+j)
      a(jb+j) = co1*t0 - si1*u0
      b(jb+j) = si1*t0 + co1*u0
      j = j + jump
  240 continue
*
      endif
*
*-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  245 continue
  250 continue
*-----(end of loop along transforms)
      kk = kk + 2 * la
  260 continue
*-----(end of loop on nonzero k)
*-----(end of radix-2 pass)
*
      la = 2 * la
      go to 400
*
*  central radix-8 pass
*  --------------------
  300 continue
      if (m8.eq.0) go to 400
      jstep = (n*inc) / (8*la)
      jstepl = jstep - ninc
      mu = mod(inq,8)
      if (isign.eq.-1) mu = 8 - mu
      c1 = 1.0
      if (mu.eq.3.or.mu.eq.7) c1 = -1.0
      c2 = sqrt(0.5)
      if (mu.eq.3.or.mu.eq.5) c2 = -c2
      c3 = c1 * c2
*
*  stage 1
*  -------
      do 320 k = 0 , jstep - ink , ink
      do 315 jjj = k , (n-1)*inc , 8*jstep
      ja = istart + jjj
*
*     "transverse" loop
*     -----------------
      do 312 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      j = 0
cdir$ ivdep, shortloop
      do 310 l = 1 , nvex
      t0 = a(ja+j) - a(je+j)
      a(ja+j) = a(ja+j) + a(je+j)
      t1 = c1 * ( a(jc+j) - a(jg+j) )
      a(je+j) = a(jc+j) + a(jg+j)
      t2 = a(jb+j) - a(jf+j)
      a(jc+j) = a(jb+j) + a(jf+j)
      t3 = a(jd+j) - a(jh+j)
      a(jg+j) = a(jd+j) + a(jh+j)
      a(jb+j) = t0
      a(jf+j) = t1
      a(jd+j) = c2 * ( t2 - t3 )
      a(jh+j) = c3 * ( t2 + t3 )
      u0 = b(ja+j) - b(je+j)
      b(ja+j) = b(ja+j) + b(je+j)
      u1 = c1 * ( b(jc+j) - b(jg+j) )
      b(je+j) = b(jc+j) + b(jg+j)
      u2 = b(jb+j) - b(jf+j)
      b(jc+j) = b(jb+j) + b(jf+j)
      u3 = b(jd+j) - b(jh+j)
      b(jg+j) = b(jd+j) + b(jh+j)
      b(jb+j) = u0
      b(jf+j) = u1
      b(jd+j) = c2 * ( u2 - u3 )
      b(jh+j) = c3 * ( u2 + u3 )
      j = j + jump
  310 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  312 continue
  315 continue
  320 continue
*
*  stage 2
*  -------
*
*  k=0 (no twiddle factors)
*  ------------------------
      do 330 jjj = 0 , (n-1)*inc , 8*jstep
      ja = istart + jjj
*
*     "transverse" loop
*     -----------------
      do 328 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      j = 0
cdir$ ivdep, shortloop
      do 325 l = 1 , nvex
      t0 = a(ja+j) + a(je+j)
      t2 = a(ja+j) - a(je+j)
      t1 = a(jc+j) + a(jg+j)
      t3 = c1 * ( a(jc+j) - a(jg+j) )
      u0 = b(ja+j) + b(je+j)
      u2 = b(ja+j) - b(je+j)
      u1 = b(jc+j) + b(jg+j)
      u3 = c1 * ( b(jc+j) - b(jg+j ) )
      a(ja+j) = t0 + t1
      a(je+j) = t0 - t1
      b(ja+j) = u0 + u1
      b(je+j) = u0 - u1
      a(jc+j) = t2 - u3
      a(jg+j) = t2 + u3
      b(jc+j) = u2 + t3
      b(jg+j) = u2 - t3
      t0 = a(jb+j) + a(jd+j)
      t2 = a(jb+j) - a(jd+j)
      t1 = a(jf+j) - a(jh+j)
      t3 = a(jf+j) + a(jh+j)
      u0 = b(jb+j) + b(jd+j)
      u2 = b(jb+j) - b(jd+j)
      u1 = b(jf+j) - b(jh+j)
      u3 = b(jf+j) + b(jh+j)
      a(jb+j) = t0 - u3
      a(jh+j) = t0 + u3
      b(jb+j) = u0 + t3
      b(jh+j) = u0 - t3
      a(jd+j) = t2 + u1
      a(jf+j) = t2 - u1
      b(jd+j) = u2 - t1
      b(jf+j) = u2 + t1
      j = j + jump
  325 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  328 continue
  330 continue
*
      if (n2.eq.8) go to 490
*
*  loop on nonzero k
*  -----------------
      kk = 2 * la
*
      do 350 k = ink , jstep - ink , ink
*
      co1 = trigs(kk+1)
      si1 = s * trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s * trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s * trigs(3*kk+2)
      co4 = trigs(4*kk+1)
      si4 = s * trigs(4*kk+2)
      co5 = trigs(5*kk+1)
      si5 = s * trigs(5*kk+2)
      co6 = trigs(6*kk+1)
      si6 = s * trigs(6*kk+2)
      co7 = trigs(7*kk+1)
      si7 = s * trigs(7*kk+2)
*
      do 345 jjj = k , (n-1)*inc , 8*jstep
      ja = istart + jjj
*
*     "transverse" loop
*     -----------------
      do 342 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      j = 0
cdir$ ivdep, shortloop
      do 340 l = 1 , nvex
      t0 = a(ja+j) + a(je+j)
      t2 = a(ja+j) - a(je+j)
      t1 = a(jc+j) + a(jg+j)
      t3 = c1 * ( a(jc+j) - a(jg+j) )
      u0 = b(ja+j) + b(je+j)
      u2 = b(ja+j) - b(je+j)
      u1 = b(jc+j) + b(jg+j)
      u3 = c1 * ( b(jc+j) - b(jg+j ) )
      a(ja+j) = t0 + t1
      b(ja+j) = u0 + u1
      a(je+j) = co4*(t0-t1) - si4*(u0-u1)
      b(je+j) = si4*(t0-t1) + co4*(u0-u1)
      a(jc+j) = co2*(t2-u3) - si2*(u2+t3)
      b(jc+j) = si2*(t2-u3) + co2*(u2+t3)
      a(jg+j) = co6*(t2+u3) - si6*(u2-t3)
      b(jg+j) = si6*(t2+u3) + co6*(u2-t3)
      t0 = a(jb+j) + a(jd+j)
      t2 = a(jb+j) - a(jd+j)
      t1 = a(jf+j) - a(jh+j)
      t3 = a(jf+j) + a(jh+j)
      u0 = b(jb+j) + b(jd+j)
      u2 = b(jb+j) - b(jd+j)
      u1 = b(jf+j) - b(jh+j)
      u3 = b(jf+j) + b(jh+j)
      a(jb+j) = co1*(t0-u3) - si1*(u0+t3)
      b(jb+j) = si1*(t0-u3) + co1*(u0+t3)
      a(jh+j) = co7*(t0+u3) - si7*(u0-t3)
      b(jh+j) = si7*(t0+u3) + co7*(u0-t3)
      a(jd+j) = co3*(t2+u1) - si3*(u2-t1)
      b(jd+j) = si3*(t2+u1) + co3*(u2-t1)
      a(jf+j) = co5*(t2-u1) - si5*(u2+t1)
      b(jf+j) = si5*(t2-u1) + co5*(u2+t1)
      j = j + jump
  340 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  342 continue
  345 continue
      kk = kk + 2 * la
  350 continue
*
      la = 8 * la
*
*  loop on type ii radix-4 passes
*  ------------------------------
  400 continue
      mu = mod(inq,4)
      if (isign.eq.-1) mu = 4 - mu
      ss = 1.0
      if (mu.eq.3) ss = -1.0
*
      do 480 ipass = mh+1 , m
      jstep = (n*inc) / (4*la)
      jstepl = jstep - ninc
      laincl = la * ink - ninc
*
*  k=0 loop (no twiddle factors)
*  -----------------------------
      do 430 ll = 0 , (la-1)*ink , 4*jstep
*
      do 420 jjj = ll , (n-1)*inc , 4*la*ink
      ja = istart + jjj
*
*     "transverse" loop
*     -----------------
      do 415 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = ja + laincl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = je + laincl
      if (ji.lt.istart) ji = ji + ninc
      jj = ji + jstepl
      if (jj.lt.istart) jj = jj + ninc
      jk = jj + jstepl
      if (jk.lt.istart) jk = jk + ninc
      jl = jk + jstepl
      if (jl.lt.istart) jl = jl + ninc
      jm = ji + laincl
      if (jm.lt.istart) jm = jm + ninc
      jn = jm + jstepl
      if (jn.lt.istart) jn = jn + ninc
      jo = jn + jstepl
      if (jo.lt.istart) jo = jo + ninc
      jp = jo + jstepl
      if (jp.lt.istart) jp = jp + ninc
      j = 0
*
*  loop across transforms
*  ----------------------
cdir$ ivdep, shortloop
      do 410 l = 1 , nvex
      t0 = a(ja+j) + a(jc+j)
      t2 = a(ja+j) - a(jc+j)
      t1 = a(jb+j) + a(jd+j)
      t3 = ss * ( a(jb+j) - a(jd+j) )
      a(jc+j) = a(ji+j)
      u0 = b(ja+j) + b(jc+j)
      u2 = b(ja+j) - b(jc+j)
      u1 = b(jb+j) + b(jd+j)
      u3 = ss * ( b(jb+j) - b(jd+j) )
      a(jb+j) = a(je+j)
      a(ja+j) = t0 + t1
      a(ji+j) = t0 - t1
      b(ja+j) = u0 + u1
      b(jc+j) = u0 - u1
      b(jd+j) = b(jm+j)
      a(je+j) = t2 - u3
      a(jd+j) = t2 + u3
      b(jb+j) = u2 + t3
      b(jm+j) = u2 - t3
*----------------------
      t0 = a(jb+j) + a(jg+j)
      t2 = a(jb+j) - a(jg+j)
      t1 = a(jf+j) + a(jh+j)
      t3 = ss * ( a(jf+j) - a(jh+j) )
      a(jg+j) = a(jj+j)
      u0 = b(je+j) + b(jg+j)
      u2 = b(je+j) - b(jg+j)
      u1 = b(jf+j) + b(jh+j)
      u3 = ss * ( b(jf+j) - b(jh+j) )
      b(je+j) = b(jb+j)
      a(jb+j) = t0 + t1
      a(jj+j) = t0 - t1
      b(jg+j) = b(jj+j)
      b(jb+j) = u0 + u1
      b(jj+j) = u0 - u1
      a(jf+j) = t2 - u3
      a(jh+j) = t2 + u3
      b(jf+j) = u2 + t3
      b(jh+j) = u2 - t3
*----------------------
      t0 = a(jc+j) + a(jk+j)
      t2 = a(jc+j) - a(jk+j)
      t1 = a(jg+j) + a(jl+j)
      t3 = ss * ( a(jg+j) - a(jl+j) )
      u0 = b(ji+j) + b(jk+j)
      u2 = b(ji+j) - b(jk+j)
      a(jl+j) = a(jo+j)
      u1 = b(jg+j) + b(jl+j)
      u3 = ss * ( b(jg+j) - b(jl+j) )
      b(ji+j) = b(jc+j)
      a(jc+j) = t0 + t1
      a(jk+j) = t0 - t1
      b(jl+j) = b(jo+j)
      b(jc+j) = u0 + u1
      b(jk+j) = u0 - u1
      a(jg+j) = t2 - u3
      a(jo+j) = t2 + u3
      b(jg+j) = u2 + t3
      b(jo+j) = u2 - t3
*----------------------
      t0 = a(jm+j) + a(jl+j)
      t2 = a(jm+j) - a(jl+j)
      t1 = a(jn+j) + a(jp+j)
      t3 = ss * ( a(jn+j) - a(jp+j) )
      a(jm+j) = a(jd+j)
      u0 = b(jd+j) + b(jl+j)
      u2 = b(jd+j) - b(jl+j)
      u1 = b(jn+j) + b(jp+j)
      u3 = ss * ( b(jn+j) - b(jp+j) )
      a(jn+j) = a(jh+j)
      a(jd+j) = t0 + t1
      a(jl+j) = t0 - t1
      b(jd+j) = u0 + u1
      b(jl+j) = u0 - u1
      b(jn+j) = b(jh+j)
      a(jh+j) = t2 - u3
      a(jp+j) = t2 + u3
      b(jh+j) = u2 + t3
      b(jp+j) = u2 - t3
      j = j + jump
  410 continue
*-----( end of loop across transforms )
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  415 continue
  420 continue
  430 continue
*-----( end of double loop for k=0 )
*
*  finished if last pass
*  ---------------------
      if (ipass.eq.m) go to 490
*
      kk = 2*la
*
*     loop on nonzero k
*     -----------------
      do 470 k = ink , jstep-ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s*trigs(3*kk+2)
*
*  double loop along first transform in block
*  ------------------------------------------
      do 460 ll = k , (la-1)*ink , 4*jstep
*
      do 450 jjj = ll , (n-1)*inc , 4*la*ink
      ja = istart + jjj
*
*     "transverse" loop
*     -----------------
      do 445 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = ja + laincl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = je + laincl
      if (ji.lt.istart) ji = ji + ninc
      jj = ji + jstepl
      if (jj.lt.istart) jj = jj + ninc
      jk = jj + jstepl
      if (jk.lt.istart) jk = jk + ninc
      jl = jk + jstepl
      if (jl.lt.istart) jl = jl + ninc
      jm = ji + laincl
      if (jm.lt.istart) jm = jm + ninc
      jn = jm + jstepl
      if (jn.lt.istart) jn = jn + ninc
      jo = jn + jstepl
      if (jo.lt.istart) jo = jo + ninc
      jp = jo + jstepl
      if (jp.lt.istart) jp = jp + ninc
      j = 0
*
*  loop across transforms
*  ----------------------
cdir$ ivdep, shortloop
      do 440 l = 1 , nvex
      t0 = a(ja+j) + a(jc+j)
      t2 = a(ja+j) - a(jc+j)
      t1 = a(jb+j) + a(jd+j)
      t3 = ss * ( a(jb+j) - a(jd+j) )
      a(jc+j) = a(ji+j)
      u0 = b(ja+j) + b(jc+j)
      u2 = b(ja+j) - b(jc+j)
      u1 = b(jb+j) + b(jd+j)
      u3 = ss * ( b(jb+j) - b(jd+j) )
      a(jb+j) = a(je+j)
      a(ja+j) = t0 + t1
      b(ja+j) = u0 + u1
      a(je+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jb+j) = si1*(t2-u3) + co1*(u2+t3)
      b(jd+j) = b(jm+j)
      a(ji+j) = co2*(t0-t1) - si2*(u0-u1)
      b(jc+j) = si2*(t0-t1) + co2*(u0-u1)
      a(jd+j) = co3*(t2+u3) - si3*(u2-t3)
      b(jm+j) = si3*(t2+u3) + co3*(u2-t3)
*----------------------------------------
      t0 = a(jb+j) + a(jg+j)
      t2 = a(jb+j) - a(jg+j)
      t1 = a(jf+j) + a(jh+j)
      t3 = ss * ( a(jf+j) - a(jh+j) )
      a(jg+j) = a(jj+j)
      u0 = b(je+j) + b(jg+j)
      u2 = b(je+j) - b(jg+j)
      u1 = b(jf+j) + b(jh+j)
      u3 = ss * ( b(jf+j) - b(jh+j) )
      b(je+j) = b(jb+j)
      a(jb+j) = t0 + t1
      b(jb+j) = u0 + u1
      b(jg+j) = b(jj+j)
      a(jf+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jf+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jj+j) = co2*(t0-t1) - si2*(u0-u1)
      b(jj+j) = si2*(t0-t1) + co2*(u0-u1)
      a(jh+j) = co3*(t2+u3) - si3*(u2-t3)
      b(jh+j) = si3*(t2+u3) + co3*(u2-t3)
*----------------------------------------
      t0 = a(jc+j) + a(jk+j)
      t2 = a(jc+j) - a(jk+j)
      t1 = a(jg+j) + a(jl+j)
      t3 = ss * ( a(jg+j) - a(jl+j) )
      u0 = b(ji+j) + b(jk+j)
      u2 = b(ji+j) - b(jk+j)
      a(jl+j) = a(jo+j)
      u1 = b(jg+j) + b(jl+j)
      u3 = ss * ( b(jg+j) - b(jl+j) )
      b(ji+j) = b(jc+j)
      a(jc+j) = t0 + t1
      b(jc+j) = u0 + u1
      b(jl+j) = b(jo+j)
      a(jg+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jg+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jk+j) = co2*(t0-t1) - si2*(u0-u1)
      b(jk+j) = si2*(t0-t1) + co2*(u0-u1)
      a(jo+j) = co3*(t2+u3) - si3*(u2-t3)
      b(jo+j) = si3*(t2+u3) + co3*(u2-t3)
*----------------------------------------
      t0 = a(jm+j) + a(jl+j)
      t2 = a(jm+j) - a(jl+j)
      t1 = a(jn+j) + a(jp+j)
      t3 = ss * ( a(jn+j) - a(jp+j) )
      a(jm+j) = a(jd+j)
      u0 = b(jd+j) + b(jl+j)
      u2 = b(jd+j) - b(jl+j)
      a(jn+j) = a(jh+j)
      u1 = b(jn+j) + b(jp+j)
      u3 = ss * ( b(jn+j) - b(jp+j) )
      b(jn+j) = b(jh+j)
      a(jd+j) = t0 + t1
      b(jd+j) = u0 + u1
      a(jh+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jh+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jl+j) = co2*(t0-t1) - si2*(u0-u1)
      b(jl+j) = si2*(t0-t1) + co2*(u0-u1)
      a(jp+j) = co3*(t2+u3) - si3*(u2-t3)
      b(jp+j) = si3*(t2+u3) + co3*(u2-t3)
      j = j + jump
  440 continue
*-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  445 continue
  450 continue
  460 continue
*-----( end of double loop for this k )
      kk = kk + 2*la
  470 continue
*-----( end of loop over values of k )
      la = 4*la
  480 continue
*-----( end of loop on type ii radix-4 passes )
*-----( nvex transforms completed)
  490 continue
      istart = istart + nvex * jump
  500 continue
*-----( end of loop on blocks of transforms )
*
      return
      end
*     fortran version of *gpfa3* -
*     radix-3 section of self-sorting, in-place
*        generalized pfa
*
*-------------------------------------------------------------------
*
      subroutine gpfa3f(a,b,trigs,inc,jump,n,mm,lot,isign)
      dimension a(*), b(*), trigs(*)
      data sin60/0.866025403784437/
*
      n3 = 3**mm
      inq = n/n3
      jstepx = (n3-n) * inc
      ninc = n * inc
      ink = inc * inq
      mu = mod(inq,3)
      if (isign.eq.-1) mu = 3-mu
      m = mm
      mh = (m+1)/2
      s = float(isign)
      c1 = sin60
      if (mu.eq.2) c1 = -c1
*
      nblox = 1 + (lot-1)/64
      left = lot
      istart = 1
*
*  loop on blocks of 64 transforms
*  -------------------------------
      do 500 nb = 1 , nblox
*
      if (left.le.64) then
         nvex = left
      else if (left.lt.128) then
         nvex = left/2
      else
         nvex = 64
      endif
      left = left - nvex
*
      la = 1
*
*  loop on type i radix-3 passes
*  -----------------------------
      do 160 ipass = 1 , mh
      jstep = (n*inc) / (3*la)
      jstepl = jstep - ninc
*
*  k = 0 loop (no twiddle factors)
*  -------------------------------
      do 120 jjj = 0 , (n-1)*inc , 3*jstep
      ja = istart + jjj
*
*  "transverse" loop
*  -----------------
      do 115 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      j = 0
*
*  loop across transforms
*  ----------------------
cdir$ ivdep, shortloop
      do 110 l = 1 , nvex
      t1 = a(jb+j) + a(jc+j)
      t2 = a(ja+j) - 0.5 * t1
      t3 = c1 * ( a(jb+j) - a(jc+j) )
      u1 = b(jb+j) + b(jc+j)
      u2 = b(ja+j) - 0.5 * u1
      u3 = c1 * ( b(jb+j) - b(jc+j) )
      a(ja+j) = a(ja+j) + t1
      b(ja+j) = b(ja+j) + u1
      a(jb+j) = t2 - u3
      b(jb+j) = u2 + t3
      a(jc+j) = t2 + u3
      b(jc+j) = u2 - t3
      j = j + jump
  110 continue
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  115 continue
  120 continue
*
*  finished if n3 = 3
*  ------------------
      if (n3.eq.3) go to 490
      kk = 2 * la
*
*  loop on nonzero k
*  -----------------
      do 150 k = ink , jstep-ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
*
*  loop along transform
*  --------------------
      do 140 jjj = k , (n-1)*inc , 3*jstep
      ja = istart + jjj
*
*  "transverse" loop
*  -----------------
      do 135 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      j = 0
*
*  loop across transforms
*  ----------------------
cdir$ ivdep,shortloop
      do 130 l = 1 , nvex
      t1 = a(jb+j) + a(jc+j)
      t2 = a(ja+j) - 0.5 * t1
      t3 = c1 * ( a(jb+j) - a(jc+j) )
      u1 = b(jb+j) + b(jc+j)
      u2 = b(ja+j) - 0.5 * u1
      u3 = c1 * ( b(jb+j) - b(jc+j) )
      a(ja+j) = a(ja+j) + t1
      b(ja+j) = b(ja+j) + u1
      a(jb+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jb+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jc+j) = co2*(t2+u3) - si2*(u2-t3)
      b(jc+j) = si2*(t2+u3) + co2*(u2-t3)
      j = j + jump
  130 continue
*-----( end of loop across transforms )
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  135 continue
  140 continue
*-----( end of loop along transforms )
      kk = kk + 2*la
  150 continue
*-----( end of loop on nonzero k )
      la = 3*la
  160 continue
*-----( end of loop on type i radix-3 passes)
*
*  loop on type ii radix-3 passes
*  ------------------------------
  400 continue
*
      do 480 ipass = mh+1 , m
      jstep = (n*inc) / (3*la)
      jstepl = jstep - ninc
      laincl = la*ink - ninc
*
*  k=0 loop (no twiddle factors)
*  -----------------------------
      do 430 ll = 0 , (la-1)*ink , 3*jstep
*
      do 420 jjj = ll , (n-1)*inc , 3*la*ink
      ja = istart + jjj
*
*  "transverse" loop
*  -----------------
      do 415 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = ja + laincl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jd + laincl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = jh + jstepl
      if (ji.lt.istart) ji = ji + ninc
      j = 0
*
*  loop across transforms
*  ----------------------
cdir$ ivdep, shortloop
      do 410 l = 1 , nvex
      t1 = a(jb+j) + a(jc+j)
      t2 = a(ja+j) - 0.5 * t1
      t3 = c1 * ( a(jb+j) - a(jc+j) )
      a(jb+j) = a(jd+j)
      u1 = b(jb+j) + b(jc+j)
      u2 = b(ja+j) - 0.5 * u1
      u3 = c1 * ( b(jb+j) - b(jc+j) )
      b(jb+j) = b(jd+j)
      a(ja+j) = a(ja+j) + t1
      b(ja+j) = b(ja+j) + u1
      a(jd+j) = t2 - u3
      b(jd+j) = u2 + t3
      a(jc+j) = t2 + u3
      b(jc+j) = u2 - t3
*----------------------
      t1 = a(je+j) + a(jf+j)
      t2 = a(jb+j) - 0.5 * t1
      t3 = c1 * ( a(je+j) - a(jf+j) )
      a(jf+j) = a(jh+j)
      u1 = b(je+j) + b(jf+j)
      u2 = b(jb+j) - 0.5 * u1
      u3 = c1 * ( b(je+j) - b(jf+j) )
      b(jf+j) = b(jh+j)
      a(jb+j) = a(jb+j) + t1
      b(jb+j) = b(jb+j) + u1
      a(je+j) = t2 - u3
      b(je+j) = u2 + t3
      a(jh+j) = t2 + u3
      b(jh+j) = u2 - t3
*----------------------
      t1 = a(jf+j) + a(ji+j)
      t2 = a(jg+j) - 0.5 * t1
      t3 = c1 * ( a(jf+j) - a(ji+j) )
      t1 = a(jg+j) + t1
      a(jg+j) = a(jc+j)
      u1 = b(jf+j) + b(ji+j)
      u2 = b(jg+j) - 0.5 * u1
      u3 = c1 * ( b(jf+j) - b(ji+j) )
      u1 = b(jg+j) + u1
      b(jg+j) = b(jc+j)
      a(jc+j) = t1
      b(jc+j) = u1
      a(jf+j) = t2 - u3
      b(jf+j) = u2 + t3
      a(ji+j) = t2 + u3
      b(ji+j) = u2 - t3
      j = j + jump
  410 continue
*-----( end of loop across transforms )
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  415 continue
  420 continue
  430 continue
*-----( end of double loop for k=0 )
*
*  finished if last pass
*  ---------------------
      if (ipass.eq.m) go to 490
*
      kk = 2*la
*
*     loop on nonzero k
*     -----------------
      do 470 k = ink , jstep-ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
*
*  double loop along first transform in block
*  ------------------------------------------
      do 460 ll = k , (la-1)*ink , 3*jstep
*
      do 450 jjj = ll , (n-1)*inc , 3*la*ink
      ja = istart + jjj
*
*  "transverse" loop
*  -----------------
      do 445 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = ja + laincl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = je + jstepl
      if (jf.lt.istart) jf = jf + ninc
      jg = jd + laincl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = jh + jstepl
      if (ji.lt.istart) ji = ji + ninc
      j = 0
*
*  loop across transforms
*  ----------------------
cdir$ ivdep, shortloop
      do 440 l = 1 , nvex
      t1 = a(jb+j) + a(jc+j)
      t2 = a(ja+j) - 0.5 * t1
      t3 = c1 * ( a(jb+j) - a(jc+j) )
      a(jb+j) = a(jd+j)
      u1 = b(jb+j) + b(jc+j)
      u2 = b(ja+j) - 0.5 * u1
      u3 = c1 * ( b(jb+j) - b(jc+j) )
      b(jb+j) = b(jd+j)
      a(ja+j) = a(ja+j) + t1
      b(ja+j) = b(ja+j) + u1
      a(jd+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jd+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jc+j) = co2*(t2+u3) - si2*(u2-t3)
      b(jc+j) = si2*(t2+u3) + co2*(u2-t3)
*----------------------
      t1 = a(je+j) + a(jf+j)
      t2 = a(jb+j) - 0.5 * t1
      t3 = c1 * ( a(je+j) - a(jf+j) )
      a(jf+j) = a(jh+j)
      u1 = b(je+j) + b(jf+j)
      u2 = b(jb+j) - 0.5 * u1
      u3 = c1 * ( b(je+j) - b(jf+j) )
      b(jf+j) = b(jh+j)
      a(jb+j) = a(jb+j) + t1
      b(jb+j) = b(jb+j) + u1
      a(je+j) = co1*(t2-u3) - si1*(u2+t3)
      b(je+j) = si1*(t2-u3) + co1*(u2+t3)
      a(jh+j) = co2*(t2+u3) - si2*(u2-t3)
      b(jh+j) = si2*(t2+u3) + co2*(u2-t3)
*----------------------
      t1 = a(jf+j) + a(ji+j)
      t2 = a(jg+j) - 0.5 * t1
      t3 = c1 * ( a(jf+j) - a(ji+j) )
      t1 = a(jg+j) + t1
      a(jg+j) = a(jc+j)
      u1 = b(jf+j) + b(ji+j)
      u2 = b(jg+j) - 0.5 * u1
      u3 = c1 * ( b(jf+j) - b(ji+j) )
      u1 = b(jg+j) + u1
      b(jg+j) = b(jc+j)
      a(jc+j) = t1
      b(jc+j) = u1
      a(jf+j) = co1*(t2-u3) - si1*(u2+t3)
      b(jf+j) = si1*(t2-u3) + co1*(u2+t3)
      a(ji+j) = co2*(t2+u3) - si2*(u2-t3)
      b(ji+j) = si2*(t2+u3) + co2*(u2-t3)
      j = j + jump
  440 continue
*-----(end of loop across transforms)
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  445 continue
  450 continue
  460 continue
*-----( end of double loop for this k )
      kk = kk + 2*la
  470 continue
*-----( end of loop over values of k )
      la = 3*la
  480 continue
*-----( end of loop on type ii radix-3 passes )
*-----( nvex transforms completed)
  490 continue
      istart = istart + nvex * jump
  500 continue
*-----( end of loop on blocks of transforms )
*
      return
      end
*     fortran version of *gpfa5* -
*     radix-5 section of self-sorting, in-place,
*        generalized pfa
*
*-------------------------------------------------------------------
*
      subroutine gpfa5f(a,b,trigs,inc,jump,n,mm,lot,isign)
      dimension a(*), b(*), trigs(*)
      data sin36/0.587785252292473/, sin72/0.951056516295154/,
     *      qrt5/0.559016994374947/
*
      n5 = 5 ** mm
      inq = n / n5
      jstepx = (n5-n) * inc
      ninc = n * inc
      ink = inc * inq
      mu = mod(inq,5)
      if (isign.eq.-1) mu = 5 - mu
*
      m = mm
      mh = (m+1)/2
      s = float(isign)
      c1 = qrt5
      c2 = sin72
      c3 = sin36
      if (mu.eq.2.or.mu.eq.3) then
         c1 = -c1
         c2 = sin36
         c3 = sin72
      endif
      if (mu.eq.3.or.mu.eq.4) c2 = -c2
      if (mu.eq.2.or.mu.eq.4) c3 = -c3
*
      nblox = 1 + (lot-1)/64
      left = lot
      istart = 1
*
*  loop on blocks of 64 transforms
*  -------------------------------
      do 500 nb = 1 , nblox
*
      if (left.le.64) then
         nvex = left
      else if (left.lt.128) then
         nvex = left/2
      else
         nvex = 64
      endif
      left = left - nvex
*
      la = 1
*
*  loop on type i radix-5 passes
*  -----------------------------
      do 160 ipass = 1 , mh
      jstep = (n*inc) / (5*la)
      jstepl = jstep - ninc
      kk = 0
*
*  loop on k
*  ---------
      do 150 k = 0 , jstep-ink , ink
*
      if (k.gt.0) then
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s*trigs(3*kk+2)
      co4 = trigs(4*kk+1)
      si4 = s*trigs(4*kk+2)
      endif
*
*  loop along transform
*  --------------------
      do 140 jjj = k , (n-1)*inc , 5*jstep
      ja = istart + jjj
*
*     "transverse" loop
*     -----------------
      do 135 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      j = 0
*
*  loop across transforms
*  ----------------------
      if (k.eq.0) then
*
cdir$ ivdep, shortloop
      do 110 l = 1 , nvex
      t1 = a(jb+j) + a(je+j)
      t2 = a(jc+j) + a(jd+j)
      t3 = a(jb+j) - a(je+j)
      t4 = a(jc+j) - a(jd+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ja+j) - 0.25 * t5
      a(ja+j) = a(ja+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jb+j) + b(je+j)
      u2 = b(jc+j) + b(jd+j)
      u3 = b(jb+j) - b(je+j)
      u4 = b(jc+j) - b(jd+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ja+j) - 0.25 * u5
      b(ja+j) = b(ja+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jb+j) = t8 - u11
      b(jb+j) = u8 + t11
      a(je+j) = t8 + u11
      b(je+j) = u8 - t11
      a(jc+j) = t9 - u10
      b(jc+j) = u9 + t10
      a(jd+j) = t9 + u10
      b(jd+j) = u9 - t10
      j = j + jump
  110 continue
*
      else
*
cdir$ ivdep,shortloop
      do 130 l = 1 , nvex
      t1 = a(jb+j) + a(je+j)
      t2 = a(jc+j) + a(jd+j)
      t3 = a(jb+j) - a(je+j)
      t4 = a(jc+j) - a(jd+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ja+j) - 0.25 * t5
      a(ja+j) = a(ja+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jb+j) + b(je+j)
      u2 = b(jc+j) + b(jd+j)
      u3 = b(jb+j) - b(je+j)
      u4 = b(jc+j) - b(jd+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ja+j) - 0.25 * u5
      b(ja+j) = b(ja+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jb+j) = co1*(t8-u11) - si1*(u8+t11)
      b(jb+j) = si1*(t8-u11) + co1*(u8+t11)
      a(je+j) = co4*(t8+u11) - si4*(u8-t11)
      b(je+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jc+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jc+j) = si2*(t9-u10) + co2*(u9+t10)
      a(jd+j) = co3*(t9+u10) - si3*(u9-t10)
      b(jd+j) = si3*(t9+u10) + co3*(u9-t10)
      j = j + jump
  130 continue
*
      endif
*
*-----( end of loop across transforms )
*
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  135 continue
  140 continue
*-----( end of loop along transforms )
      kk = kk + 2*la
  150 continue
*-----( end of loop on nonzero k )
      la = 5*la
  160 continue
*-----( end of loop on type i radix-5 passes)
*
      if (n.eq.5) go to 490
*
*  loop on type ii radix-5 passes
*  ------------------------------
  400 continue
*
      do 480 ipass = mh+1 , m
      jstep = (n*inc) / (5*la)
      jstepl = jstep - ninc
      laincl = la * ink - ninc
      kk = 0
*
*     loop on k
*     ---------
      do 470 k = 0 , jstep-ink , ink
*
      if (k.gt.0) then
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s*trigs(3*kk+2)
      co4 = trigs(4*kk+1)
      si4 = s*trigs(4*kk+2)
      endif
*
*  double loop along first transform in block
*  ------------------------------------------
      do 460 ll = k , (la-1)*ink , 5*jstep
*
      do 450 jjj = ll , (n-1)*inc , 5*la*ink
      ja = istart + jjj
*
*     "transverse" loop
*     -----------------
      do 445 nu = 1 , inq
      jb = ja + jstepl
      if (jb.lt.istart) jb = jb + ninc
      jc = jb + jstepl
      if (jc.lt.istart) jc = jc + ninc
      jd = jc + jstepl
      if (jd.lt.istart) jd = jd + ninc
      je = jd + jstepl
      if (je.lt.istart) je = je + ninc
      jf = ja + laincl
      if (jf.lt.istart) jf = jf + ninc
      jg = jf + jstepl
      if (jg.lt.istart) jg = jg + ninc
      jh = jg + jstepl
      if (jh.lt.istart) jh = jh + ninc
      ji = jh + jstepl
      if (ji.lt.istart) ji = ji + ninc
      jj = ji + jstepl
      if (jj.lt.istart) jj = jj + ninc
      jk = jf + laincl
      if (jk.lt.istart) jk = jk + ninc
      jl = jk + jstepl
      if (jl.lt.istart) jl = jl + ninc
      jm = jl + jstepl
      if (jm.lt.istart) jm = jm + ninc
      jn = jm + jstepl
      if (jn.lt.istart) jn = jn + ninc
      jo = jn + jstepl
      if (jo.lt.istart) jo = jo + ninc
      jp = jk + laincl
      if (jp.lt.istart) jp = jp + ninc
      jq = jp + jstepl
      if (jq.lt.istart) jq = jq + ninc
      jr = jq + jstepl
      if (jr.lt.istart) jr = jr + ninc
      js = jr + jstepl
      if (js.lt.istart) js = js + ninc
      jt = js + jstepl
      if (jt.lt.istart) jt = jt + ninc
      ju = jp + laincl
      if (ju.lt.istart) ju = ju + ninc
      jv = ju + jstepl
      if (jv.lt.istart) jv = jv + ninc
      jw = jv + jstepl
      if (jw.lt.istart) jw = jw + ninc
      jx = jw + jstepl
      if (jx.lt.istart) jx = jx + ninc
      jy = jx + jstepl
      if (jy.lt.istart) jy = jy + ninc
      j = 0
*
*  loop across transforms
*  ----------------------
      if (k.eq.0) then
*
cdir$ ivdep, shortloop
      do 410 l = 1 , nvex
      t1 = a(jb+j) + a(je+j)
      t2 = a(jc+j) + a(jd+j)
      t3 = a(jb+j) - a(je+j)
      t4 = a(jc+j) - a(jd+j)
      a(jb+j) = a(jf+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ja+j) - 0.25 * t5
      a(ja+j) = a(ja+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jc+j) = a(jk+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jb+j) + b(je+j)
      u2 = b(jc+j) + b(jd+j)
      u3 = b(jb+j) - b(je+j)
      u4 = b(jc+j) - b(jd+j)
      b(jb+j) = b(jf+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ja+j) - 0.25 * u5
      b(ja+j) = b(ja+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jc+j) = b(jk+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jf+j) = t8 - u11
      b(jf+j) = u8 + t11
      a(je+j) = t8 + u11
      b(je+j) = u8 - t11
      a(jk+j) = t9 - u10
      b(jk+j) = u9 + t10
      a(jd+j) = t9 + u10
      b(jd+j) = u9 - t10
*----------------------
      t1 = a(jg+j) + a(jj+j)
      t2 = a(jh+j) + a(ji+j)
      t3 = a(jg+j) - a(jj+j)
      t4 = a(jh+j) - a(ji+j)
      a(jh+j) = a(jl+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jb+j) - 0.25 * t5
      a(jb+j) = a(jb+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(ji+j) = a(jq+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jg+j) + b(jj+j)
      u2 = b(jh+j) + b(ji+j)
      u3 = b(jg+j) - b(jj+j)
      u4 = b(jh+j) - b(ji+j)
      b(jh+j) = b(jl+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jb+j) - 0.25 * u5
      b(jb+j) = b(jb+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(ji+j) = b(jq+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jg+j) = t8 - u11
      b(jg+j) = u8 + t11
      a(jj+j) = t8 + u11
      b(jj+j) = u8 - t11
      a(jl+j) = t9 - u10
      b(jl+j) = u9 + t10
      a(jq+j) = t9 + u10
      b(jq+j) = u9 - t10
*----------------------
      t1 = a(jh+j) + a(jo+j)
      t2 = a(jm+j) + a(jn+j)
      t3 = a(jh+j) - a(jo+j)
      t4 = a(jm+j) - a(jn+j)
      a(jn+j) = a(jr+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jc+j) - 0.25 * t5
      a(jc+j) = a(jc+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jo+j) = a(jw+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jh+j) + b(jo+j)
      u2 = b(jm+j) + b(jn+j)
      u3 = b(jh+j) - b(jo+j)
      u4 = b(jm+j) - b(jn+j)
      b(jn+j) = b(jr+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jc+j) - 0.25 * u5
      b(jc+j) = b(jc+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jo+j) = b(jw+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jh+j) = t8 - u11
      b(jh+j) = u8 + t11
      a(jw+j) = t8 + u11
      b(jw+j) = u8 - t11
      a(jm+j) = t9 - u10
      b(jm+j) = u9 + t10
      a(jr+j) = t9 + u10
      b(jr+j) = u9 - t10
*----------------------
      t1 = a(ji+j) + a(jt+j)
      t2 = a(jn+j) + a(js+j)
      t3 = a(ji+j) - a(jt+j)
      t4 = a(jn+j) - a(js+j)
      a(jt+j) = a(jx+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jp+j) - 0.25 * t5
      ax = a(jp+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jp+j) = a(jd+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      a(jd+j) = ax
      u1 = b(ji+j) + b(jt+j)
      u2 = b(jn+j) + b(js+j)
      u3 = b(ji+j) - b(jt+j)
      u4 = b(jn+j) - b(js+j)
      b(jt+j) = b(jx+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jp+j) - 0.25 * u5
      bx = b(jp+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jp+j) = b(jd+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      b(jd+j) = bx
      a(ji+j) = t8 - u11
      b(ji+j) = u8 + t11
      a(jx+j) = t8 + u11
      b(jx+j) = u8 - t11
      a(jn+j) = t9 - u10
      b(jn+j) = u9 + t10
      a(js+j) = t9 + u10
      b(js+j) = u9 - t10
*----------------------
      t1 = a(jv+j) + a(jy+j)
      t2 = a(jo+j) + a(jt+j)
      t3 = a(jv+j) - a(jy+j)
      t4 = a(jo+j) - a(jt+j)
      a(jv+j) = a(jj+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ju+j) - 0.25 * t5
      ax = a(ju+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(ju+j) = a(je+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      a(je+j) = ax
      u1 = b(jv+j) + b(jy+j)
      u2 = b(jo+j) + b(jt+j)
      u3 = b(jv+j) - b(jy+j)
      u4 = b(jo+j) - b(jt+j)
      b(jv+j) = b(jj+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ju+j) - 0.25 * u5
      bx = b(ju+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(ju+j) = b(je+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      b(je+j) = bx
      a(jj+j) = t8 - u11
      b(jj+j) = u8 + t11
      a(jy+j) = t8 + u11
      b(jy+j) = u8 - t11
      a(jo+j) = t9 - u10
      b(jo+j) = u9 + t10
      a(jt+j) = t9 + u10
      b(jt+j) = u9 - t10
      j = j + jump
  410 continue
*
      else
*
cdir$ ivdep, shortloop
      do 440 l = 1 , nvex
      t1 = a(jb+j) + a(je+j)
      t2 = a(jc+j) + a(jd+j)
      t3 = a(jb+j) - a(je+j)
      t4 = a(jc+j) - a(jd+j)
      a(jb+j) = a(jf+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ja+j) - 0.25 * t5
      a(ja+j) = a(ja+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jc+j) = a(jk+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jb+j) + b(je+j)
      u2 = b(jc+j) + b(jd+j)
      u3 = b(jb+j) - b(je+j)
      u4 = b(jc+j) - b(jd+j)
      b(jb+j) = b(jf+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ja+j) - 0.25 * u5
      b(ja+j) = b(ja+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jc+j) = b(jk+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jf+j) = co1*(t8-u11) - si1*(u8+t11)
      b(jf+j) = si1*(t8-u11) + co1*(u8+t11)
      a(je+j) = co4*(t8+u11) - si4*(u8-t11)
      b(je+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jk+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jk+j) = si2*(t9-u10) + co2*(u9+t10)
      a(jd+j) = co3*(t9+u10) - si3*(u9-t10)
      b(jd+j) = si3*(t9+u10) + co3*(u9-t10)
*----------------------
      t1 = a(jg+j) + a(jj+j)
      t2 = a(jh+j) + a(ji+j)
      t3 = a(jg+j) - a(jj+j)
      t4 = a(jh+j) - a(ji+j)
      a(jh+j) = a(jl+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jb+j) - 0.25 * t5
      a(jb+j) = a(jb+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(ji+j) = a(jq+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jg+j) + b(jj+j)
      u2 = b(jh+j) + b(ji+j)
      u3 = b(jg+j) - b(jj+j)
      u4 = b(jh+j) - b(ji+j)
      b(jh+j) = b(jl+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jb+j) - 0.25 * u5
      b(jb+j) = b(jb+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(ji+j) = b(jq+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jg+j) = co1*(t8-u11) - si1*(u8+t11)
      b(jg+j) = si1*(t8-u11) + co1*(u8+t11)
      a(jj+j) = co4*(t8+u11) - si4*(u8-t11)
      b(jj+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jl+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jl+j) = si2*(t9-u10) + co2*(u9+t10)
      a(jq+j) = co3*(t9+u10) - si3*(u9-t10)
      b(jq+j) = si3*(t9+u10) + co3*(u9-t10)
*----------------------
      t1 = a(jh+j) + a(jo+j)
      t2 = a(jm+j) + a(jn+j)
      t3 = a(jh+j) - a(jo+j)
      t4 = a(jm+j) - a(jn+j)
      a(jn+j) = a(jr+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jc+j) - 0.25 * t5
      a(jc+j) = a(jc+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jo+j) = a(jw+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      u1 = b(jh+j) + b(jo+j)
      u2 = b(jm+j) + b(jn+j)
      u3 = b(jh+j) - b(jo+j)
      u4 = b(jm+j) - b(jn+j)
      b(jn+j) = b(jr+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jc+j) - 0.25 * u5
      b(jc+j) = b(jc+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jo+j) = b(jw+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      a(jh+j) = co1*(t8-u11) - si1*(u8+t11)
      b(jh+j) = si1*(t8-u11) + co1*(u8+t11)
      a(jw+j) = co4*(t8+u11) - si4*(u8-t11)
      b(jw+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jm+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jm+j) = si2*(t9-u10) + co2*(u9+t10)
      a(jr+j) = co3*(t9+u10) - si3*(u9-t10)
      b(jr+j) = si3*(t9+u10) + co3*(u9-t10)
*----------------------
      t1 = a(ji+j) + a(jt+j)
      t2 = a(jn+j) + a(js+j)
      t3 = a(ji+j) - a(jt+j)
      t4 = a(jn+j) - a(js+j)
      a(jt+j) = a(jx+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(jp+j) - 0.25 * t5
      ax = a(jp+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(jp+j) = a(jd+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      a(jd+j) = ax
      u1 = b(ji+j) + b(jt+j)
      u2 = b(jn+j) + b(js+j)
      u3 = b(ji+j) - b(jt+j)
      u4 = b(jn+j) - b(js+j)
      b(jt+j) = b(jx+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(jp+j) - 0.25 * u5
      bx = b(jp+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(jp+j) = b(jd+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      b(jd+j) = bx
      a(ji+j) = co1*(t8-u11) - si1*(u8+t11)
      b(ji+j) = si1*(t8-u11) + co1*(u8+t11)
      a(jx+j) = co4*(t8+u11) - si4*(u8-t11)
      b(jx+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jn+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jn+j) = si2*(t9-u10) + co2*(u9+t10)
      a(js+j) = co3*(t9+u10) - si3*(u9-t10)
      b(js+j) = si3*(t9+u10) + co3*(u9-t10)
*----------------------
      t1 = a(jv+j) + a(jy+j)
      t2 = a(jo+j) + a(jt+j)
      t3 = a(jv+j) - a(jy+j)
      t4 = a(jo+j) - a(jt+j)
      a(jv+j) = a(jj+j)
      t5 = t1 + t2
      t6 = c1 * ( t1 - t2 )
      t7 = a(ju+j) - 0.25 * t5
      ax = a(ju+j) + t5
      t8 = t7 + t6
      t9 = t7 - t6
      a(ju+j) = a(je+j)
      t10 = c3 * t3 - c2 * t4
      t11 = c2 * t3 + c3 * t4
      a(je+j) = ax
      u1 = b(jv+j) + b(jy+j)
      u2 = b(jo+j) + b(jt+j)
      u3 = b(jv+j) - b(jy+j)
      u4 = b(jo+j) - b(jt+j)
      b(jv+j) = b(jj+j)
      u5 = u1 + u2
      u6 = c1 * ( u1 - u2 )
      u7 = b(ju+j) - 0.25 * u5
      bx = b(ju+j) + u5
      u8 = u7 + u6
      u9 = u7 - u6
      b(ju+j) = b(je+j)
      u10 = c3 * u3 - c2 * u4
      u11 = c2 * u3 + c3 * u4
      b(je+j) = bx
      a(jj+j) = co1*(t8-u11) - si1*(u8+t11)
      b(jj+j) = si1*(t8-u11) + co1*(u8+t11)
      a(jy+j) = co4*(t8+u11) - si4*(u8-t11)
      b(jy+j) = si4*(t8+u11) + co4*(u8-t11)
      a(jo+j) = co2*(t9-u10) - si2*(u9+t10)
      b(jo+j) = si2*(t9-u10) + co2*(u9+t10)
      a(jt+j) = co3*(t9+u10) - si3*(u9-t10)
      b(jt+j) = si3*(t9+u10) + co3*(u9-t10)
      j = j + jump
  440 continue
*
      endif
*
*-----(end of loop across transforms)
*
      ja = ja + jstepx
      if (ja.lt.istart) ja = ja + ninc
  445 continue
  450 continue
  460 continue
*-----( end of double loop for this k )
      kk = kk + 2*la
  470 continue
*-----( end of loop over values of k )
      la = 5*la
  480 continue
*-----( end of loop on type ii radix-5 passes )
*-----( nvex transforms completed)
  490 continue
      istart = istart + nvex * jump
  500 continue
*-----( end of loop on blocks of transforms )
*
      return
      end
      subroutine primfac(n,ifax)
      integer n,m,f2,f3,f5,ifax(*)
      m=n
      f2=0
      f3=0
      f5=0
   2  if(mod(m,2).eq.0)then
         f2=f2+1
         m=m/2
         go to 2
      endif
   3  if(mod(m,3).eq.0)then
         f3=f3+1
         m=m/3
         go to 3
      endif
   5  if(mod(m,5).eq.0)then
         f5=f5+1
         m=m/5
         go to 5
      endif
      if(m.ne.1)then
         print *,'primfac : bad prime factor in n'
         stop
      endif
      ifax(2)=f2
      ifax(3)=f3
      ifax(5)=f5
      return
      end
************************************************************************
      subroutine set99(ex,ifax,n)
      call fftrig(ex,n,3)
      call fax(ifax,n,3)
      return
      end
*     subroutine 'trig37' - sets up trigs for gpfa routines
*
*     n = (2**ip) * (3**iq) * (5**ir)
*
*----------------------------------------------------------------------
*
      subroutine trig37(trigs,n,ip,iq,ir)
      dimension trigs(*)
      dimension nj(3)
*
      nj(1) = 2**ip
      nj(2) = 3**iq
      nj(3) = 5**ir
      twopi = 4.0 * asin(1.0)
      ll = 1
*
      do 20 i = 1 , 3
      ni = nj(i)
      if (ni.eq.1) go to 20
*
      del = twopi / float(ni)
      irot = n / ni
      kink = mod(irot,ni)
      kk = 0
*
      do 10 k = 1 , ni
      angle = float(kk) * del
      trigs(ll) = cos(angle)
      trigs(ll+1) = sin(angle)
      ll = ll + 2
      kk = kk + kink
      if (kk.gt.ni) kk = kk - ni
   10 continue
   20 continue
*
      return
      end
************************************************************************
      subroutine vpassm(a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la)
c
c     routine 'vpassm' - multiple version of 'vpassa'
c     performs one pass through data
c     as part of multiple complex fft routine
c     a is first real input vector
c     b is first imaginary input vector
c     c is first real output vector
c     d is first imaginary output vector
c     trigs is precalculated table of sines ' cosines
c     inc1 is addressing increment for a and b
c     inc2 is addressing increment for c and d
c     inc3 is addressing increment between a's & b's
c     inc4 is addressing increment between c's & d's
c     lot is the number of vectors
c     n is length of vectors
c     ifac is current factor of n
c     la is product of previous factors
c
c     routine vpassm(a,b,c,d,trigs,inc1,inc2,inc3,inc4,lot,n,ifac,la)
c
c     end
      dimension a(n),b(n),c(n),d(n),trigs(n)
      data sin36/0.587785252292473/,cos36/0.809016994374947/,
     *     sin72/0.951056516295154/,cos72/0.309016994374947/,
     *     sin60/0.866025403784437/
c
      m=n/ifac
      iink=m*inc1
      jink=la*inc2
      jump=(ifac-1)*jink
      ibase=0
      jbase=0
      igo=ifac-1
c     check factors are correct - ensure non-negative
      if (igo.le.0) goto 998
      if (igo.gt.4) go to 999
      go to (10,50,90,130),igo
c
c     coding for factor 2
c
   10 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      do 20 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 15 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=a(ia+i)-a(ib+i)
      d(jb+j)=b(ia+i)-b(ib+i)
      i=i+inc3
      j=j+inc4
   15 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   20 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 40 k=la1,m,la
      kb=k+k-2
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      do 30 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 25 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
      d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
      i=i+inc3
      j=j+inc4
   25 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   30 continue
      jbase=jbase+jump
   40 continue
      return
c
c     coding for factor 3
c
   50 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      do 60 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 55 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
      c(jc+j)=(a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
      d(jb+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
      d(jc+j)=(b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
      i=i+inc3
      j=j+inc4
   55 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   60 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 80 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      do 70 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 65 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=
     *    c1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i))))
     *   -s1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      d(jb+j)=
     *    s1*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i))))
     *   +c1*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      c(jc+j)=
     *    c2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i))))
     *   -s2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      d(jc+j)=
     *    s2*((a(ia+i)-0.5*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i))))
     *   +c2*((b(ia+i)-0.5*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      i=i+inc3
      j=j+inc4
   65 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   70 continue
      jbase=jbase+jump
   80 continue
      return
c
c     coding for factor 4
c
   90 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      do 100 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 95 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      d(jc+j)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
      c(jb+j)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
      c(jd+j)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
      d(jb+j)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
      d(jd+j)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
      i=i+inc3
      j=j+inc4
   95 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  100 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 120 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      do 110 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 105 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      c(jc+j)=
     *    c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      d(jc+j)=
     *    s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      c(jb+j)=
     *    c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))
     *   -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      d(jb+j)=
     *    s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))
     *   +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      c(jd+j)=
     *    c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))
     *   -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      d(jd+j)=
     *    s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))
     *   +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  105 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  110 continue
      jbase=jbase+jump
  120 continue
      return
c
c     coding for factor 5
c
  130 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      ie=id+iink
      je=jd+jink
      do 140 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 135 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *  -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *  +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *  +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *  -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *  -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *  +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *  +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *  -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  135 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  140 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 160 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      ke=kd+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      c4=trigs(ke+1)
      s4=trigs(ke+2)
      do 150 l=1,la
      i=ibase
      j=jbase
cdir$ ivdep
      do 145 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=
     *    c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(jb+j)=
     *    s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(je+j)=
     *    c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(je+j)=
     *    s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(jc+j)=
     *    c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jc+j)=
     *    s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      c(jd+j)=
     *    c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jd+j)=
     *    s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      i=i+inc3
      j=j+inc4
  145 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  150 continue
      jbase=jbase+jump
  160 continue
      return
c  ** error - factor less than 1  not allowed **
998   print *,' fft99: factors are incorrect '
      call abort
c  ** error - factor higher than 5 not allowed **
999   print *,' fft99: factors higher than 5 are not supported '
      call abort
      end
**************************************************************************
      subroutine fftcc(a,w,ex,ifax,inc,jump,n,nft,isign)
      real a(*)
      call cfftmlt(a(1),a(2),w,ex,ifax,inc,jump,n,nft,isign)
      return
      end  

