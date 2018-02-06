.SUFFIXES: .o .f .F

.F.o:
	$(CPP) $(CPPFLAGS) -P $*.F > $*.f
	$(CFT) -c $(FFLAGS) $*.f

.f.o:
	$(CFT) -c $(FFLAGS) $<

.F.f:
	$(CPP) $(CPPFLAGS) -P $*.F > $*.f

.F:
	$(CFT) -o $@ $(FFLAGS) $<

.f:
	$(CFT) -o $@ $(FFLAGS) $<

          RM = /bin/rm
         CPP = /usr/bin/cpp
    CPPFLAGS = -traditional
         CFT = gfortran -fdefault-real-8  -O3  -march=native -ffast-math -funroll-loops -finline-functions -mssse3 -mcmodel=medium -flto -fwhole-program -fno-protect-parens -fopenmp  # -I/opt/netcdf/netcdf3/include 
        LIBS =  # -I/opt/netcdf/netcdf3/include  -L/opt/netcdf/netcdf3/lib  -lnetcdff -lnetcdf

RBIN = ns3d

SRCS =	dcopy.f				\
	defk.F				\
	divergence.F			\
	eulero.F			\
	fft.F				\
	force_velocity.F     		\
	ifft.F				\
	initial.F			\
	initfft.F			\
	leapfo.F			\
	navisto.F			\
	nonlidns.F			\
	precal.F			\
	project.F			\
	rando.F				\
	rdmfld.F			\
	rotv.F				\
	runge.F				\
	statistics.F			\
	trunck.F			\
	vecpro.F			\
	zerkam.F			\
	zero.F				\
	fftx.F				\
	ffty.F				\
	fftz.F				\
	energy.f			\
	temperton.real.F

#	write_stat.F			\
#	write_field.F			\

OBJS=$(SRCS:.F=.o)

$(RBIN): ${OBJS}
		${CFT} ${FFLAGS} -o $(RBIN) ${OBJS} ${LIBS}

clean:
		$(RM) -f *.o


