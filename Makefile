# --- External configuration ----------------------------------
CC       = g++
FC       = gfortran
CCFLAGS  = -O -Wall -fPIC -m64
CPPFLAGS = -I${ROOTSYS}/include

ROOTLIB  = $(shell root-config --libs)
CERNLIB  = -L/home/cupsoft/Products/CERNLIB/V2006/lib -lmathlib -lkernlib -lpacklib
# --- External configuration ----------------------------------

EXE = flat

default: exe

music-sr.o: music-sr.f
	${FC} -c music-sr.f

${EXE}.o: ${EXE}.cc
	${CC} ${CCFLAGS} ${CPPFLAGS} -c ${EXE}.cc -o ${EXE}.o

${EXE}: ${EXE}.o music-sr.o
	${CC} ${CCFLAGS} -o ${EXE} ${EXE}.o music-sr.o -lgfortran -lm ${CERNLIB} ${ROOTLIB}

exe: ${EXE}
	@rm -f *.o

clean:
	@rm -f *.o ${EXE}
