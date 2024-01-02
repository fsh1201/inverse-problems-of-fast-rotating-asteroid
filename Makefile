MAKE=make
CC=gcc
LD=ld
AR=ar

OPTFLAGS=-O3 -g
CFLAGS=$(OPTFLAGS) -Wall
LDFLAGS=-lm
MY_LIBS=trifac.o areanorm.o sphfunc.o ellfit.o covsrt.o ludcmp.o lubksb.o mrqmin.o mrqcof.o\
        curv.o blmatrix.o conv.o gauss.o phasec.o matrix.o bright.o memory.o\
	dot_product.o 

all: convexinv period_scan spin_axis minkowski standardtri shape2obj

libs: $(MY_LIBS)

convexinv: convexinv.c $(MY_LIBS)
	$(CC) $(CFLAGS) -o $@ $< $(MY_LIBS) $(LDFLAGS) -g

period_scan: period_scan.c $(MY_LIBS)
	$(CC) $(CFLAGS) -o $@ $< $(MY_LIBS) $(LDFLAGS) -g

spin_axis: spin_axis.cpp
	g++ spin_axis.cpp -o spin_axis

%.o: %.c
	$(CC) $(CFLAGS) -c $<

minkowski: minkowski.f
	gfortran ./minkowski.f  -o minkowski

standardtri: standardtri.f
	gfortran ./standardtri.f -o standardtri

shape2obj: shape2obj.cpp
	g++ shape2obj.cpp -o shape2obj

clean:
	rm -f *.o period_scan convexinv spin_axis minkowski standardtri shape2obj
