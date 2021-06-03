

BASEDIR = $(shell dirname $$(dirname $$(dirname $$(which nvcc))))
CUFFTDIR = $(BASEDIR)/math_libs/10.2/lib64
CUDALIBDIR = $(BASEDIR)/cuda/11.1/targets/x86_64-linux/lib
CULIBS = -lcufft -lcuda -lcudart
CUINC = $(BASEDIR)/math_libs/10.2/include

all : test.exe

getaccstrm.o : getaccstrm.c
	nvc -acc -c -I$(CUINC) -o $@ $<

cufftwrapper.o : cufftwrapper.c
	nvcc -c -I$(CUINC) -o $@ $<

getplan.o : getplan.c
	nvcc -c -I$(CUINC) -o $@ $<

fortranfft.o : fortranfft.f90
	nvfortran -c -r8 -acc -o $@ $^

test.exe : test.f90 fortranfft.o cufftwrapper.o getaccstrm.o getplan.o
	nvfortran -r8 -acc -o $@ $^ -L$(CUFFTDIR) -L$(CUDALIBDIR) $(CULIBS)

tst.exe: tst.c
	nvcc -o $@ $^ -I$(CUINC) -L$(CUFFTDIR) $(CULIBS)

clean:
	rm -f *.o *.exe
