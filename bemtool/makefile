GCC  = g++ -std=c++14 -O3 
INCLUDE = -I/usr/include/eigen3/ -I./bemtool/ -I../

GCCFLAGS = -L/usr/lib64/ -DQR=${QUAD} -DQR2=${QUAD2}
LIBS = -lcblas -llapack

all: testEM

######################################################


testFullD: testFullD.o
	$(GCC) $(GCCFLAGS) -fopenmp testFullD.o -o testFullD $(LIBS)

testFullD.o: testFullD.cpp
	$(GCC) $(GCCFLAGS) $(INCLUDE) -fopenmp -c testFullD.cpp -o testFullD.o $(LIBS)


testFEM: testFEM.o
	$(GCC) $(GCCFLAGS) -fopenmp testFEM.o -o testFEM $(LIBS)

testFEM.o: testFEM.cpp
	$(GCC) $(GCCFLAGS) $(INCLUDE) -fopenmp -c testFEM.cpp -o testFEM.o $(LIBS)

testEM: testEM.o
	$(GCC) $(GCCFLAGS) -fopenmp testEM.o -o testEM $(LIBS)

testEM.o: testEM.cpp
	$(GCC) $(GCCFLAGS) $(INCLUDE) -fopenmp -c testEM.cpp -o testEM.o $(LIBS)

testBEM: testBEM.o
	$(GCC) $(GCCFLAGS) -fopenmp $^ -o $@ $(LIBS)

testBEM.o: testBEM.cpp
	$(GCC) $(GCCFLAGS) $(INCLUDE) -fopenmp -c $^ -o $@ $(LIBS)

clean:
	rm *.o
	rm output/*.txt
