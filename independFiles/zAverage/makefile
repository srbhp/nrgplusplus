all :
	icpx -qmkl -qopenmp -xHost -g  -lmpfr -lgmp   -lhdf5 -lhdf5_cpp -fp-model=precise main.cpp hoppingElement.cpp
	./a.out

#icpc -mkl -qopenmp solver.cpp -xHost -g  -Wcheck -check-pointers=rw -traceback -std=c++0x  -lmpfr -lgmp   -lhdf5 -lhdf5_cpp -check-pointers-dangling=all -debug -fp-model consistent
#icpc -mkl -qopenmp solver.cpp  -xHost -rdynamic -g -DREPORT -Wcheck -check-pointers=rw -traceback -std=c++0x  -lmpfr -lgmp   -lhdf5 -lhdf5_cpp -fp-model consistent
#./a.out
#python3 plot.py

