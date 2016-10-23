all:
	g++ solver.cpp -llapacke -llapack -I/usr/include/lapacke
	#icpc -mkl -qopenmp solver.cpp
	./a.out
	#python one.py  

