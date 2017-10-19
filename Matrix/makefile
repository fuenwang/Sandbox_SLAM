all:
	#g++ Matrix.cpp MAT.cpp VEC.cpp -lboost_numpy -lboost_python -lpython2.7
	g++ -c -fPIC Matrix.cpp MAT.cpp VEC.cpp
	g++ -shared -o Matrix.so Matrix.o MAT.o VEC.o  -lboost_numpy -lboost_python -lpython2.7
	rm *.o
