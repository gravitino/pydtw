CC      = g++
CFLAGS  = -Wall -g -O3 -fPIC -march=native -std=c++11 -pthread -mavx
PYPATH  = /usr/include/python2.7

GLOBALMEASURESOBJ = global_measures.o global_measures_wrap.o

all: _libglobalmeasures.so

global_measures_wrap.cxx: global_measures.i numpy.i global_measures.cpp global_measures.hpp helper.hpp
	swig -python -c++ global_measures.i

global_measures.o: global_measures.cpp global_measures.hpp helper.hpp
	$(CC) $(CFLAGS) -c global_measures.cpp

global_measures_wrap.o: global_measures_wrap.cxx
	$(CC) $(CFLAGS) -c global_measures_wrap.cxx -I $(PYPATH) 

_libglobalmeasures.so: global_measures.o global_measures_wrap.o
	$(CC) $(CFLAGS) $(GLOBALMEASURESOBJ) -shared -o _libglobalmeasures.so
	rm -f global_measures_wrap.cxx
	rm -f $(GLOBALMEASURESOBJ)

clean:
	rm -f libglobalmeasures.py
	rm -f libglobalmeasures.pyc
	rm -f _libglobalmeasures.so


