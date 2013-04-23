CC      = g++
CFLAGS  = -Wall -g -O3 -fPIC -march=native
OBJ     = dtw.o dtw_wrap.o
PYPATH  = /usr/include/python2.7

all: _libdtw.so

dtw_wrap.cxx: dtw.i dtw.cpp dtw.hpp
	swig -python -c++ dtw.i

dtw.o: dtw.cpp dtw.hpp
	$(CC) $(CFLAGS) -c dtw.cpp

dtw_wrap.o: dtw_wrap.cxx
	$(CC) $(CFLAGS) -c dtw_wrap.cxx -I $(PYPATH) 

_libdtw.so: dtw.o dtw_wrap.o
	$(CC) $(CFLAGS) $(OBJ) -shared -o _libdtw.so
	rm -f dtw_wrap.cxx
	rm -f $(OBJ)
clean:
	rm -f libdtw.py
	rm -f libdtw.pyc
	rm -f _libdtw.so


