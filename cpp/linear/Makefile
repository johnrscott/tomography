CC = g++
LDLIBS = -I../ -I../qpp -L/usr/local/lib/ -lnlopt
CFLAGS = -std=c++11 -O3 -c
LFLAGS = -std=c++11 -O3
BINARY = linear

linear: linear.o estimation.o stats.o progress.o
	$(CC) $(LFLAGS) -o $(BINARY) linear.o estimation.o stats.o progress.o $(LDLIBS)

linear.o: linear.cpp linear.h 
	$(CC) $(CFLAGS) linear.cpp $(LDLIBS)

estimation.o: estimation.cpp estimation.h
	$(CC) $(CFLAGS) estimation.cpp $(LDLIBS)

stats.o: stats.cpp stats.h
	$(CC) $(CFLAGS) stats.cpp $(LDLIBS)

progress.o: progress.cpp progress.h
	$(CC) $(CFLAGS) progress.cpp $(LDLIBS)

clean:
	rm -f $(BINARY) *.o 
