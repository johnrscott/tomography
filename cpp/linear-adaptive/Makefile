CC = g++
LDLIBS = -I../ -L/usr/local/lib/ -lnlopt
CFLAGS = -std=c++11 -O3 -c
LFLAGS = -std=c++11 -O3
BINARY = linear-adaptive

linear-adaptive: linear-adaptive.o simulation.o estimation.o stats.o proj.o progress.o meas.o linear-adaptive.h
	$(CC) $(LFLAGS) -o $(BINARY) linear-adaptive.o simulation.o estimation.o stats.o proj.o progress.o meas.o $(LDLIBS)

linear-adaptive.o: linear-adaptive.cpp simulation.h 
	$(CC) $(CFLAGS) linear-adaptive.cpp $(LDLIBS)

simulation.o: simulation.cpp simulation.h
	$(CC) $(CFLAGS) simulation.cpp $(LDLIBS)

estimation.o: estimation.cpp estimation.h
	$(CC) $(CFLAGS) estimation.cpp $(LDLIBS)

stats.o: stats.cpp stats.h
	$(CC) $(CFLAGS) stats.cpp $(LDLIBS)

proj.o: proj.cpp proj.h
	$(CC) $(CFLAGS) -DDEBUG proj.cpp $(LDLIBS)

progress.o: progress.cpp progress.h
	$(CC) $(CFLAGS) progress.cpp $(LDLIBS)

meas.o: meas.cpp meas.h
	$(CC) $(CFLAGS) meas.cpp $(LDLIBS)

clean:
	rm -f $(BINARY) *.o 
