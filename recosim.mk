# Makefile for the new ntuple processing
main: histset.o recosim.o ParTreeProcessing.C
	g++ -o compiledThreads ParTreeProcessing.C -pthread `root-config --cflags --libs`

histset.o: histset.C recosim.o Hungarian.h
	g++ -c -pthread histset.C `root-config --cflags --libs`

recosim.o: recosim.C recosim.h
	g++ -c -pthread recosim.C `root-config --cflags --libs`

hung.o: Hungarian.cpp Hungarian.h
	g++ -c -pthread Hungarian.cpp `root-config --cflags --libs`

clean:
	rm *.o
