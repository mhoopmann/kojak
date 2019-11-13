#Set these variables if needed
C = gcc
CC = g++
FLAGS = -O3 -static -std=c++11 -D_NOSQLITE -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DGCC

#Paths to supporting software
MSTOOLKITPATH = ../MSToolkit
HARDKLORPATH = ../Hardklor

#Do not touch these variables
LIBPATH = -L$(MSTOOLKITPATH) -L$(HARDKLORPATH)
LIBS = -lmstoolkitlite -lhardklor -lpthread
INCLUDE = -I$(MSTOOLKITPATH)/include -I$(HARDKLORPATH)


#Do not touch these variables
KOJAK = KojakManager.o KParams.o KAnalysis.o KData.o KDB.o KPrecursor.o KSpectrum.o KIons.o KIonSet.o KLog.o KTopPeps.o Threading.o CometDecoys.o


#Make statements
kojak : Kojak.cpp $(KOJAK)
	$(CC) $(FLAGS) $(INCLUDE) $(KOJAK) Kojak.cpp $(LIBPATH) $(LIBS) -o kojak

clean:
	rm *.o kojak


#Kojak objects
KojakManager.o : KojakManager.cpp
	$(CC) $(FLAGS) $(INCLUDE) KojakManager.cpp -c
	
KParams.o : KParams.cpp
	$(CC) $(FLAGS) $(INCLUDE) KParams.cpp -c

KAnalysis.o : KAnalysis.cpp
	$(CC) $(FLAGS) $(INCLUDE) KAnalysis.cpp -c

KData.o : KData.cpp
	$(CC) $(FLAGS) $(INCLUDE) KData.cpp -c

KDB.o : KDB.cpp
	$(CC) $(FLAGS) $(INCLUDE) KDB.cpp -c

KPrecursor.o : KPrecursor.cpp
	$(CC) $(FLAGS) $(INCLUDE) KPrecursor.cpp -c

KSpectrum.o : KSpectrum.cpp
	$(CC) $(FLAGS) $(INCLUDE) KSpectrum.cpp -c

KIons.o : KIons.cpp
	$(CC) $(FLAGS) $(INCLUDE) KIons.cpp -c
		
KIonSet.o : KIonSet.cpp
	$(CC) $(FLAGS) $(INCLUDE) KIonSet.cpp -c
	
KLog.o : KLog.cpp
	$(CC) $(FLAGS) $(INCLUDE) KLog.cpp -c
	
KTopPeps.o : KTopPeps.cpp
	$(CC) $(FLAGS) $(INCLUDE) KTopPeps.cpp -c
		
Threading.o : Threading.cpp
	$(CC) $(FLAGS) $(INCLUDE) Threading.cpp -c
	
CometDecoys.o : CometDecoys.cpp
	$(CC) $(FLAGS) $(INCLUDE) CometDecoys.cpp -c

