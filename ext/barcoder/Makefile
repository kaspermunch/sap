LIBS     =   -lm

CFLAGS   =  -O4 

CC       =   g++

OBJECTS  =   main.o alignment.o condlikes.o constraints.o iomanager.o MbBitfield.o MbEigensystem.o MbMath.o MbRandom.o MbTransitionMatrix.o mcmc.o model.o parm.o parm_basefreqs.o parm_shape.o parm_subrates.o parm_tree.o settings.o tiprobs.o

PROGS    = bc1

all:		$(PROGS)

bc1:		$(OBJECTS)
		$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o bc1
		
main.o:	main.cpp
		$(CC) $(CFLAGS) -c main.cpp

alignment.o:	alignment.cpp
		$(CC) $(CFLAGS) -c alignment.cpp

condlikes.o:	condlikes.cpp
		$(CC) $(CFLAGS) -c condlikes.cpp

constraints.o:	constraints.cpp
		$(CC) $(CFLAGS) -c constraints.cpp

iomanager.o:	iomanager.cpp
		$(CC) $(CFLAGS) -c iomanager.cpp

MbBitfield.o:	MbBitfield.cpp
		$(CC) $(CFLAGS) -c MbBitfield.cpp

MbEigensystem.o:	MbEigensystem.cpp
		$(CC) $(CFLAGS) -c MbEigensystem.cpp

MbMath.o:	MbMath.cpp
		$(CC) $(CFLAGS) -c MbMath.cpp

MbRandom.o:	MbRandom.cpp
		$(CC) $(CFLAGS) -c MbRandom.cpp

MbTransitionMatrix.o:	MbTransitionMatrix.cpp
		$(CC) $(CFLAGS) -c MbTransitionMatrix.cpp

mcmc.o:	mcmc.cpp
		$(CC) $(CFLAGS) -c mcmc.cpp

model.o:	model.cpp
		$(CC) $(CFLAGS) -c model.cpp

parm.o:	parm.cpp
		$(CC) $(CFLAGS) -c parm.cpp

parm_basefreqs.o:	parm_basefreqs.cpp
		$(CC) $(CFLAGS) -c parm_basefreqs.cpp

parm_shape.o:	parm_shape.cpp
		$(CC) $(CFLAGS) -c parm_shape.cpp

parm_subrates.o:	parm_subrates.cpp
		$(CC) $(CFLAGS) -c parm_subrates.cpp

parm_tree.o:	parm_tree.cpp
		$(CC) $(CFLAGS) -c parm_tree.cpp

settings.o:	settings.cpp
		$(CC) $(CFLAGS) -c settings.cpp

tiprobs.o:	tiprobs.cpp
		$(CC) $(CFLAGS) -c tiprobs.cpp


clean:		
		rm -f *.o
