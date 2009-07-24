LIBS     =   -lm

CFLAGS   =  -O4 

CC       =   g++

OBJECTS  =   cConstrainedNJlib.o utils.o

PROGS    = cnj

all:		$(PROGS)

cnj:		$(OBJECTS)
		$(CC) $(CFLAGS) $(OBJECTS) $(LIBS) -o cnj
		
cConstrainedNJlib.o:	cConstrainedNJlib.cpp
		$(CC) $(CFLAGS) -c cConstrainedNJlib.cpp

utils.o:	utils.cpp
		$(CC) $(CFLAGS) -c utils.cpp

clean:		
		rm -f *.o

