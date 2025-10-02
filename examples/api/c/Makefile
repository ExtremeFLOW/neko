CFLAGS  = `pkg-config --cflags jansson,neko`
LIBS    = `pkg-config --libs jansson,neko`
CC      = mpicc       # for compiling C code
FC      = mpifort     # for linking (Fortran libs need this)

DEST    = cylinder
SRC     = cylinder.c
OBJ     = ${SRC:.c=.o}

all: $(DEST)

clean:
	-rm -f *.o core *.core $(OBJ) $(DEST)

$(DEST): ${OBJ}
	$(FC) $(CFLAGS) ${OBJ} -o $@ $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<


