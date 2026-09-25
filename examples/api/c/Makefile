CFLAGS  = `pkg-config --cflags jansson,neko`
LIBS    = `pkg-config --libs jansson,neko`
CC      = mpicc
FC      = mpif90

DEST    = cylinder
SRC	= cylinder.c
OBJ	= ${SRC:.c=.o}

all: $(DEST)

install:

clean:
	-rm -f *.o core *.core $(OBJ) $(DEST)

$(DEST): ${OBJ}
	${FC} $(CFLAGS) ${OBJ} -o $@  $(LIBS)
