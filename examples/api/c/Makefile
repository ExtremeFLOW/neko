CFLAGS  = `pkg-config --cflags jansson,neko`
LIBS    = `pkg-config --libs jansson,neko`
CC      = mpicc

DEST    = demo
SRC	= driver.c
OBJ	= ${SRC:.c=.o}

all: $(DEST)

install:

clean:
	-rm -f *.o core *.core $(OBJ) $(DEST)

$(DEST): ${OBJ}
	$(CC) $(CFLAGS) ${OBJ} -o $@  $(LIBS)

%.o: %.f90
	${CC} ${CFLAGS} -c $<


