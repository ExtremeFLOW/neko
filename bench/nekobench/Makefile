FCLAGS  = `pkg-config --cflags neko`
LIBS    = `pkg-config --libs neko`
FC      = `pkg-config --variable=compiler neko`

DEST    = nekobench
SRC	= driver.f90
OBJ	= ${SRC:.f90=.o}

all: $(DEST)

install:

clean:
	-rm -f *.o core *.core $(OBJ) $(DEST)

$(DEST): ${OBJ}
	$(FC) $(FCLAGS) ${OBJ} -o $@  $(LIBS)

%.o: %.f90
	${FC} ${FCLAGS} -c $<


