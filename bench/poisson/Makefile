FFLAGS  = `pkg-config --cflags neko`
LIBS    = `pkg-config --libs neko`
FC      = `pkg-config --variable=compiler neko`

DEST    = poisson
SRC	= setup.f90 driver.f90
OBJ	= ${SRC:.f90=.o}

all: $(DEST)

install:

clean:
	-rm -f *.o core *.core $(OBJ) $(DEST)

$(DEST): ${OBJ}
	$(FC) $(FFLAGS) ${OBJ} -o $@  $(LIBS)

%.o: %.f90
	${FC} ${FFLAGS} -c $<


