FFLAGS  = `pkg-config --cflags neko`
LIBS    = `pkg-config --libs neko`
FC      = `pkg-config --variable=compiler neko`

DEST    = mathbench
SRC	= driver.f90
OBJ	= ${SRC:.f90=.o}

all: $(DEST)

install:

clean:
	-rm -f *.o core *.core *.mod $(OBJ) $(DEST)

$(DEST): ${OBJ}
	$(FC) $(FFLAGS) ${OBJ} -o $@  $(LIBS)

%.o: %.f90
	${FC} ${FFLAGS} -c $<
