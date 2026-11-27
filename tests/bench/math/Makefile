FFLAGS  = `pkg-config --cflags neko`
LIBS    = `pkg-config --libs neko`
FC      = `pkg-config --variable=compiler neko`
BCKND   = `pkg-config --variable=backend neko`

FFLAGS  += -DNEKO_BCKND_$(BCKND)

DEST    = mathbench
SRC	= driver.F90
OBJ	= ${SRC:.F90=.o}

all: $(DEST)

install:

clean:
	-rm -f *.o core *.core *.mod $(OBJ) $(DEST)

$(DEST): ${OBJ}
	$(FC) $(FFLAGS) ${OBJ} -o $@  $(LIBS)

%.o: %.F90
	${FC} ${FFLAGS} -c $<
