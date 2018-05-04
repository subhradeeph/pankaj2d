#IDIR =../include
IDIR =.
CC=gcc -O3
CFLAGS=-I$(IDIR) -std=gnu99 -g

#ODIR=obj
ODIR=.
#LDIR =../lib
LDIR =.

LIBS=-lfftw3 -lm

_DEPS = binary.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o evolve.o get_input.o init_conf.o out_conf.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

spinodal: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 
