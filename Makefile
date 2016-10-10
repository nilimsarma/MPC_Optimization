IDIR = ./
CC = gcc
CFLAGS = -I$(IDIR)

ODIR = ./
LDIR = 

LIBS = -lm

_DEPS = ip.hpp ip_iter.hpp ip_primal_dual_dir.hpp mpc_discretize.hpp mpc_formulate.hpp
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = ip_iter.cpp ip_primal_dual_dir.cpp mpc_discretize.cpp mpc_formulate.cpp
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

hellomake: $(OBJ)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~ 