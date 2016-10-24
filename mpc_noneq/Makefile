CC      = g++
CFLAGS  = -g -I /usr/include/lapacke
LDFLAGS = -L /usr/lib64 -lm -llapack -llapacke

my_mpc: ip_iter.o ip_primal_dual_dir.o mpc_discretize.o mpc_formulate.o utils.o
	$(CC) -o $@ $^ $(LDFLAGS)
	
ip_iter.o: ip_iter.cpp ip_iter.hpp ip.hpp
	$(CC) -c $(CFLAGS) $<

ip_primal_dual_dir.o: ip_primal_dual_dir.cpp ip_primal_dual_dir.hpp ip.hpp
	$(CC) -c $(CFLAGS) $<

mpc_discretize.o: mpc_discretize.cpp mpc_discretize.hpp ip.hpp
	$(CC) -c $(CFLAGS) $<
	
mpc_formulate.o: mpc_formulate.cpp mpc_formulate.hpp ip.hpp
	$(CC) -c $(CFLAGS) $<

utils.o: utils.cpp utils.hpp 
	$(CC) -c $(CFLAGS) $<
	
.PHONY: clean 

clean:
	rm -rf *.o my_mpc 