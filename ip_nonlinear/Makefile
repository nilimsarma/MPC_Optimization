CC      = g++
CFLAGS  = -g -I /usr/include/lapacke
LDFLAGS = -L /usr/lib64 -lm -llapack -llapacke

my_mpc: ip_iter.o ip_primal_dual_dir.o ip_formulate.o utils.o
	$(CC) -o $@ $^ $(LDFLAGS)
	
ip_iter.o: ip_iter.cpp ip_iter.hpp ip.hpp
	$(CC) -c $(CFLAGS) $<

ip_primal_dual_dir.o: ip_primal_dual_dir.cpp ip_primal_dual_dir.hpp ip.hpp
	$(CC) -c $(CFLAGS) $<
	
ip_formulate.o: ip_formulate.cpp ip_formulate.hpp ip.hpp
	$(CC) -c $(CFLAGS) $<

utils.o: utils.cpp utils.hpp 
	$(CC) -c $(CFLAGS) $<
	
.PHONY: clean 

clean:
	rm -rf *.o my_mpc 
