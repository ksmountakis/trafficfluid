CFLAGS=-O3 -Wall -fopenmp 
LDFLAGS=-lm

all: sim


sim:sim.c sim_forces.c Makefile
	$(CC) $(CFLAGS) sim.c -o sim $(LDFLAGS)
