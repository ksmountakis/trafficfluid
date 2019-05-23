CFLAGS=-O3 -Wall -fopenmp 
LDFLAGS=-lm
UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
	CC=gcc-8
endif

all: sim


sim:sim.c sim_forces.c Makefile
	$(CC) $(CFLAGS) sim.c -o sim $(LDFLAGS)
