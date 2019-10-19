
FC=gfortran
FCFLAGS=-std=legacy -g -O0 -fno-bounds-check

SRCS=D1.FOR D2.FOR D3.FOR D4.FOR D5.FOR D6.FOR
OBJS=$(SRCS:.FOR=.o)
PROG=AD

all: $(SRCS) $(PROG)

$(PROG): $(OBJS)
	$(FC) $(OBJS) -o $(PROG)

$(OBJS): $(SRCS)
	$(FC) $(FCFLAGS) -c $^

clean:
	rm -f $(PROG) *.o

