
FC=gfortran
FCFLAGS=-g -O0 -fno-bounds-check

SRCS=PD1.FOR PD2.FOR PD3.FOR PD4.FOR PD5.FOR PD6.FOR PD7.FOR PD8.FOR PD9.FOR
OBJS=$(SRCS:.FOR=.o)
PROG=PD

all: $(SRCS) $(PROG)

$(PROG): $(OBJS)
	$(FC) $(OBJS) -o $(PROG)

$(OBJS): $(SRCS)
	$(FC) $(FCFLAGS) -c $^

clean:
	rm -f $(PROG) *.o

