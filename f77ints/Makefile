FFLAGS = -Wall

all: fints.o test.e

test:
	./test.e

clean:
	rm -f *.o *.e *~

test.e: fints.o test.o
	$(FC) $(LFLAGS) $^ -o $@
