# Mac OS X
CFLAGS = -Wall 
veclib_incl = -I/System/Library/Frameworks/vecLib.Framework/Headers/
veclib = -Wl,-framework -Wl,veclib  

all: test_linalg.e h2.e

clean:
	$(RM) *~ *.o *.e

test: test_linalg.e h2.e
	./test_linalg.e
	./h2.e

test_linalg.e: test_linalg.o linalg.o
	$(CC) $(veclib_incl) -o $@ $^ $(veclib)

h2o.e: h2o.o cints.o qc.o linalg.o
	$(CC) -o $@ $^ 

h2.e: h2.o cints.o qc.o linalg.o
	$(CC) -o $@ $^ 

h2_sto.e: h2_sto.o cints.o qc.o linalg.o
	$(CC) -o $@ $^ 
