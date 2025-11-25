LIBS = -lm
CC = cc -L$(HOME)/lib -I$(HOME)/include
EXE = flexcalc

$(EXE) : flexcalc.o
	$(CC) -o $@ $< $(LIBS)

.c.o :
	$(CC) -c -o $@ $<

clean :
	\rm -f *.o $(EXE)
