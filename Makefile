EXE = csf2det.exe

CFLAGS = -O2 -Wall

LIB = -lm -largtable2

%.exe: %.o
	$(CC) -o $@ $< $(LIB)

default: $(EXE)

.PHONY: clean
clean:
	rm $(EXE)
