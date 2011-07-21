CC=gcc
CCFLAGS=-Wall -dynamiclib -arch i386 -arch x86_64

libenvironhack.dylib: libenvironhack.c
	$(CC) $(CCFLAGS) -init _init_environ $< -o $@

all: libenvironhack.dylib

clean:
	rm -f libenvironhack.dylib
