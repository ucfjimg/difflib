diff: diff.o main.o

diff.o: diff.c diff.h
main.o: main.c diff.h
