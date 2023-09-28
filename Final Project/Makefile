CC = gcc

CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors -std=c99

LIBS = -lm

symnmf: symnmf.o symnmf.h
	$(CC) -o symnmf $(CFLAGS) symnmf.o $(LIBS)

symnmf.o: symnmf.c
	$(CC) -c $(CFLAGS) symnmf.c $(LIBS)