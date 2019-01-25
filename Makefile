# B145208 #

MF=	Makefile

#morar
CC=	mpicc
CFLAGS = -O3

LFLAGS=	-lm

EXE	=	image

SRC= \
	par_image.c \
  src/arralloc.c \
	src/pgmio.c \
  src/boundary.c \
  src/average.c \

INC=\
  src/arralloc.h \
	src/pgmio.h \
  src/boundary.h \
  src/average.h \
  
OUT=output

#
# No need to edit below this line
#


.SUFFIXES:
.SUFFIXES: .c .o

OBJ=	$(SRC:.c=.o)

.c.o:
	$(CC) $(CFLAGS) -o $(<:.c=.o) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(INC)

$(OBJ):	$(MF)

clean:
	rm -rf $(OBJ) $(EXE) $(OUT)/* core
