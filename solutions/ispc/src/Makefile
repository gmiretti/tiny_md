OLEVEL  = 3
COMMONFLAGS = -g
CC      = gcc
CFLAGS	= -O$(OLEVEL)
WFLAGS	= -std=c11 -Wall -Wextra
LDFLAGS = -lm
ISPC    = ispc
ISPCFLAGS  = -O$(OLEVEL) --pic

TARGETS	= tiny_md viz
SOURCES	= $(shell echo *.c)
SOURCES_ISPC = $(shell echo *.ispc)
HEADERS_ISPC = $(SOURCES_ISPC:%.ispc=%_ispc.h)
OBJECTS = core.o wtime.o core_ispc.o

all: $(TARGETS)

viz: viz.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lGL -lGLU -lglut

tiny_md: tiny_md.o $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

%_ispc.o: %.ispc
	$(ISPC) $(COMMONFLAGS) $(ISPCFLAGS) $(PARAMFLAGS) -o $@ $<

%.o: %.c
	$(CC) $(WFLAGS) $(PARAMFLAGS) $(CFLAGS) -c $<

%_ispc.h: %.ispc
	$(ISPC) $(COMMONFLAGS) $(ISPCFLAGS) $(PARAMFLAGS) -h $(basename $<)_ispc.h $<

.depend_cc: $(SOURCES)
	$(CC) -MM $^ > $@

-include .depend_cc

.depend_ispc: $(SOURCES_ISPC) $(HEADERS_ISPC)
	$(ISPC) -M $(SOURCES_ISPC) > $@

-include .depend_ispc

clean:
	rm -f $(TARGETS) *.o .depend* *_ispc.h thermo.log trajectory.xyz

.PHONY: clean all
