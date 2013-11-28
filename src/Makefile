GCC=gcc -std=gnu99
C_FLAGS=-fpic -O2 -pipe -g -c
SOURCES= r_ehh.c r_ehhs.c ehh_utils.c r_scan_hh.c
LD_FLAGS=-I/usr/share/R/include -I/usr/include/gperftools
LIBS= -lR -lprofiler
OBJECTS=$(SOURCES:.c=.o)
SHARED_LIBRARY=rehh.so

all: $(SOURCES) $(SHARED_LIBRARY)

$(SHARED_LIBRARY): $(OBJECTS) 
	$(GCC) $(LD_FLAGS) -shared $(OBJECTS) -o $(SHARED_LIBRARY) $(LIBS)

.c.o:
	$(GCC) $(C_FLAGS) $(LD_FLAGS) $<
clean:
	rm -f $(OBJECTS) $(SHARED_LIBRARY)
