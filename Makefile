# libintervalmap -- see doc/SPEC.md
CC      ?= cc
CFLAGS  ?= -std=c99 -O2 -g -Wall -Wextra
CPPFLAGS += -D_POSIX_C_SOURCE=200809L
LDLIBS  += -lz

SRC  = src/imap_codec.c src/imap_io.c src/imap_file.c
OBJ  = $(SRC:.c=.o)
LIB  = libintervalmap.a
TESTBINS = bin/xcheck bin/roundtrip bin/robust bin/query bin/fuzz_reader

all: $(LIB) $(TESTBINS)

$(LIB): $(OBJ)
	$(AR) rcs $@ $(OBJ)

bin/%: tests/%.c $(LIB) | bin
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(LIB) $(LDLIBS)

bin:
	mkdir -p bin

# Round-trip and robustness run on generated data; xcheck needs runs from a caller.
check: $(TESTBINS)
	@bin/robust $${IMAP_TEST_FILE:-/dev/null} 2>/dev/null || \
	  echo "check: set IMAP_TEST_FILE to an .imap file (see tests/README)"

clean:
	rm -f $(OBJ) $(LIB); rm -rf bin

.PHONY: all check clean
