# liftmap -- see doc/SPEC.md
CC      ?= cc
CFLAGS  ?= -std=c99 -O2 -g -Wall -Wextra
CPPFLAGS += -D_POSIX_C_SOURCE=200809L -D_GNU_SOURCE
LDLIBS  += -lz

SRC  = src/lmap_codec.c src/lmap_io.c src/lmap_file.c
OBJ  = $(SRC:.c=.o)
LIB  = libliftmap.a
TESTBINS = bin/xcheck bin/roundtrip bin/robust bin/query bin/fuzz_reader bin/sections

all: $(LIB) $(TESTBINS)

$(LIB): $(OBJ)
	$(AR) rcs $@ $(OBJ)

bin/%: tests/%.c $(LIB) | bin
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(LIB) $(LDLIBS)

bin:
	mkdir -p bin

# Round-trip and robustness run on generated data; xcheck needs runs from a caller.
check: $(TESTBINS)
	@bin/robust $${LMAP_TEST_FILE:-/dev/null} 2>/dev/null || \
	  echo "check: set LMAP_TEST_FILE to an .lmap file (see tests/README)"

clean:
	rm -f $(OBJ) $(LIB); rm -rf bin

.PHONY: all check clean
