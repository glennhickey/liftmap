# liftmap -- see doc/SPEC.md
CC      ?= cc
CFLAGS  ?= -std=c99 -O2 -g -Wall -Wextra
CPPFLAGS += -D_POSIX_C_SOURCE=200809L -D_GNU_SOURCE
LDLIBS  += -lz -lpthread

SRC  = src/lmap_codec.c src/lmap_io.c src/lmap_io_curl.c src/lmap_file.c src/lmap_build.c src/lmap_chain.c

# make HTTP=1: an HTTP(S) backend over libcurl, so lmap_open (and bin/liftmap) read URLs.
ifeq ($(HTTP),1)
CPPFLAGS += -DLMAP_HTTP
LDLIBS   += -lcurl
endif
OBJ  = $(SRC:.c=.o)
LIB  = libliftmap.a
TESTBINS = bin/xcheck bin/roundtrip bin/robust bin/query bin/fuzz_reader bin/sections bin/remote

all: $(LIB) bin/liftmap $(TESTBINS)

$(LIB): $(OBJ)
	$(AR) rcs $@ $(OBJ)

bin/%: tests/%.c $(LIB) | bin
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(LIB) $(LDLIBS)

bin/liftmap: tools/liftmap.c $(LIB) | bin
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $< $(LIB) $(LDLIBS)

bin:
	mkdir -p bin

# Self-contained: generated data only.  roundtrip/query/xcheck take real runs from a caller
# (see their usage lines); robust also runs inside test_tools on a generated file.
check: all
	bin/sections $${TMPDIR:-/tmp}/liftmap_sections_test.lmap
	python3 tests/test_tools.py bin

clean:
	rm -f $(OBJ) $(LIB); rm -rf bin

.PHONY: all check clean
