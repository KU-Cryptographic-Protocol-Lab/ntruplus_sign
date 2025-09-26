CC ?= /usr/bin/cc
#CFLAGS += -Wall -Wextra -Wpedantic -Wmissing-prototypes -Wredundant-decls \
#  -Wshadow -Wpointer-arith -O3 -fomit-frame-pointer
CFLAGS += -Wall -Wextra -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wpointer-arith -O3 -fomit-frame-pointer  
NISTFLAGS += -Wno-unused-result -O3 -fomit-frame-pointer
RM = /bin/rm

LDFLAGS = -lcrypto -lssl

SOURCES = poly.c ntt.c fft.c fips202.c packing.c sampler.c sign.c symmetric-shake.c randombytes.c decompose.c gaussian.c rejection.c encoding.c
HEADERS = params.h poly.h ntt.h fft.h sign.h api.h packing.h sampler.h fips202.h symmetric.h randombytes.h decompose.h gaussian.h rejection.h info.h encoding.h

.PHONY: all test clean

all: main

PQCgenKAT_kem: $(HEADERS) rng.h $(SOURCES) rng.c
	$(CC) $(NISTFLAGS) -o $@ PQCgenKAT_kem.c rng.c $(SOURCES) $(LDFLAGS)

main: $(HEADERS) $(SOURCES) main.c
	$(CC) $(CFLAGS) -o $@ main.c $(SOURCES) $(LDFLAGS)
	./main

test1: $(HEADERS) $(SOURCES) test1.c
	$(CC) $(CFLAGS) -o $@ test1.c $(SOURCES) $(LDFLAGS)
	
test2: $(HEADERS) $(SOURCES) test2.c
	$(CC) $(CFLAGS) -o $@ test2.c $(SOURCES) $(LDFLAGS)	

prec_encoding: $(HEADERS) $(SOURCES) prec_encoding.c
	$(CC) $(CFLAGS) -o $@ prec_encoding.c $(SOURCES) $(LDFLAGS)	

clean:
	-$(RM) main
