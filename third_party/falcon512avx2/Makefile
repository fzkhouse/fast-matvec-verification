# This Makefile compiles the implementation in this directory along with
# the known answer tests generator located in the
# ../../../KAT/generator/ directory.

.POSIX:

CC      = c99
LD      = c99
CFLAGS  = -W -Wall -O2 -mavx2
LDFLAGS =
LIBS    =

ALGNAME = falcon512avx2

BUILD   = build

HEAD1 = api.h inner.h fpr.h
HEAD2 = api.h

# Core Falcon objects (no nist.c, so no randombytes dependency)
OBJ_CORE = \
	$(BUILD)/codec.o \
	$(BUILD)/common.o \
	$(BUILD)/fft.o \
	$(BUILD)/fpr.o \
	$(BUILD)/keygen.o \
	$(BUILD)/rng.o \
	$(BUILD)/shake.o \
	$(BUILD)/sign.o \
	$(BUILD)/vrfy.o

# KAT extra objects
OBJ_NIST = $(BUILD)/nist.o
OBJ_KAT  = $(BUILD)/PQCgenKAT_sign.o $(BUILD)/katrng.o

# Bench object
OBJ_BENCH_NATIVE = $(BUILD)/bench_matvec_falcon_native.o

# KAT link set (includes nist)
OBJ1 = $(OBJ_CORE) $(OBJ_NIST)
OBJ2 = $(OBJ_KAT)

all: $(BUILD) $(BUILD)/kat512avx2 $(BUILD)/bench_matvec_falcon_native

$(BUILD):
	mkdir -p $(BUILD)

# ===== KAT target =====
$(BUILD)/kat512avx2: $(OBJ1) $(OBJ2)
	$(LD) $(LDFLAGS) -o $@ $(OBJ1) $(OBJ2) $(LIBS)

# ===== Bench target =====
# IMPORTANT: do not link nist.o here (avoids undefined randombytes)
$(BUILD)/bench_matvec_falcon_native: $(OBJ_CORE) $(OBJ_BENCH_NATIVE)
	$(LD) $(LDFLAGS) -o $@ $(OBJ_CORE) $(OBJ_BENCH_NATIVE) $(LIBS)

# ===== Core compilation rules =====
$(BUILD)/codec.o: codec.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o $@ codec.c

$(BUILD)/common.o: common.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o $@ common.c

$(BUILD)/fft.o: fft.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o $@ fft.c

$(BUILD)/fpr.o: fpr.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o $@ fpr.c

$(BUILD)/keygen.o: keygen.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o $@ keygen.c

$(BUILD)/nist.o: nist.c $(HEAD2)
	$(CC) $(CFLAGS) -c -o $@ nist.c

$(BUILD)/rng.o: rng.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o $@ rng.c

$(BUILD)/shake.o: shake.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o $@ shake.c

$(BUILD)/sign.o: sign.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o $@ sign.c

$(BUILD)/vrfy.o: vrfy.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o $@ vrfy.c

# ===== KAT generator rules =====
$(BUILD)/PQCgenKAT_sign.o: ../../../KAT/generator/PQCgenKAT_sign.c $(HEAD2)
	$(CC) $(CFLAGS) -I . -DALGNAME=$(ALGNAME) -c -o $@ ../../../KAT/generator/PQCgenKAT_sign.c

$(BUILD)/katrng.o: ../../../KAT/generator/katrng.c
	$(CC) $(CFLAGS) -I . -c -o $@ ../../../KAT/generator/katrng.c

# ===== Bench rules =====
$(BUILD)/bench_matvec_falcon_native.o: bench_matvec_falcon_native.c $(HEAD1)
	$(CC) $(CFLAGS) -c -o $@ bench_matvec_falcon_native.c

clean:
	-rm -f \
		$(BUILD)/kat512avx2 \
    	$(BUILD)/bench_matvec_falcon_native \
    	$(OBJ_CORE) $(OBJ_NIST) $(OBJ_KAT) $(OBJ_BENCH_NATIVE)

.PHONY: all clean