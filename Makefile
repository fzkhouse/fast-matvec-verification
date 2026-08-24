CC ?= cc
USE_AVX2 ?= 0
CFLAGS ?= -O3 -std=c11 -Wall -Wextra -Wshadow -Wconversion -pedantic
CPPFLAGS ?=
LDFLAGS ?=

ifeq ($(USE_AVX2),1)
CPPFLAGS += -DUSE_AVX2=1
CFLAGS += -mavx2
endif

.PHONY: all clean test
all: build/falcon-reference build/falcon-baseline build/kyber-reference build/kyber-baseline build/dilithium-reference build/dilithium-baseline build/dilithium5-gen build/dilithium2-reference build/dilithium2-baseline build/dilithium2-gen build/rpt-reference build/rpt-baseline

build:
	mkdir -p build

FALCON_DIR = third_party/falcon512avx2
FALCON_CORE = $(FALCON_DIR)/codec.c $(FALCON_DIR)/common.c \
	$(FALCON_DIR)/fft.c $(FALCON_DIR)/fpr.c $(FALCON_DIR)/keygen.c \
	$(FALCON_DIR)/rng.c $(FALCON_DIR)/shake.c $(FALCON_DIR)/sign.c \
	$(FALCON_DIR)/vrfy.c

# This is the complete original Falcon-native matvec baseline. It is kept
# separate from the portable verifier and is intentionally not simplified.
build/falcon-baseline: $(FALCON_CORE) $(FALCON_DIR)/bench_matvec_falcon_native.c | build
	$(CC) -O2 -mavx2 -std=c99 -I$(FALCON_DIR) -o $@ $(FALCON_CORE) $(FALCON_DIR)/bench_matvec_falcon_native.c

FALCON_REFERENCE = reference/falcon/main.c reference/falcon/ntt.c reference/falcon/poly.c
build/falcon-reference: $(FALCON_REFERENCE) | build
	$(CC) -O3 -std=c11 -D_POSIX_C_SOURCE=200809L -mavx2 -Ireference/falcon -o $@ $(FALCON_REFERENCE) -lm

KYBER_DIR = third_party/kyber512/NIST-PQ-Submission-Kyber-20201001/Additional_Implementations/avx2/crypto_kem/kyber512
KYBER_INCLUDES = -Ireference/kyber -I$(KYBER_DIR) -I$(KYBER_DIR)/keccak4x
KYBER_CORE = $(KYBER_DIR)/poly.c $(KYBER_DIR)/cbd.c $(KYBER_DIR)/symmetric-shake.c $(KYBER_DIR)/fips202.c $(KYBER_DIR)/fips202x4.c $(KYBER_DIR)/keccak4x/KeccakP-1600-times4-SIMD256.c $(KYBER_DIR)/consts.c $(KYBER_DIR)/verify.c $(KYBER_DIR)/ntt.S $(KYBER_DIR)/invntt.S $(KYBER_DIR)/basemul.S $(KYBER_DIR)/fq.S $(KYBER_DIR)/shuffle.S
KYBER_FLAGS = -O3 -std=c11 -mavx2 -mfma -mbmi2 -mpopcnt -maes
USE_KYBER_ASM ?= 1

ifeq ($(USE_KYBER_ASM),1)
KYBER_REFERENCE_ASM = reference/kyber/accum_u32_avx2.S reference/kyber/accum_r_u64_u32.S
KYBER_REFERENCE_FLAGS = -DKYBER_REFERENCE_USE_ASM=1
endif

build/kyber-reference: reference/kyber/main_kyber512.c reference/kyber/brlt_core.c reference/kyber/kyber_avx2_adapter.c $(KYBER_REFERENCE_ASM) $(KYBER_CORE) | build
	$(CC) $(KYBER_FLAGS) $(KYBER_REFERENCE_FLAGS) $(KYBER_INCLUDES) -o $@ reference/kyber/main_kyber512.c reference/kyber/brlt_core.c reference/kyber/kyber_avx2_adapter.c $(KYBER_REFERENCE_ASM) $(KYBER_CORE) -lm

# This is the unmodified native baseline selected by the original AC26 script.
build/kyber-baseline: baseline/kyber/bench_kyber_native_mv.c $(KYBER_CORE) | build
	$(CC) $(KYBER_FLAGS) -I$(KYBER_DIR) -I$(KYBER_DIR)/keccak4x -o $@ baseline/kyber/bench_kyber_native_mv.c $(KYBER_CORE) -lm

DILITHIUM_DIR = third_party/dilithium5
DILITHIUM_INCLUDES = -Ireference/dilithium -I$(DILITHIUM_DIR)
DILITHIUM_CORE = $(DILITHIUM_DIR)/consts.c $(DILITHIUM_DIR)/ntt.S $(DILITHIUM_DIR)/invntt.S $(DILITHIUM_DIR)/pointwise.S
DILITHIUM_COMMON = reference/dilithium/poly.c reference/dilithium/io.c reference/dilithium/table.c
DILITHIUM_FLAGS = -O3 -std=c11 -mavx2 -DDILITHIUM_MODE=5

build/dilithium-reference: reference/dilithium/main_dilithium.c reference/dilithium/brlt_core.c $(DILITHIUM_COMMON) $(DILITHIUM_CORE) | build
	$(CC) $(DILITHIUM_FLAGS) $(DILITHIUM_INCLUDES) -o $@ reference/dilithium/main_dilithium.c reference/dilithium/brlt_core.c $(DILITHIUM_COMMON) $(DILITHIUM_CORE)

# `main_baseline.c` is copied verbatim from the baseline invoked by AC26.
build/dilithium-baseline: baseline/dilithium/main_baseline.c $(DILITHIUM_COMMON) $(DILITHIUM_CORE) | build
	$(CC) $(DILITHIUM_FLAGS) $(DILITHIUM_INCLUDES) -o $@ baseline/dilithium/main_baseline.c $(DILITHIUM_COMMON) $(DILITHIUM_CORE)

build/dilithium5-gen: data/generator/dilithium5_gen.c $(DILITHIUM_COMMON) $(DILITHIUM_CORE) | build
	$(CC) $(DILITHIUM_FLAGS) $(DILITHIUM_INCLUDES) -o $@ data/generator/dilithium5_gen.c $(DILITHIUM_COMMON) $(DILITHIUM_CORE)

# The original Dilithium benchmark and reference code are parameter-generic;
# this target selects the upstream Dilithium-2 constants without changing the
# implementation or weakening its native baseline.
DILITHIUM2_FLAGS = -O3 -std=c11 -mavx2 -DDILITHIUM_MODE=2
build/dilithium2-reference: reference/dilithium/main_dilithium.c reference/dilithium/brlt_core.c $(DILITHIUM_COMMON) $(DILITHIUM_CORE) | build
	$(CC) $(DILITHIUM2_FLAGS) $(DILITHIUM_INCLUDES) -o $@ reference/dilithium/main_dilithium.c reference/dilithium/brlt_core.c $(DILITHIUM_COMMON) $(DILITHIUM_CORE)

build/dilithium2-baseline: baseline/dilithium/main_baseline.c $(DILITHIUM_COMMON) $(DILITHIUM_CORE) | build
	$(CC) $(DILITHIUM2_FLAGS) $(DILITHIUM_INCLUDES) -o $@ baseline/dilithium/main_baseline.c $(DILITHIUM_COMMON) $(DILITHIUM_CORE)

build/dilithium2-gen: data/generator/dilithium5_gen.c $(DILITHIUM_COMMON) $(DILITHIUM_CORE) | build
	$(CC) $(DILITHIUM2_FLAGS) $(DILITHIUM_INCLUDES) -o $@ data/generator/dilithium5_gen.c $(DILITHIUM_COMMON) $(DILITHIUM_CORE)

RPT_FLAGS = -O3 -std=c11 -mavx2
build/rpt-reference: reference/rpt/rpt_verify.c | build
	$(CC) $(RPT_FLAGS) -DRPT_REFERENCE=1 -o $@ reference/rpt/rpt_verify.c

build/rpt-baseline: reference/rpt/rpt_verify.c | build
	$(CC) $(RPT_FLAGS) -DRPT_REFERENCE=0 -o $@ reference/rpt/rpt_verify.c

test: all
	./tests/test.sh

clean:
	rm -rf build data/generated results
