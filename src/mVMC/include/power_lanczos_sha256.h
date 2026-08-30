#ifndef MVMC_POWER_LANCZOS_SHA256_H
#define MVMC_POWER_LANCZOS_SHA256_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include <stddef.h>

#define MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES 32
#define MVMC_POWER_LANCZOS_SHA256_HEX_CAPACITY 65

int mvmc_power_lanczos_sha256_bytes(
    const void *input, size_t input_size,
    unsigned char digest[MVMC_POWER_LANCZOS_SHA256_DIGEST_BYTES]);

int mvmc_power_lanczos_sha256_hex(
    const void *input, size_t input_size,
    char hex[MVMC_POWER_LANCZOS_SHA256_HEX_CAPACITY]);

#endif

#endif
