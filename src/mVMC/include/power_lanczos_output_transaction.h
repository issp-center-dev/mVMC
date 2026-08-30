#ifndef MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_H
#define MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_H

#if defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)

#include "power_lanczos_sha256.h"

#include <stddef.h>
#include <stdint.h>

#define MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_VERSION UINT64_C(1)
#define MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_MAX_FILES 4
#define MVMC_POWER_LANCZOS_OUTPUT_BASENAME_CAPACITY 192

typedef enum {
  MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_OK = 0,
  MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_INVALID_ARGUMENT,
  MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_PATH_REJECTED,
  MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_TARGET_EXISTS,
  MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_IO_FAILURE,
  MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_VERIFY_FAILURE,
  MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_COMMIT_FAILURE,
  MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_ROLLBACK_FAILURE
} MVMCPowerLanczosOutputTransactionStatus;

typedef struct {
  const char *target_basename;
  const unsigned char *contents;
  size_t content_size;
} MVMCPowerLanczosOutputFile;

#if defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
typedef enum {
  MVMC_POWER_LANCZOS_OUTPUT_FAULT_NONE = 0,
  MVMC_POWER_LANCZOS_OUTPUT_FAULT_WRITE,
  MVMC_POWER_LANCZOS_OUTPUT_FAULT_FILE_FSYNC,
  MVMC_POWER_LANCZOS_OUTPUT_FAULT_VERIFY,
  MVMC_POWER_LANCZOS_OUTPUT_FAULT_CLOSE,
  MVMC_POWER_LANCZOS_OUTPUT_FAULT_PUBLISH,
  MVMC_POWER_LANCZOS_OUTPUT_FAULT_DIRECTORY_FSYNC
} MVMCPowerLanczosOutputFault;
#endif

typedef struct {
  const char *directory;
  const MVMCPowerLanczosOutputFile *files;
  size_t file_count;
  size_t manifest_index;
#if defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
  MVMCPowerLanczosOutputFault injected_fault;
  size_t injected_fault_file_index;
#endif
} MVMCPowerLanczosOutputTransactionConfig;

typedef struct {
  int valid;
  MVMCPowerLanczosOutputTransactionStatus status;
  uint64_t version;
  size_t file_count;
  size_t manifest_index;
  size_t failed_file_index;
  int committed;
  int rollback_complete;
  char sha256[MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_MAX_FILES]
             [MVMC_POWER_LANCZOS_SHA256_HEX_CAPACITY];
} MVMCPowerLanczosOutputTransactionResult;

int mvmc_power_lanczos_output_basename_valid(const char *basename);

MVMCPowerLanczosOutputTransactionStatus
mvmc_power_lanczos_output_transaction_commit(
    const MVMCPowerLanczosOutputTransactionConfig *config,
    MVMCPowerLanczosOutputTransactionResult *result);

const char *mvmc_power_lanczos_output_transaction_status_string(
    MVMCPowerLanczosOutputTransactionStatus status);

#endif /* bounded engine */

#endif /* MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_H */
