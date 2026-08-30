#include "power_lanczos_output_transaction.h"

#if !defined(MVMC_ENABLE_POWER_LANCZOS_BOUNDED_ENGINE)
#error "power_lanczos_output_transaction.c requires bounded Krylov"
#endif

#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#ifndef O_CLOEXEC
#define O_CLOEXEC 0
#endif

#ifndef O_NOFOLLOW
#define O_NOFOLLOW 0
#endif

typedef struct {
  char temporary[MVMC_POWER_LANCZOS_OUTPUT_BASENAME_CAPACITY];
  int descriptor;
  int temporary_exists;
  int target_exists;
} StagedFile;

int mvmc_power_lanczos_output_basename_valid(const char *basename) {
  size_t length;
  size_t index;
  if (basename == NULL) return 0;
  length = strlen(basename);
  if (length == 0 || length >= MVMC_POWER_LANCZOS_OUTPUT_BASENAME_CAPACITY ||
      basename[0] == '.' || strstr(basename, "..") != NULL) {
    return 0;
  }
  for (index = 0; index < length; ++index) {
    const unsigned char value = (unsigned char)basename[index];
    const int alpha_numeric =
        (value >= 'a' && value <= 'z') ||
        (value >= 'A' && value <= 'Z') ||
        (value >= '0' && value <= '9');
    if (!alpha_numeric && value != '_' && value != '-' && value != '.') {
      return 0;
    }
  }
  return 1;
}

#if defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
static int injected_file_fault(
    const MVMCPowerLanczosOutputTransactionConfig *config,
    MVMCPowerLanczosOutputFault fault,
    size_t file_index) {
  return config->injected_fault == fault &&
         config->injected_fault_file_index == file_index;
}
#endif

static int write_all(
    int descriptor, const unsigned char *contents, size_t content_size,
    const MVMCPowerLanczosOutputTransactionConfig *config,
    size_t file_index) {
  size_t offset = 0;
#if !defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
  (void)config;
  (void)file_index;
#endif
  while (offset < content_size) {
    size_t request = content_size - offset;
    ssize_t written;
    if (request > (size_t)SSIZE_MAX) request = (size_t)SSIZE_MAX;
#if defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
    if (injected_file_fault(config, MVMC_POWER_LANCZOS_OUTPUT_FAULT_WRITE,
                            file_index)) {
      if (offset == 0 && request > 1) {
        written = write(descriptor, contents, request / 2);
        if (written > 0) offset += (size_t)written;
      }
      errno = EIO;
      return 0;
    }
#endif
    written = write(descriptor, contents + offset, request);
    if (written < 0 && errno == EINTR) continue;
    if (written <= 0) return 0;
    offset += (size_t)written;
  }
  return 1;
}

static int verify_contents(int descriptor, const unsigned char *contents,
                           size_t content_size,
                           const MVMCPowerLanczosOutputTransactionConfig *config,
                           size_t file_index) {
  unsigned char buffer[4096];
  size_t offset = 0;
#if defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
  if (injected_file_fault(config, MVMC_POWER_LANCZOS_OUTPUT_FAULT_VERIFY,
                          file_index)) {
    return 0;
  }
#else
  (void)config;
  (void)file_index;
#endif
  while (offset < content_size) {
    const size_t remaining = content_size - offset;
    const size_t request = remaining < sizeof(buffer) ? remaining
                                                       : sizeof(buffer);
    ssize_t received = pread(descriptor, buffer, request, (off_t)offset);
    if (received < 0 && errno == EINTR) continue;
    if (received <= 0 || (size_t)received > request ||
        memcmp(buffer, contents + offset, (size_t)received) != 0) {
      return 0;
    }
    offset += (size_t)received;
  }
  {
    unsigned char extra;
    ssize_t received;
    do {
      received = pread(descriptor, &extra, 1, (off_t)content_size);
    } while (received < 0 && errno == EINTR);
    if (received != 0) return 0;
  }
  return 1;
}

static int create_temporary(int directory_descriptor, size_t file_index,
                            StagedFile *staged) {
  unsigned int attempt;
  const long process_id = (long)getpid();
  for (attempt = 0; attempt < 1024; ++attempt) {
    const int written = snprintf(
        staged->temporary, sizeof(staged->temporary),
        ".mvmc-pl-%ld-%zu-%u.tmp", process_id, file_index, attempt);
    if (written <= 0 || (size_t)written >= sizeof(staged->temporary)) {
      return 0;
    }
    staged->descriptor = openat(
        directory_descriptor, staged->temporary,
        O_CREAT | O_EXCL | O_RDWR | O_CLOEXEC | O_NOFOLLOW, S_IRUSR | S_IWUSR);
    if (staged->descriptor >= 0) {
      staged->temporary_exists = 1;
      return 1;
    }
    if (errno != EEXIST) return 0;
  }
  errno = EEXIST;
  return 0;
}

static int target_state(int directory_descriptor, const char *basename) {
  struct stat status;
  if (fstatat(directory_descriptor, basename, &status,
              AT_SYMLINK_NOFOLLOW) == 0) {
    return 0;
  }
  return errno == ENOENT ? 1 : -1;
}

static int publish_file(
    int directory_descriptor, StagedFile *staged, const char *target,
    const MVMCPowerLanczosOutputTransactionConfig *config,
    size_t file_index) {
#if defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
  if (injected_file_fault(config, MVMC_POWER_LANCZOS_OUTPUT_FAULT_PUBLISH,
                          file_index)) {
    errno = EIO;
    return 0;
  }
#else
  (void)config;
  (void)file_index;
#endif
  if (linkat(directory_descriptor, staged->temporary, directory_descriptor,
             target,
             0) != 0) {
    return 0;
  }
  staged->target_exists = 1;
  if (unlinkat(directory_descriptor, staged->temporary, 0) != 0) {
    return 0;
  }
  staged->temporary_exists = 0;
  return 1;
}

static int close_staged(
    StagedFile *staged,
    const MVMCPowerLanczosOutputTransactionConfig *config,
    size_t file_index) {
  int close_status;
  if (staged->descriptor < 0) return 1;
  close_status = close(staged->descriptor);
  staged->descriptor = -1;
#if defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
  if (injected_file_fault(config, MVMC_POWER_LANCZOS_OUTPUT_FAULT_CLOSE,
                          file_index)) {
    errno = EIO;
    return 0;
  }
#else
  (void)config;
  (void)file_index;
#endif
  return close_status == 0;
}

MVMCPowerLanczosOutputTransactionStatus
mvmc_power_lanczos_output_transaction_commit(
    const MVMCPowerLanczosOutputTransactionConfig *config,
    MVMCPowerLanczosOutputTransactionResult *result) {
  MVMCPowerLanczosOutputTransactionResult candidate;
  StagedFile staged[MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_MAX_FILES];
  int directory_descriptor = -1;
  MVMCPowerLanczosOutputTransactionStatus status =
      MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_OK;
  size_t file_index;
  size_t commit_ordinal;
  int rollback_ok = 1;
  if (result == NULL) {
    return MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_INVALID_ARGUMENT;
  }
  memset(result, 0, sizeof(*result));
  result->status = MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_INVALID_ARGUMENT;
  result->failed_file_index = SIZE_MAX;
  memset(&candidate, 0, sizeof(candidate));
  candidate.valid = 1;
  candidate.version = MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_VERSION;
  candidate.failed_file_index = SIZE_MAX;
  memset(staged, 0, sizeof(staged));
  for (file_index = 0;
       file_index < MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_MAX_FILES;
       ++file_index) {
    staged[file_index].descriptor = -1;
  }
  if (config == NULL || config->directory == NULL ||
      config->files == NULL || config->file_count == 0 ||
      config->file_count > MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_MAX_FILES ||
      config->manifest_index >= config->file_count) {
    return MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_INVALID_ARGUMENT;
  }
  candidate.file_count = config->file_count;
  candidate.manifest_index = config->manifest_index;
  for (file_index = 0; file_index < config->file_count; ++file_index) {
    size_t other;
    const MVMCPowerLanczosOutputFile *file = &config->files[file_index];
    if (!mvmc_power_lanczos_output_basename_valid(file->target_basename)) {
      candidate.status = MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_PATH_REJECTED;
      candidate.failed_file_index = file_index;
      *result = candidate;
      return candidate.status;
    }
    if (file->contents == NULL || file->content_size == 0 ||
        file->content_size > (size_t)INT64_MAX) {
      candidate.status = MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_INVALID_ARGUMENT;
      candidate.failed_file_index = file_index;
      *result = candidate;
      return candidate.status;
    }
    for (other = 0; other < file_index; ++other) {
      if (strcmp(file->target_basename,
                 config->files[other].target_basename) == 0) {
        candidate.status =
            MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_PATH_REJECTED;
        candidate.failed_file_index = file_index;
        *result = candidate;
        return candidate.status;
      }
    }
  }
  directory_descriptor = open(config->directory,
                              O_RDONLY | O_CLOEXEC | O_NOFOLLOW | O_DIRECTORY);
  if (directory_descriptor < 0) {
    candidate.status = MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_IO_FAILURE;
    *result = candidate;
    return candidate.status;
  }
  for (file_index = 0; file_index < config->file_count; ++file_index) {
    const int target = target_state(
        directory_descriptor, config->files[file_index].target_basename);
    if (target != 1) {
      status = target == 0
                   ? MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_TARGET_EXISTS
                   : MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_IO_FAILURE;
      candidate.failed_file_index = file_index;
      goto rollback;
    }
  }
  for (file_index = 0; file_index < config->file_count; ++file_index) {
    const MVMCPowerLanczosOutputFile *file = &config->files[file_index];
    if (!create_temporary(directory_descriptor, file_index,
                          &staged[file_index]) ||
        !write_all(staged[file_index].descriptor, file->contents,
                   file->content_size, config, file_index)) {
      status = MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_IO_FAILURE;
      candidate.failed_file_index = file_index;
      goto rollback;
    }
#if defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
    if (injected_file_fault(config,
                            MVMC_POWER_LANCZOS_OUTPUT_FAULT_FILE_FSYNC,
                            file_index)) {
      errno = EIO;
      status = MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_IO_FAILURE;
      candidate.failed_file_index = file_index;
      goto rollback;
    }
#endif
    if (fsync(staged[file_index].descriptor) != 0) {
      status = MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_IO_FAILURE;
      candidate.failed_file_index = file_index;
      goto rollback;
    }
    if (!verify_contents(staged[file_index].descriptor, file->contents,
                         file->content_size, config, file_index)) {
      status = MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_VERIFY_FAILURE;
      candidate.failed_file_index = file_index;
      goto rollback;
    }
    if (!mvmc_power_lanczos_sha256_hex(file->contents, file->content_size,
                                        candidate.sha256[file_index])) {
      status = MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_VERIFY_FAILURE;
      candidate.failed_file_index = file_index;
      goto rollback;
    }
    if (!close_staged(&staged[file_index], config, file_index)) {
      status = MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_IO_FAILURE;
      candidate.failed_file_index = file_index;
      goto rollback;
    }
  }
  for (commit_ordinal = 0; commit_ordinal < config->file_count;
       ++commit_ordinal) {
    file_index = commit_ordinal == config->file_count - 1
                     ? config->manifest_index
                     : commit_ordinal < config->manifest_index
                           ? commit_ordinal
                           : commit_ordinal + 1;
    if (!publish_file(directory_descriptor, &staged[file_index],
                      config->files[file_index].target_basename, config,
                      file_index)) {
      status = MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_COMMIT_FAILURE;
      candidate.failed_file_index = file_index;
      goto rollback;
    }
  }
#if defined(MVMC_ENABLE_POWER_LANCZOS_P6_TESTING)
  if (config->injected_fault ==
      MVMC_POWER_LANCZOS_OUTPUT_FAULT_DIRECTORY_FSYNC) {
    errno = EIO;
    status = MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_COMMIT_FAILURE;
    goto rollback;
  }
#endif
  if (fsync(directory_descriptor) != 0) {
    status = MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_COMMIT_FAILURE;
    goto rollback;
  }
  (void)close(directory_descriptor);
  candidate.status = MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_OK;
  candidate.committed = 1;
  candidate.rollback_complete = 1;
  *result = candidate;
  return candidate.status;

rollback:
  for (file_index = config->file_count; file_index > 0; --file_index) {
    const size_t index = file_index - 1;
    if (!close_staged(&staged[index], config, index)) rollback_ok = 0;
    if (staged[index].temporary_exists &&
        unlinkat(directory_descriptor, staged[index].temporary, 0) != 0 &&
        errno != ENOENT) {
      rollback_ok = 0;
    }
    if (staged[index].target_exists &&
        unlinkat(directory_descriptor,
                 config->files[index].target_basename, 0) != 0 &&
        errno != ENOENT) {
      rollback_ok = 0;
    }
  }
  if (fsync(directory_descriptor) != 0) rollback_ok = 0;
  if (close(directory_descriptor) != 0) rollback_ok = 0;
  candidate.status = rollback_ok
                         ? status
                         : MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_ROLLBACK_FAILURE;
  candidate.committed = 0;
  candidate.rollback_complete = rollback_ok;
  *result = candidate;
  return candidate.status;
}

const char *mvmc_power_lanczos_output_transaction_status_string(
    MVMCPowerLanczosOutputTransactionStatus status) {
  switch (status) {
    case MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_OK:
      return "ok";
    case MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_INVALID_ARGUMENT:
      return "invalid_argument";
    case MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_PATH_REJECTED:
      return "path_rejected";
    case MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_TARGET_EXISTS:
      return "target_exists";
    case MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_IO_FAILURE:
      return "io_failure";
    case MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_VERIFY_FAILURE:
      return "verify_failure";
    case MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_COMMIT_FAILURE:
      return "commit_failure";
    case MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_ROLLBACK_FAILURE:
      return "rollback_failure";
    default:
      return "invalid";
  }
}
