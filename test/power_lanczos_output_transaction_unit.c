#include "power_lanczos_output_transaction.h"

#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

static int failures = 0;

#define CHECK(condition, ...)                                               \
  do {                                                                      \
    if (!(condition)) {                                                     \
      fprintf(stderr, "PowerLanczosOutputTransaction_Unit FAIL: ");       \
      fprintf(stderr, __VA_ARGS__);                                         \
      fprintf(stderr, " (line %d)\n", __LINE__);                         \
      ++failures;                                                           \
    }                                                                       \
  } while (0)

static int directory_entry_count(const char *directory) {
  DIR *stream = opendir(directory);
  struct dirent *entry;
  int count = 0;
  if (stream == NULL) return -1;
  while ((entry = readdir(stream)) != NULL) {
    if (strcmp(entry->d_name, ".") != 0 &&
        strcmp(entry->d_name, "..") != 0) {
      ++count;
    }
  }
  if (closedir(stream) != 0) return -1;
  return count;
}

static int read_matches(const char *directory, const char *basename,
                        const unsigned char *expected, size_t size) {
  char path[512];
  unsigned char buffer[256];
  FILE *stream;
  size_t received;
  int trailing;
  if (size > sizeof(buffer) ||
      snprintf(path, sizeof(path), "%s/%s", directory, basename) <= 0) {
    return 0;
  }
  stream = fopen(path, "rb");
  if (stream == NULL) return 0;
  received = fread(buffer, 1, sizeof(buffer), stream);
  trailing = ferror(stream);
  if (fclose(stream) != 0) return 0;
  return !trailing && received == size &&
         memcmp(buffer, expected, size) == 0;
}

static int remove_target(const char *directory, const char *basename) {
  char path[512];
  if (snprintf(path, sizeof(path), "%s/%s", directory, basename) <= 0) {
    return 0;
  }
  return unlink(path) == 0;
}

static MVMCPowerLanczosOutputTransactionConfig fixture_config(
    const char *directory, MVMCPowerLanczosOutputFile files[4]) {
  static const unsigned char numeric[] = "v4 1.25 -0.5\n";
  static const unsigned char matrix[] = "S 0 0 1 0\nK 0 0 2 0\n";
  static const unsigned char support[] = "support finite\n";
  static const unsigned char metadata[] = "{\"schema_version\":4}\n";
  MVMCPowerLanczosOutputTransactionConfig config;
  memset(files, 0, 4 * sizeof(*files));
  files[0].target_basename = "case_pl_out_7.dat";
  files[0].contents = numeric;
  files[0].content_size = sizeof(numeric) - 1;
  files[1].target_basename = "case_pl_matrix_7.dat";
  files[1].contents = matrix;
  files[1].content_size = sizeof(matrix) - 1;
  files[2].target_basename = "case_pl_support_7.dat";
  files[2].contents = support;
  files[2].content_size = sizeof(support) - 1;
  files[3].target_basename = "case_pl_meta_7.json";
  files[3].contents = metadata;
  files[3].content_size = sizeof(metadata) - 1;
  memset(&config, 0, sizeof(config));
  config.directory = directory;
  config.files = files;
  config.file_count = 4;
  config.manifest_index = 3;
  return config;
}

static void check_no_fixture_targets(const char *directory,
                                     const MVMCPowerLanczosOutputFile *files) {
  size_t index;
  for (index = 0; index < 4; ++index) {
    char path[512];
    struct stat status;
    (void)snprintf(path, sizeof(path), "%s/%s", directory,
                   files[index].target_basename);
    CHECK(lstat(path, &status) != 0 && errno == ENOENT,
          "rollback left target %s", files[index].target_basename);
  }
  CHECK(directory_entry_count(directory) == 0,
        "rollback left temporary entries");
}

static void test_basename_grammar(void) {
  static const char *invalid[] = {
      "", ".hidden", "../escape", "a..b", "a/b", "a\\b",
      "a:b", "file://x", "line\nname", "space name", NULL};
  size_t index;
  CHECK(mvmc_power_lanczos_output_basename_valid("x_pl_meta_0.json"),
        "valid basename rejected");
  CHECK(mvmc_power_lanczos_output_basename_valid("A-1_b.2.dat"),
        "valid grammar rejected");
  for (index = 0; invalid[index] != NULL; ++index) {
    CHECK(!mvmc_power_lanczos_output_basename_valid(invalid[index]),
          "invalid basename accepted: %s", invalid[index]);
  }
}

static void test_success_and_precondition(const char *directory) {
  MVMCPowerLanczosOutputFile files[4];
  MVMCPowerLanczosOutputTransactionConfig config =
      fixture_config(directory, files);
  MVMCPowerLanczosOutputTransactionResult result;
  size_t index;
  CHECK(mvmc_power_lanczos_output_transaction_commit(&config, &result) ==
                MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_OK &&
            result.valid && result.committed && result.rollback_complete &&
            result.version ==
                MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_VERSION &&
            result.file_count == 4 && result.manifest_index == 3 &&
            result.failed_file_index == SIZE_MAX,
        "four-file transaction status=%s",
        mvmc_power_lanczos_output_transaction_status_string(result.status));
  CHECK(directory_entry_count(directory) == 4,
        "successful transaction left temporary entries");
  for (index = 0; index < 4; ++index) {
    CHECK(read_matches(directory, files[index].target_basename,
                       files[index].contents, files[index].content_size),
          "published payload mismatch %zu", index);
    CHECK(strlen(result.sha256[index]) == 64,
          "published SHA missing %zu", index);
  }
  CHECK(mvmc_power_lanczos_output_transaction_commit(&config, &result) ==
                MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_TARGET_EXISTS &&
            !result.committed && result.rollback_complete,
        "existing public target precondition");
  CHECK(directory_entry_count(directory) == 4,
        "precondition failure changed public bundle");
  for (index = 0; index < 4; ++index) {
    CHECK(remove_target(directory, files[index].target_basename),
          "fixture cleanup %zu", index);
  }
}

static void test_metadata_only(const char *directory) {
  static const unsigned char metadata[] = "{\"failure\":true}\n";
  MVMCPowerLanczosOutputFile file;
  MVMCPowerLanczosOutputTransactionConfig config;
  MVMCPowerLanczosOutputTransactionResult result;
  memset(&file, 0, sizeof(file));
  file.target_basename = "case_pl_meta_8.json";
  file.contents = metadata;
  file.content_size = sizeof(metadata) - 1;
  memset(&config, 0, sizeof(config));
  config.directory = directory;
  config.files = &file;
  config.file_count = 1;
  config.manifest_index = 0;
  CHECK(mvmc_power_lanczos_output_transaction_commit(&config, &result) ==
                MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_OK &&
            result.committed && directory_entry_count(directory) == 1 &&
            read_matches(directory, file.target_basename, file.contents,
                         file.content_size),
        "metadata-only transaction");
  CHECK(remove_target(directory, file.target_basename),
        "metadata-only cleanup");
}

static void test_fault_rollback(const char *directory) {
  static const MVMCPowerLanczosOutputFault faults[] = {
      MVMC_POWER_LANCZOS_OUTPUT_FAULT_WRITE,
      MVMC_POWER_LANCZOS_OUTPUT_FAULT_FILE_FSYNC,
      MVMC_POWER_LANCZOS_OUTPUT_FAULT_VERIFY,
      MVMC_POWER_LANCZOS_OUTPUT_FAULT_CLOSE,
      MVMC_POWER_LANCZOS_OUTPUT_FAULT_PUBLISH};
  MVMCPowerLanczosOutputFile files[4];
  size_t fault_index;
  for (fault_index = 0; fault_index < sizeof(faults) / sizeof(faults[0]);
       ++fault_index) {
    size_t file_index;
    for (file_index = 0; file_index < 4; ++file_index) {
      MVMCPowerLanczosOutputTransactionConfig config =
          fixture_config(directory, files);
      MVMCPowerLanczosOutputTransactionResult result;
      config.injected_fault = faults[fault_index];
      config.injected_fault_file_index = file_index;
      CHECK(mvmc_power_lanczos_output_transaction_commit(&config, &result) !=
                    MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_OK &&
                result.valid && !result.committed &&
                result.rollback_complete &&
                result.failed_file_index == file_index,
            "fault %zu file %zu status=%s rollback=%d", fault_index,
            file_index,
            mvmc_power_lanczos_output_transaction_status_string(result.status),
            result.rollback_complete);
      check_no_fixture_targets(directory, files);
    }
  }
  {
    MVMCPowerLanczosOutputTransactionConfig config =
        fixture_config(directory, files);
    MVMCPowerLanczosOutputTransactionResult result;
    config.injected_fault =
        MVMC_POWER_LANCZOS_OUTPUT_FAULT_DIRECTORY_FSYNC;
    CHECK(mvmc_power_lanczos_output_transaction_commit(&config, &result) ==
                  MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_COMMIT_FAILURE &&
              result.valid && !result.committed && result.rollback_complete &&
              result.failed_file_index == SIZE_MAX,
          "directory fsync fault status=%s rollback=%d",
          mvmc_power_lanczos_output_transaction_status_string(result.status),
          result.rollback_complete);
    check_no_fixture_targets(directory, files);
  }
}

static void test_no_follow_boundaries(const char *directory) {
  MVMCPowerLanczosOutputFile files[4];
  MVMCPowerLanczosOutputTransactionConfig config =
      fixture_config(directory, files);
  MVMCPowerLanczosOutputTransactionResult result;
  char target_path[512];
  char directory_link[512];
  CHECK(snprintf(target_path, sizeof(target_path), "%s/%s", directory,
                 files[0].target_basename) > 0,
        "target symlink path");
  CHECK(symlink("sentinel-outside-transaction", target_path) == 0,
        "target symlink fixture");
  CHECK(mvmc_power_lanczos_output_transaction_commit(&config, &result) ==
                MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_TARGET_EXISTS &&
            result.rollback_complete,
        "target symlink no-clobber");
  CHECK(unlink(target_path) == 0, "target symlink cleanup");

  CHECK(snprintf(directory_link, sizeof(directory_link), "%s-link",
                 directory) > 0,
        "directory symlink path");
  CHECK(symlink(directory, directory_link) == 0,
        "directory symlink fixture");
  config.directory = directory_link;
  CHECK(mvmc_power_lanczos_output_transaction_commit(&config, &result) ==
                MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_IO_FAILURE &&
            !result.committed,
        "directory symlink rejected");
  CHECK(unlink(directory_link) == 0, "directory symlink cleanup");
  CHECK(directory_entry_count(directory) == 0,
        "no-follow tests changed output directory");
}

static void test_invalid_config(const char *directory) {
  MVMCPowerLanczosOutputFile files[4];
  MVMCPowerLanczosOutputTransactionConfig config =
      fixture_config(directory, files);
  MVMCPowerLanczosOutputTransactionResult result;
  files[1].target_basename = "../escape";
  CHECK(mvmc_power_lanczos_output_transaction_commit(&config, &result) ==
            MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_PATH_REJECTED,
        "path mutation rejected before open");
  files[1].target_basename = files[0].target_basename;
  CHECK(mvmc_power_lanczos_output_transaction_commit(&config, &result) ==
            MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_PATH_REJECTED,
        "duplicate target rejected");
  config = fixture_config(directory, files);
  files[0].content_size = 0;
  CHECK(mvmc_power_lanczos_output_transaction_commit(&config, &result) ==
            MVMC_POWER_LANCZOS_OUTPUT_TRANSACTION_INVALID_ARGUMENT,
        "empty payload rejected");
  CHECK(directory_entry_count(directory) == 0,
        "invalid config opened output files");
}

int main(void) {
  char directory[] = "./power_lanczos_output_transaction_XXXXXX";
  if (mkdtemp(directory) == NULL) {
    perror("mkdtemp");
    return 1;
  }
  test_basename_grammar();
  test_success_and_precondition(directory);
  test_metadata_only(directory);
  test_fault_rollback(directory);
  test_no_follow_boundaries(directory);
  test_invalid_config(directory);
  CHECK(directory_entry_count(directory) == 0, "final directory not empty");
  CHECK(rmdir(directory) == 0, "fixture directory cleanup");
  if (failures != 0) {
    fprintf(stderr, "%d output transaction checks failed\n", failures);
    return 1;
  }
  puts("power-Lanczos output transaction checks passed");
  return 0;
}
