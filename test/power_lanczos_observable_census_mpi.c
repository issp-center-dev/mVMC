#include <limits.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "power_lanczos_observable_census.h"

static int WriteObservable(const char *path) {
  static const char bytes[] =
      "====================\n"
      "NCisAjs 3\n"
      "====================\n"
      "observable rows\n"
      "====================\n"
      "0 0 0 0\n"
      "0 0 1 0\n"
      "1 0 0 0\n";
  FILE *stream = fopen(path, "wb");
  if (stream == NULL) return 0;
  if (fwrite(bytes, 1, sizeof(bytes) - 1, stream) != sizeof(bytes) - 1) {
    fclose(stream);
    return 0;
  }
  return fclose(stream) == 0;
}

int main(int argc, char **argv) {
  MVMCPowerLanczosObservablePlan root_plan;
  MVMCPowerLanczosObservablePlan received_plan;
  MVMCPowerLanczosObservableCensusStatus status =
      MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK;
  unsigned char *wire = NULL;
  unsigned long long wire_size_wire = 0;
  size_t wire_size = 0;
  size_t packed_size = 0;
  char diagnostic[512] = {0};
  char expected_sha[65] = {0};
  char expected_semantic_sha[65] = {0};
  char input_path[256] = {0};
  int rank = 0;
  int world_size = 0;
  int status_code = 0;
  int local_failure = 0;
  int global_failure = 0;
  (void)argv;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  mvmc_power_lanczos_observable_plan_init(&root_plan);
  mvmc_power_lanczos_observable_plan_init(&received_plan);

  if (rank == 0) {
    const char *paths[MVMC_POWER_LANCZOS_OBSERVABLE_FAMILY_COUNT] = {
        input_path, NULL, NULL};
    snprintf(input_path, sizeof(input_path),
             "power_lanczos_observable_census_mpi.%ld.def", (long)getpid());
    if (!WriteObservable(input_path)) {
      status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_IO_ERROR;
    } else {
      status = mvmc_power_lanczos_observable_plan_build_from_files(
          4, 2, paths, &root_plan, diagnostic, sizeof(diagnostic));
    }
    if (status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) {
      status = mvmc_power_lanczos_observable_plan_wire_size(&root_plan,
                                                           &wire_size);
    }
    if (status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK &&
        wire_size > (size_t)INT_MAX) {
      status = MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_SIZE_OVERFLOW;
    }
    if (status == MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK) {
      strcpy(expected_sha, root_plan.raw_observable_census_sha256);
      strcpy(expected_semantic_sha,
             root_plan.semantic_observable_census_sha256);
      wire_size_wire = (unsigned long long)wire_size;
    }
    status_code = (int)status;
  }
  MPI_Bcast(&status_code, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&wire_size_wire, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(expected_sha, (int)sizeof(expected_sha), MPI_CHAR, 0,
            MPI_COMM_WORLD);
  MPI_Bcast(expected_semantic_sha, (int)sizeof(expected_semantic_sha),
            MPI_CHAR, 0, MPI_COMM_WORLD);
  if (status_code != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK ||
      wire_size_wire == 0 || wire_size_wire > (unsigned long long)INT_MAX) {
    local_failure = 1;
  } else {
    wire_size = (size_t)wire_size_wire;
    wire = (unsigned char *)malloc(wire_size);
    if (wire == NULL) local_failure = 1;
  }
  MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  if (global_failure == 0) {
    if (rank == 0) {
      status = mvmc_power_lanczos_observable_plan_pack(
          &root_plan, wire, wire_size, &packed_size);
      if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK ||
          packed_size != wire_size) {
        local_failure = 1;
      }
    }
    MPI_Bcast(&local_failure, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (local_failure == 0) {
      MPI_Bcast(wire, (int)wire_size, MPI_BYTE, 0, MPI_COMM_WORLD);
      status = mvmc_power_lanczos_observable_plan_unpack(
          wire, wire_size, &received_plan, diagnostic, sizeof(diagnostic));
      if (status != MVMC_POWER_LANCZOS_OBSERVABLE_CENSUS_OK ||
          strcmp(received_plan.raw_observable_census_sha256, expected_sha) !=
              0 ||
          strcmp(received_plan.semantic_observable_census_sha256,
                 expected_semantic_sha) != 0 ||
          received_plan.record_count != 3 || received_plan.records == NULL ||
          received_plan.records[1].adjoint_class !=
              MVMC_POWER_LANCZOS_OBSERVABLE_ADJOINT_REQUESTED) {
        local_failure = 1;
      }
    }
  }
  MPI_Allreduce(&local_failure, &global_failure, 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  if (rank == 0 && global_failure != 0) {
    fprintf(stderr,
            "PowerLanczosObservableCensus_MPI FAIL: status=%s diagnostic=%s\n",
            mvmc_power_lanczos_observable_census_status_string(status),
            diagnostic);
  }

  free(wire);
  mvmc_power_lanczos_observable_plan_destroy(&received_plan);
  mvmc_power_lanczos_observable_plan_destroy(&root_plan);
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0 && input_path[0] != '\0') unlink(input_path);
  if (rank == 0 && global_failure == 0) {
    printf("PowerLanczosObservableCensus_MPI: PASS world_size=%d raw_sha256=%s "
           "semantic_sha256=%s\n",
           world_size, expected_sha, expected_semantic_sha);
  }
  MPI_Finalize();
  return global_failure == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
