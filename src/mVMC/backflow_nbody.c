#include "./include/backflow_nbody.h"

#include "./include/global.h"
#include "./include/locgrn.h"
#include "./include/locgrn_fsz.h"
#include "./include/matrix.h"
#include "./include/projection.h"
#include "./include/sector_projection.h"
#include "./include/slater.h"

#include <string.h>

int GetBFNBodyScratchSizes(int maxOrder, int useFsz,
                           BFNBodyScratchSizes *sizes) {
  BFNBodyScratchDimensions dimensions;

  if(sizes == NULL || maxOrder < 1 || (useFsz != 0 && useFsz != 1)) {
    return BF_NBODY_WORKSPACE_ERROR;
  }

  memset(&dimensions, 0, sizeof(dimensions));
  dimensions.maxOrder = maxOrder;
  dimensions.useFsz = useFsz;
  dimensions.nsize = Nsize;
  dimensions.nsite2 = Nsite2;
  dimensions.nproj = NProj;
  dimensions.nsite = Nsite;
  dimensions.nrange = Nrange;
  dimensions.nqpFull = NQPFull;
  dimensions.lapackLWork = LapackLWork;

  if(useFsz
     && (GetGreenFuncBF_fsz_buffer_work_size(
             &dimensions.bfFszGreenBufferCount,
             &dimensions.pfUpdateIntCount,
             &dimensions.pfUpdateDoubleCount) != 0
         || GetSlaterElmBF_fsz_hop_int_work_size(
             &dimensions.bfHopIntSize) != 0)) {
    return BF_NBODY_WORKSPACE_ERROR;
  }

  return ComputeBFNBodyScratchSizes(&dimensions, sizes);
}
