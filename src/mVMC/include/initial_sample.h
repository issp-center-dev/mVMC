#ifndef MVMC_INITIAL_SAMPLE_H
#define MVMC_INITIAL_SAMPLE_H

int InitialSamplePlaceLocalSpin(int *eleIdx, int *eleCfg, int *eleSpn);
int InitialSamplePlaceLocalSpinTry(int *eleIdx, int *eleCfg, int *eleSpn,
                                   const int nTryMax);

#endif /* MVMC_INITIAL_SAMPLE_H */
