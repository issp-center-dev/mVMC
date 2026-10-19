#ifndef MVMC_BACKFLOW_ACCEPTANCE_H
#define MVMC_BACKFLOW_ACCEPTANCE_H
#include <math.h>

/* Returns validity separately from the physical accept/reject decision.
 * A checked zero candidate is represented by newLog == -INFINITY. */
static inline int BFLogAcceptance(double projection, double newLog, double oldLog,
    double u, int *accept) {
  *accept=0;
  if(!isfinite(projection) || (!isfinite(newLog) && newLog != -INFINITY)
      || !isfinite(oldLog) || !isfinite(u) || u < 0.0 || u >= 1.0) return 0;
  if(newLog == -INFINITY) return 1;
  const long double ratio=2.0L*((long double)projection+(long double)newLog-(long double)oldLog);
  *accept=logl((long double)u) < fminl(0.0L,ratio);
  return 1;
}
#endif
