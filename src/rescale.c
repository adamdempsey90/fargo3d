
#include "fargo3d.h"

void rescale () {
DT *= sqrt(R0*R0*R0/G/MSTAR);
OMEGAFRAME *= sqrt(G*MSTAR/(R0*R0*R0));
SIGMA0 *= MSTAR/(R0*R0);
}
