
#include "fargo3d.h"

void rescale () {
OMEGAFRAME *= sqrt(G*MSTAR/(R0*R0*R0));
DT *= sqrt(R0*R0*R0/G/MSTAR);
NU *= sqrt(G*MSTAR*R0);
SIGMA0 *= MSTAR/(R0*R0);
}
