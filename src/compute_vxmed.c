#include "fargo3d.h"

void ComputeVymed(Field *V) {
  int j,k;
  int ll2D;
  

  reduction_SUM(V, 0, Ny+2*NGHY, 0, Nz+2*NGHZ);

  INPUT2D(Reduction2D);
  OUTPUT2D(Vx0_avg);

  k = j = 0;

#ifdef Z
  for (k = 0; k < Nz+2*NGHZ; k++) {
#endif
#ifdef Y
    for (j = 0; j < Ny+2*NGHY; j++) {
#endif
      ll2D = l2D;
      Vx0_avg->field_cpu[ll2D] = Reduction2D->field_cpu[ll2D]/(real)Nx;
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
}
