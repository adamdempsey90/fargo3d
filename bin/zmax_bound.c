//<FLAGS>
//#define __GPU
//#define __NOPROTO
//<\FLAGS>

//<INCLUDES>
#include "fargo3d.h"
#define LEFT  0
#define RIGHT 1
#define DOWN  2
#define UP    3
//<\INCLUDES>

void boundary_zmax_cpu () {

//<USER_DEFINED>
INPUT(Density);
INPUT(Vx);
INPUT(Vy);
INPUT(Vz);
OUTPUT(Density);
OUTPUT(Vx);
OUTPUT(Vy);
OUTPUT(Vz);
//<\USER_DEFINED>

//<INTERNAL>
  int i;
  int j;
  int k;
  int jact;
  int jgh;
  int kact;
  int kgh;
  int lgh;
  int lghs;
  int lact;
  int lacts;
  int lacts_null;
//<\INTERNAL>

//<EXTERNAL>
  real* density = Density->field_cpu;
  real* vx = Vx->field_cpu;
  real* vy = Vy->field_cpu;
  real* vz = Vz->field_cpu;
  int size_x = Nx+2*NGHX;
  int size_y = Ny+2*NGHY;
  int size_z = NGHZ;
  int nx = Nx;
  int ny = Ny;
  int nz = Nz;
  int nghy = NGHY;
  int nghz = NGHZ;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  real dx = Dx;
//<\EXTERNAL>

//<CONSTANT>
// real xmin(Nx+2*NGHX+1);
// real ymin(Ny+2*NGHY+1);
// real zmin(Nz+2*NGHZ+1);
//<\CONSTANT>

//<MAIN_LOOP>

  i = j = k = 0;

#ifdef Z
  for(k=0; k<size_z; k++) {
#endif
#ifdef Y
    for(j=0; j<size_y; j++) {
#endif
#ifdef X
      for(i=0; i<size_x; i++) {
#endif
//<#>

	lgh = i + j*pitch + (nz+nghz+k)*stride;
	lghs = i + j*pitch + (nz+nghz+1+k)*stride;
	lact = i + j*pitch + (nz+nghz-1-k)*stride;
	lacts = i + j*pitch + (nz+nghz-1-k)*stride;
	lacts_null = i + j*pitch + (nz+nghz)*stride;
	kgh = (nz+nghz+k);
	kact = (nz+nghz-1-k);

	density[lgh] = density[lact];
	vx[lgh] = vx[lact];
	vy[lgh] = vy[lact];
	if (k<size_z-1)
		vz[lghs] = -vz[lacts];
	vz[lacts_null] = 0;
//<\#>
#ifdef X
      }
#endif
#ifdef Y
    }
#endif
#ifdef Z
  }
#endif
//<\MAIN_LOOP>
}
