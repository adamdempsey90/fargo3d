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

void boundary_zmin_cpu () {

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
  real r0 = R0;
  real aspectratio = ASPECTRATIO;
  real flaringindex = FLARINGINDEX;
  real sigmaslope = SIGMASLOPE;
  real omegaframe = OMEGAFRAME;
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

	lgh = l;
	lghs = l;
	lact = i + j*pitch + (2*nghz-k-1)*stride;
	lacts = i + j*pitch + (2*nghz-k)*stride;
	lacts_null = i + j*pitch + nghz*stride;
	kgh = k;
	kact = (2*nghz-k-1);

	density[lgh] = density[lact]*pow(sin(zmed(kgh))/sin(zmed(kact)),flaringindex-2.-sigmaslope+1./(aspectratio*aspectratio)*pow(ymed(j)/r0,-2.*flaringindex));
	vx[lgh] = (vx[lact]+omegaframe*ymed(j)*sin(zmed(kact)))*(1.+flaringindex*cos(.5*(zmed(kgh)+zmed(kact)))*(zmed(kact)-zmed(kgh))) -ymed(j)*omegaframe*sin(zmed(kgh));
	vy[lgh] = vy[lact];
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
