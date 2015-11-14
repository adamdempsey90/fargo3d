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

void boundary_ymax_cpu () {

//<USER_DEFINED>
INPUT(Density);
INPUT(Vx);
INPUT(Vy);
OUTPUT(Density);
OUTPUT(Vx);
OUTPUT(Vy);
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
  int size_x = Nx+2*NGHX;
  int size_y = NGHY;
  int size_z = Nz+2*NGHZ;
  int nx = Nx;
  int ny = Ny;
  int nz = Nz;
  int nghy = NGHY;
  int nghz = NGHZ;
  int pitch  = Pitch_cpu;
  int stride = Stride_cpu;
  real dx = Dx;
  real flaringindex = FLARINGINDEX;
  real omegaframe = OMEGAFRAME;
  real alphaviscosity = ALPHAVISCOSITY;
  real aspectratio = ASPECTRATIO;
  real sigmaout = SIGMAOUT;
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

	lgh = i + (ny+nghy+j)*pitch + k*stride;
	lghs = i + (ny+nghy+1+j)*pitch + k*stride;
	lact = i + (ny+nghy-1-j)*pitch + k*stride;
	lacts = i + (ny+nghy-1-j)*pitch + k*stride;
	lacts_null = i + (ny+nghy)*pitch + k*stride;
	jgh = (ny+nghy+j);
	jact = (ny+nghy-1-j);

	density[lgh] = sigmaout;
	vx[lgh] = (vx[lact]+ymed(jact)*omegaframe)*sqrt(ymed(jact)/ymed(jgh))-ymed(jgh)*omegaframe;
	if (j<size_y-1)
		vy[lghs] = -3*alphaviscosity*aspectratio*aspectratio*pow(ymin(jgh),-0.5+2*flaringindex)*(1+2*flaringindex + log(density[lgh]/density[lact])/(2*log(ymed(jgh)/ymed(jact))));
	vy[lacts_null] = -3*alphaviscosity*aspectratio*aspectratio*pow(ymin(jgh),-0.5+2*flaringindex)*(1+2*flaringindex + log(density[lgh]/density[lact])/(2*log(ymed(jgh)/ymed(jact))));
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
