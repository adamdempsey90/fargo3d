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

void boundary_ymin_cpu () {

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
	lact = i + (2*nghy-j-1)*pitch + k*stride;
	lacts = i + (2*nghy-j)*pitch + k*stride;
	lacts_null = i + nghy*pitch + k*stride;
	jgh = j;
	jact = (2*nghy-j-1);

	density[lgh] = density[lact]*pow(ymed(jact)/ymed(jgh),sigmaslope+flaringindex+1.)*exp(-pow(cos(zmed(k)),2.)/(aspectratio*aspectratio)*(1.-ymed(jgh)/(.5*(ymed(jact)+ymed(jgh))))*flaringindex*pow(.5*(ymed(jgh)+ymed(jact))/r0,-2.*flaringindex-1.));
	vx[lgh] = (vx[lact]+omegaframe*ymed(jact)*sin(zmed(k)))*sqrt(ymed(jact)/ymed(jgh))*(1.+(2.+sigmaslope-flaringindex)*(ymed(jact)-ymed(jgh))/r0*flaringindex*aspectratio*aspectratio*pow((ymed(jgh)+ymed(jact))/(2.*r0),2.*flaringindex-1.))-ymed(jgh)*omegaframe*sin(zmed(k));
	vy[lghs] = -vy[lacts];
	vy[lacts_null] = 0;
	vz[lgh] = vz[lact];
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
