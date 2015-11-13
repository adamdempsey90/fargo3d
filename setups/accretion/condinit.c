#include "fargo3d.h"

void InitDensPlanet() {
 printf("%lg\t%lg\n",MDOT0,MDOT);
  int i,j,k;
  real *field;
  real viscosity,fac1,fac2,fac3,nuind,nu0,r;
  field = Density->field_cpu;
  nuind = 0.5+2*FLARINGINDEX;
  nu0 = ASPECTRATIO*ASPECTRATIO*ALPHAVISCOSITY;
  boolean GhostInclude = TRUE;
  
  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = (GhostInclude ? 0 : NGHX);
  int end_i = Nx+2*NGHX-begin_i;
   for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
        r = Ymed(j);
        viscosity = nu0*pow(r,nuind);
        fac1=pow(YMAX,.5-nuind)*YMIN*SIGMAIN-pow(YMIN,.5-nuind)*YMAX*SIGMAOUT;
        fac2=-pow(YMAX,1-nuind)*YMIN*SIGMAIN+pow(YMIN,1-nuind)*YMAX*SIGMAOUT;
        fac3 = pow(YMIN,1-nuind)*pow(YMAX,0.5-nuind)-pow(YMAX,1-nuind)*pow(YMIN,0.5-nuind);
        fac1 *= pow(r,-nuind);
        fac2 *= pow(r,-.5*-nuind);
      for (i = begin_i; i<end_i; i++) {
          //
          //
//         field[l] = SIGFLOOR + SIGMA0 * exp( - (Ymed(j)-1.8)*(Ymed(j)-1.8)/.01);
        field[l]= (fac1-fac2)/fac3;
      }
    }
  }
}

void InitSoundSpeedPlanet() {

  int i,j,k;
  real *field;
  real dr, dz;
  real r, z, H, r0, rho_o, t, omega, vk;
  real rho;
  FILE *fo;
  real *d;
  real *e;

  field = Energy->field_cpu;
  d = Density->field_cpu;

  boolean GhostInclude = TRUE;
  
  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = (GhostInclude ? 0 : NGHX);
  int end_i = Nx+2*NGHX-begin_i;
  
  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
      for (i = begin_i; i<end_i; i++) {	
	r = Ymed(j);
	vk = sqrt(G*MSTAR/r);
	field[l] = ASPECTRATIO * pow(Ymed(j)/R0, FLARINGINDEX) * vk; //sqrt(G*MSTAR/Ymed(j))
#ifdef ADIABATIC
	field[l] = field[l]*field[l]*d[l]/(GAMMA-1.0);
#endif
      }
    }
  }    
}

void InitVazimPlanet() {

  int i,j,k;
  real *field;
  real dr, dz;
  real r, z, H, r0, rho_o, t;
  real rho;
  FILE *fo;
  real vt, omega;
  real *vr;
  real *cs;
  real mdot_init,nuind,nu0; 
  real vr0, fac1, fac2;
  field = Vx->field_cpu;
  vr = Vy->field_cpu;
  cs = Energy->field_cpu;
    
  boolean GhostInclude = TRUE;
  nuind = 0.5+2*FLARINGINDEX;
  nu0 = ALPHAVISCOSITY*ASPECTRATIO*ASPECTRATIO;
  


  int begin_k =(GhostInclude ? 0 : NGHZ);
  int end_k = Nz+2*NGHZ-begin_k;
  int begin_j =(GhostInclude ? 0 : NGHY);
  int end_j = Ny+2*NGHY-begin_j;
  int begin_i = (GhostInclude ? 0 : NGHX);
  int end_i = Nx+2*NGHX-begin_i;

  for (k = begin_k; k<end_k; k++) {
    for (j = begin_j; j<end_j; j++) {
       r = Ymed(j); 
        vr0=-1.5*nu0*pow(r,nuind);
        fac1 = pow(YMAX,.5-nuind)*SIGMAIN*YMIN-pow(YMIN,.5-nuind)*SIGMAOUT*YMAX;
        fac2 = fac1 + (-pow(YMAX,1-nuind)*SIGMAIN*YMIN+pow(YMIN,1-nuind)*SIGMAOUT*YMAX)/sqrt(r);

      for (i = begin_i; i<end_i; i++) {
	
	omega = sqrt(G*MSTAR/r/r/r);
	vt = omega*r*sqrt(1.0+pow(ASPECTRATIO,2)*pow(r/R0,2*FLARINGINDEX)*
			  (2.0*FLARINGINDEX - 1.0 - SIGMASLOPE));

	vt -= OMEGAFRAME*r;

	field[l] = vt*(1.+ASPECTRATIO*NOISE*(drand48()-.5));
	//vr[l] = r*omega*ASPECTRATIO*NOISE*(drand48()-.5);
//    vr[l] = -1.5*ALPHAVISCOSITY*ASPECTRATIO*ASPECTRATIO*pow(r,-0.5+2*FLARINGINDEX);
    vr[l] = vr0*fac1/fac2;
      }
    }
  }    
}

void CondInit() {
  
  OUTPUT(Density);
  OUTPUT(Energy);
  OUTPUT(Vx);
  OUTPUT(Vy);

  int i,j,k;
  int index;
  real vt;
  real *field;
  real *rho;
  real *v1;
  real *v2;
  real *e;
#ifdef PLANETS
  Sys = InitPlanetarySystem(PLANETCONFIG);
  ListPlanets();
  if(COROTATING)
    OMEGAFRAME = GetPsysInfo(FREQUENCY);
  else
#endif
    OMEGAFRAME = OMEGAFRAME;

  InitDensPlanet ();
  InitSoundSpeedPlanet ();
  InitVazimPlanet ();
}
