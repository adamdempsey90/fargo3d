#include "evolve.h"
double cfl(void) {
    int i,j,k;
    i=j=k=0;
    double soundspeed2,soundspeed,visc;
    double cfl1_a, cfl1_b, cfl1;
    double cfl7_a, cfl7_b, cfl7;
    double cfl5_a,cfl5_b,cfl5;
    double cfl2, cfl3;
    double res,fac,vxx,vxxp;
#ifdef FARGO
    double med;
#endif
    cfl5_a = 0;
    cfl5_b = 0;
    res = 1e99;
    for(j=NGHY;j<size_y-NGHY;j++) {
#ifdef FARGO
        med = 0;
        for(i=0;i<size_x;i++) {

            med += vx[l];
        }
        med /= (double)nx;
#endif
        for(i=0;i<size_x;i++) {
#ifdef FARGO
            vxx = vx[l] - med;
            vxxp = vx[lxp] - med;
#else
            vxx = vx[l];
            vxxp = vx[lxp];
#endif
            soundspeed2 = energy[l]*energy[l];
            visc = params.alpha*energy[l]*energy[l]*sqrt(ymed(j)*ymed(j)*ymed(j));
            soundspeed = sqrt(soundspeed2);
            cfl1_a = soundspeed/zone_size_x(j,k);
            cfl1_b = soundspeed/zone_size_y(j,k);
            cfl1 = fmax(cfl1_a,cfl1_b);

	        cfl2 =fmax(fabs(vxx),fabs(vxxp))/zone_size_x(j,k);
	        cfl3 = fmax(fabs(vy[l]),fabs(vy[lyp]))/zone_size_y(j,k);

	        cfl7_a = 1.0/zone_size_x(j,k);	
        	cfl7_b = 1.0/zone_size_y(j,k);
	        cfl7 = 4.0*visc*pow(fmax(cfl7_a,cfl7_b),2);
#ifdef ARTIFICIALVISCOSITY
            cfl5_a = fabs(vxxp-vxx)/zone_size_x(j,k);
            cfl5_b = fabs(vy[lyp]-vy[l])/zone_size_y(j,k);
            cfl5 = fmax(cfl5_a, cfl5_b)*4.0*CVNR;
#else
            cfl5 = 0;
#endif

	        fac = CFL/sqrt(cfl1*cfl1 + cfl2*cfl2 + 
			     cfl3*cfl3 + cfl5*cfl5 +  cfl7*cfl7);

            if (fac < res) {
                res = fac;
//                printf("%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",cfl1,cfl2,cfl3,cfl5,cfl7,fac);
            }
        }
    }
    

    return res;

}
