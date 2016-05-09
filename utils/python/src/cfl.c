#include "evolve.h"
double cfl(void) {
    int i,j,k;
    i=j=k=0;
    double soundspeed2,soundspeed,visc;
    double cfl1_a, cfl1_b, cfl1;
    double cfl7_a, cfl7_b, cfl7;
    double cfl2, cfl3;
    double res,fac;

    res = 1e99;
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {
            soundspeed2 = energy[l]*energy[l];
            visc = params.alpha*energy[l]*energy[l]*sqrt(ymed(j)*ymed(j)*ymed(j));
            soundspeed = sqrt(soundspeed2);
            cfl1_a = soundspeed/zone_size_x(j,k);
            cfl1_b = soundspeed/zone_size_y(j,k);
            cfl1 = fmax(cfl1_a,cfl1_b);

	        cfl2 =fmax(fabs(vx[l]),fabs(vx[lxp]))/zone_size_x(j,k);
	        cfl3 = fmax(fabs(vy[l]),fabs(vy[lyp]))/zone_size_y(j,k);

	        cfl7_a = 1.0/zone_size_x(j,k);	
        	cfl7_b = 1.0/zone_size_y(j,k);
	        cfl7 = 4.0*visc*pow(fmax(cfl7_a,cfl7_b),2);


	        fac = CFL/sqrt(cfl1*cfl1 + cfl2*cfl2 + 
			     cfl3*cfl3 + cfl7*cfl7);

            if (fac < res) {
                res = fac;
            }
        }
    }

    return res;

}
