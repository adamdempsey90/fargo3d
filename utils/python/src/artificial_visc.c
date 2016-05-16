#include "evolve.h"


void artificial_visc(void) {
    int i,j,k;
    i=j=k=0;
    double dvx,dvy;
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y-1;j++) {
            for(i=0;i<size_x;i++) {

	            dvx = vx[lxp]-vx[l];
	            if (dvx < 0.0) {
	                Piym[l] = CVNR*CVNR*dens[l]*dvx*dvx;
	            }
	            else {
	                Piym[l] = 0.0;
	            }
                dvy = vy[lyp]-vy[l];
                if (dvy < 0.0) {
                  Piyp[l] = CVNR*CVNR*dens[l]*dvy*dvy;
                }
                else {
                  Piyp[l] = 0.0;
                }
            }
        }
    }


    for(k=0;k<size_z;k++) {
        for(j=1;j<size_y-1;j++) {
            for(i=0;i<size_x;i++) {

	            vx_temp[l] += - 2.0*(Piym[l]-Piym[lxm])/(dens[l]+dens[lxm])*dt/zone_size_x(j,k);
                vy_temp[l] += - 2.0*(Piyp[l]-Piyp[lym])/(dens[l]+dens[lym])*dt/zone_size_y(j,k);           
            }
        }
    }

    return;

}
