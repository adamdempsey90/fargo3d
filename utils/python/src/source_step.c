#include "evolve.h"
void source_step(void) {
    int i,j,k;
    i=j=k=0;
    double vxc;
    compute_Pres();
    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {
            // X
            vx_temp[l] = vx[l] -2*dt/(dens[l]+dens[lxm]) *(Pres[l]-Pres[lxm])/zone_size_x(j,k);
            vx_temp[l] -= dt*(Pot[l]-Pot[lxm])/zone_size_x(j,k);
            if (IndirectTerm) {
                vx_temp[l] -= dt*(indPot[l]-indPot[lxm])/zone_size_x(j,k);
            }
            // Y
            vy_temp[l] = vy[l] -2*dt/(dens[l]+dens[lym])*(Pres[l]-Pres[lym])/(ymed(j)-ymed(j-1));
            vxc = .25*(vx[l]+vx[lxp]+vx[lym]+vx[lxp-pitch])+omf*ymin(j);
            vy_temp[l] += dt*vxc*vxc/ymin(j);
            vy_temp[l] -= dt*(Pot[l]-Pot[lym])/(ymed(j)-ymed(j-1));
            if (IndirectTerm) {
                vy_temp[l] -= dt*(indPot[l]-indPot[lym])/(ymed(j)-ymed(j-1));
            }
        }
    }


    return;
}
