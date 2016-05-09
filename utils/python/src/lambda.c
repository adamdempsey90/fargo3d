#include "evolve.h"
void set_Lamdep(void) {
    int j;
    for(j=NGHY;j<size_y - NGHY;j++) {
        dtLt[j] += (Lt[j+size_y]-Lt[j])/dt;
        dtLd[j] += (Ld[j+size_y]-Ld[j])/dt;
        dtLw[j] = dtLt[j] - dtLd[j];
        Lamdep[j] = dtLd[j] + drFd[j];
        drFw[j] = drFt[j] - drFd[j];

    }

    return;
}
void set_Lamex(void) {
    int i,j,k;
    i=j=k=0;
    
    double res,resi;
    for(j=NGHY;j<size_y-NGHY;j++) {
        res = 0;
        resi = 0;
        for(i=0;i<size_x;i++) {
            res -= dens[l]*(Pot[lxp]-Pot[lxm])/(2*dx);
            resi -= dens[l]*(indPot[lxp]-indPot[lxm])/(2*dx);
        }
        res /=(double)nx;
        resi /=(double)nx;
        Lamex[j] += res;
        Lamex[j + size_y] += resi;

    }

    return;

}
