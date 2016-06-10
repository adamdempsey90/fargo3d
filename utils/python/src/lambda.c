#include "evolve.h"
void set_Lamdep(void) {
    int j;
#ifdef _OPENMP
    #pragma omp parallel for private(j)
#endif
    for(j=NGHY;j<size_y - NGHY;j++) {
        dtLt[j] += (Lt[j+size_y]-Lt[j]);
        dtLd[j] += (Ld[j+size_y]-Ld[j]);
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
#ifdef _OPENMP
    #pragma omp parallel for private(i,j,res,resi)
#endif
    for(j=NGHY;j<size_y-NGHY;j++) {
        res = 0;
        resi = 0;
        
        convolution_deriv(dens,Pot,Lamex,-dt,j,j,k);

        for(i=0;i<size_x;i++) {
            res -= dens[l]*(Pot[lxp]-Pot[lxm])/(2*dx);
            resi -= dens[l]*(indPot[lxp]-indPot[lxm])/(2*dx);
        }
        res /=(double)nx;
        resi /=(double)nx;
        Lamex[j] += dt*res;
        Lamex[j + size_y] += dt*resi;

    }

    return;

}
