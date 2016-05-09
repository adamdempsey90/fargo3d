#include "evolve.h"
void compute_Pres(void) {
    int i,j,k;
    i=j=k=0;
    double cs2;
    for(j=0;j<size_y;j++) {
        cs2 = Cs(ymed(j));
        cs2 *= cs2;
        for(i=0;i<size_x;i++) {
            Pres[l] = cs2*dens[l];
        }
    }
    return;
}
void compute_energy(void) {
    int i,j,k;
    i=j=k=0;
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            energy[l] = params.h * pow(ymed(j),params.flaringindex) * sqrt(1./ymed(j));
        }
    }
    return;
}
