#include "evolve.h"
void time_avg(void) {
    int j;

    for(j=0;j<size_y;j++) {
        dbart[j]/=(double)nsteps;
        //Lt[j]/=(double)nsteps;
        //Lt[j+size_y]/=(double)nsteps;
        //Ld[j]/=(double)nsteps;
        //Ld[j+size_y]/=(double)nsteps;
        Lw[j]/=(double)nsteps;
        Lw[j+size_y]/=(double)nsteps;
        drFt[j]/=(double)nsteps;
        drFd[j]/=(double)nsteps;
        drFw[j]/=(double)nsteps;
        Lamex[j]/=(double)nsteps;
        Lamex[j+size_y]/=(double)nsteps;
        Lamdep[j]/=(double)nsteps;
        dtLt[j]/=(double)nsteps;
        dtLd[j]/=(double)nsteps;
        dtLw[j]/=(double)nsteps;


    }
    return;
}
