#include "evolve.h"
void set_avg(int p) {
    int i,j,k;
    i=j=k=0;
    double resx,resy,resd;
    double resL;
    for(j=0;j<size_y;j++) {
        resx = 0;
        resy = 0;
        resd = 0;
        resL = 0;

        convolution(dens,vx,Lt,ymed(j),j,j,k);
        for(i=0;i<size_x;i++) {
            resx += vx[l];
            resy += vy[l];
            resd += dens[l];
            resL += ymed(j)*dens[l]*(.5*(vx[l] + vx[lxp]) + omf*ymed(j));
        }
        vxbar[j] = resx/(double)nx;
        vybar[j] = resy/(double)nx;
        dbar[j] = resd/(double)nx;
        Lt[j + p*size_y] = resL/(double)nx;
        if (p==1) dbart[j] += dbar[j]*dt;

        Ld[j + p*size_y] = ymed(j)*( vxbar[j] + omf*ymed(j))*dbar[j];
        Lw[j + p*size_y] = Lt[j + p*size_y] - Ld[j + p*size_y];
    }
    return;
}

void calc_modes(double *q, double *qmr, double *qmi, int m_max) {

    int i,j,k;
    int mi,lm; 
    i = j = k = 0;


#ifdef _OPENMP
#pragma omp parallel for collapse(2) private(lm,i,j,mi)
#endif
    for(j=0;j<size_y;j++) {
        for(mi=0;mi<m_max;mi++) {
            lm = mi + j*m_max; 
            for(i=0;i<size_x;i++) {
                qmr[lm] += q[l] * cos(mi*xmed(i));
                qmi[lm] += q[l] * sin(mi*xmed(i));
            }
            qmr[lm] /= (double)nx;
            qmi[lm] /= (double)nx;
        }
    }

    return;
}
