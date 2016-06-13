#include "evolve.h"
void updateX(double *q, double *qs,double dt,double *vxt) {
    int i,j,k;
    i=j=k=0;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {

                q[l] += ((vxt[l]*qs[l]*denstar[l]-vxt[lxp]*qs[lxp]*denstar[lxp])*SurfX(j,k)*dt*InvVol(j,k));

            }
        }
    }

    return;
}
void update_flux_avg(double *qs) {
    int i,j,k;
    int indx;
    i=j=k=0;
    double res;
#ifdef _OPENMP
    #pragma omp parallel for private(i,j,k,res)
#endif
    for(j=0;j<size_y-1;j++) {
        i = 0;
        convolution(&vy_temp[l],&denstar[l],drFd,-dt*qs[j]*SurfY(j,k)*InvVol(j,k),j,size_y);
        convolution(&vy_temp[lyp],&denstar[lyp],drFd,dt*qs[j+1]*SurfY(j+1,k)*InvVol(j,k),j,size_y);
    
        
        res = 0;
        for(i=0;i<size_x;i++) {    
        
            res += ((vy_temp[l]*qs[j]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*qs[j+1]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));
        }
        res /=(double)nx;

        drFd[j] -= dt*res;
    }
    return;
}
void update_flux(double *qs) {
    int i,j,k;
    i=j=k=0;
    double res;
#ifdef _OPENMP
    #pragma omp parallel for private(i,j,res)
#endif
    for(j=0;j<size_y-1;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {    
            res += ((vy_temp[l]*qs[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*qs[lyp]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));
        }
        res /=(double)nx;
        drFt[j] -= dt*res*.5;
    }
    return;
}
void updateY(double *q, double *qs,double dt) {
    int i,j,k;
    i=j=k=0;
#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y-1;j++) {
            for(i=0;i<size_x;i++) {
                
                q[l] += dt*((vy_temp[l]*qs[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*qs[lyp]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));
            }
        }
    }
    return;
}
void update_density_Y(double dt) {
    int i,j,k;
    i=j=k=0;
    double fac,res;
#ifdef _OPENMP
    #pragma omp parallel for collapse(2) private(i,j,k,res,fac)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y-1;j++) {
            res = 0;
            fac = 0;
            for(i=0;i<size_x;i++) {

                fac = ((vy_temp[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));
                dens[l] += fac*dt;
                res += fac;

            }
            res /= (double)nx;
            mdotavg[j] += res;
        }
    }
    return;
}
void update_density_X(double dt,double *vxt) {
    int i,j,k;
    i = j= k =0;

#ifdef _OPENMP
    #pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {

                dens[l] += ((vxt[l]*denstar[l]-vxt[lxp]*denstar[lxp])*SurfX(j,k)*dt*InvVol(j,k));

            }
        }
    }
    return;
}
