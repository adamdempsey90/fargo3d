#include "evolve.h"
void updateX(double *q, double *qs,double dt) {
    int i,j,k;
    i=j=k=0;
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {

            q[l] += ((vx_temp[l]*qs[l]*denstar[l]-vx_temp[lxp]*qs[lxp]*denstar[lxp])*SurfX(j,k)*dt*InvVol(j,k));

        }
    }

    return;
}
void update_flux_avg(double *qs) {
    int i,j,k;
    i=j=k=0;
    double res;
    for(j=0;j<size_y-1;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {    
            res += ((vy_temp[l]*qs[j]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*qs[j+1]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));
        }
        res /=(double)nx;

        drFd[j] -= res;
    }
    return;
}
void update_flux(double *qs) {
    int i,j,k;
    i=j=k=0;
    double res;
    for(j=0;j<size_y-1;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {    
            res += ((vy_temp[l]*qs[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*qs[lyp]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));
        }
        res /=(double)nx;
        drFt[j] -= res*.5;
    }
    return;
}
void updateY(double *q, double *qs,double dt) {
    int i,j,k;
    i=j=k=0;
    for(j=0;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {
            
            q[l] += dt*((vy_temp[l]*qs[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*qs[lyp]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));
        }
    }
    return;
}
void update_density_Y(double dt) {
    int i,j,k;
    i=j=k=0;

    for(j=0;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {

            dens[l] += ((vy_temp[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*denstar[lyp]*SurfY(j+1,k))*dt*InvVol(j,k));

        }
    }
    return;
}
void update_density_X(double dt) {
    int i,j,k;
    i = j= k =0;

    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {

            dens[l] += ((vx_temp[l]*denstar[l]-vx_temp[lxp]*denstar[lxp])*SurfX(j,k)*dt*InvVol(j,k));

        }
    }
    return;
}
