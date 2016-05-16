#include "evolve.h"
void viscosity(void) {
    int i,j,k;
    i=j=k=0;
    double visc, viscm,div_v;
    double res,fac,facp;
    visc = 0;
    viscm = 0;
    compute_energy();
    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {

            visc = params.alpha*energy[l]*energy[l]*sqrt(ymed(j)*ymed(j)*ymed(j));
            viscm = params.alpha*.5*(energy[l]*energy[l]+energy[lym]*energy[lym])*sqrt(ymin(j)*ymin(j)*ymin(j));
            div_v = 0.0;
            div_v += (vx[lxp]-vx[l])*SurfX(j,k);
            div_v += (vy[lyp]*SurfY(j+1,k)-vy[l]*SurfY(j,k));
            div_v *= 2.0/3.0*InvVol(j,k);

            tauxx[l] = visc*dens[l]*(2.0*(vx[lxp]-vx[l])/zone_size_x(j,k) - div_v);
            tauxx[l] += visc*dens[l]*(vy[lyp]+vy[l])/ymed(j);
            tauyy[l] = visc*dens[l]*(2.0*(vy[lyp]-vy[l])/(ymin(j+1)-ymin(j)) - div_v);
            tauxy[l] = viscm*.25*(dens[l]+dens[lxm]+dens[lym]+dens[lxm-pitch])*((vy[l]-vy[lxm])/(dx*ymin(j)) + (vx[l]-vx[lym])/(ymed(j)-ymed(j-1))-.5*(vx[l]+vx[lym])/ymin(j)); //centered on left, inner vertical edge in z


        }
        tauxyavg[j] = viscm*.5*(dbar[j]+dbar[j-1])*( (vxbar[j]-vxbar[j-1])/(ymed(j)-ymed(j-1))-.5*(vxbar[j]+vxbar[j-1])/ymin(j)); //centered on left, inner vertical edge in z
    }

    for(j=1;j<size_y-2;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
            // X

            vx_temp[l] += 2.0*(tauxx[l]-tauxx[lxm])/(zone_size_x(j,k)*(dens[l]+dens[lxm]))*dt;
            fac =  (ymin(j+1)*ymin(j+1)*tauxy[lyp]-ymin(j)*ymin(j)*tauxy[l])/((ymin(j+1)-ymin(j))*ymed(j));
            facp =  (ymin(j+1)*ymin(j+1)*tauxy[lxp+pitch]-ymin(j)*ymin(j)*tauxy[lxp])/((ymin(j+1)-ymin(j))*ymed(j));
            vx_temp[l] += dt*fac*2.0/(ymed(j)*(dens[l]+dens[lxm]));
            // Y
            vy_temp[l] += 2.0*(ymed(j)*tauyy[l]-ymed(j-1)*tauyy[lym])/((ymed(j)-ymed(j-1))*(dens[l]+dens[lym])*ymin(j))*dt;
            vy_temp[l] += 2.0*(tauxy[lxp]-tauxy[l])/(dx*ymin(j)*(dens[l]+dens[lym]))*dt;
            vy_temp[l] -= (tauxx[l]+tauxx[lym])/(ymin(j)*(dens[l]+dens[lym]))*dt;
            res += .5*(fac+facp);
        }
         drFt[j] += -dt*res/(double)nx;
        //drFd[j] = -(ymin(j+1)*ymin(j+1)*tauxyavg[j+1]*SurfY(j+1,k) - ymin(j)*ymin(j)*tauxyavg[j]*SurfY(j,k))*InvVol(j,k);
        drFd[j]  +=  -dt*(ymin(j+1)*ymin(j+1)*tauxyavg[j+1]-ymin(j)*ymin(j)*tauxyavg[j])/((ymin(j+1)-ymin(j))*ymed(j));

    }
    return;
}
