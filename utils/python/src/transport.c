#include "evolve.h"
void vel_to_temp(void) {
    int i,j,k;
    i = j = k =0;
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                vx_temp[l] = vx[l];
                vy_temp[l] = vy[l];
            }
        }
    }
    return;
}
void temp_to_vel(void) {
    int i,j,k;
    i = j = k = 0;
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                vx[l] = vx_temp[l];
                vy[l] = vy_temp[l];
            }
        }
    }
    return;
}
void transport_step(void) {
    set_momenta();    
    transportY();
    transportX();
    set_vel();
    return;
}

void set_momenta(void) {
    int i,j,k;
    i=j=k=0;
    for(j=0;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {
            Pixm[l] = ymed(j)*(vx_temp[l]+omf*ymed(j))*dens[l];
            Pixp[l] = ymed(j)*(vx_temp[lxp]+omf*ymed(j))*dens[l];
            Piym[l] = vy_temp[l]*dens[l];
            Piyp[l] = vy_temp[lyp]*dens[l];
        }
    }

    return;
}
void transportY(void) {
    // Y direction

    vanleer_y_a(dens);
    vanleer_y_b(dens,denstar,dt);

    DividebyRho(Pixm);
    vanleer_y_a(divrho);
    vanleer_y_b(divrho,Qs,dt);
    update_flux(Qs);
    updateY(Pixm,Qs,dt);

    DividebyRho(Pixp);
    vanleer_y_a(divrho);
    vanleer_y_b(divrho,Qs,dt);
    update_flux(Qs);
    updateY(Pixp,Qs,dt);

    DividebyRho(Piym);
    vanleer_y_a(divrho);
    vanleer_y_b(divrho,Qs,dt);
    updateY(Piym,Qs,dt);

    DividebyRho(Piyp);
    vanleer_y_a(divrho);
    vanleer_y_b(divrho,Qs,dt);
    updateY(Piyp,Qs,dt);

    update_density_Y(dt);


    vanleer_y_a_avg(dbar);
    vanleer_y_b_avg(dbar,dbarstar,dt);

    DividebyRhoavg(Ld);
    vanleer_y_a_avg(divrho);
    vanleer_y_b_avg(divrho,Qs,dt);
    update_flux_avg(Qs);



    return;

}
void DividebyRho(double *q) {
    int i,j,k;
    i=j=k=0;
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            divrho[l] = q[l]/dens[l];
        }
    }
    return;
}
void DividebyRhoavg(double *q) {
    int j;
    for(j=0;j<size_y;j++) {
        divrho[j] = q[j]/dbar[j];
    }
    return;
}
void transportX(void) {
    // X direction

    vanleer_x_a(dens);
    vanleer_x_b(dens,denstar,dt);

    DividebyRho(Pixm);
    vanleer_x_a(divrho);
    vanleer_x_b(divrho,Qs,dt);
    updateX(Pixm,Qs,dt);

    DividebyRho(Pixp);
    vanleer_x_a(divrho);
    vanleer_x_b(divrho,Qs,dt);
    updateX(Pixp,Qs,dt);

    DividebyRho(Piym);
    vanleer_x_a(divrho);
    vanleer_x_b(divrho,Qs,dt);
    updateX(Piym,Qs,dt);

    DividebyRho(Piyp);
    vanleer_x_a(divrho);
    vanleer_x_b(divrho,Qs,dt);
    updateX(Piyp,Qs,dt);

    update_density_X(dt);

    return;

}
void set_vel(void) {

    int i,j,k;
    i=j=k=0;
    for(j=1;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            vy[l] = (Piym[l] + Piyp[lym])/(dens[l]+dens[lym]);
            vx[l] = (Pixm[l] + Pixp[lxm])/(ymed(j)*(dens[l]+dens[lxm])) - omf*ymed(j);
        }
    }
/*
    for(j=1;j<size_y;j++) {
        resv = 0;
        for(i=0;i<size_x;i++) {
            resv += .5*(Pixm[l] + Pixp[l]);
            //resv += dens[l]*(.5*(vx[lxp]+vx[l])+omf*ymed(j))*ymed(j);
        }
        resv /= (double)nx;
        dtLt[j] = -resv;
        dtLd[j] = -dbar[j]*(vx[j] + omf*ymed(j))*ymed(j);
    }
*/
    return;
}
