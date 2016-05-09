#include "evolve.h"
void get_accel(double *q, double *k, double dtn) {
    double rp,coeff;
    rp = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
    coeff = -G*MSTAR*pow(rp,-3.0);
    k[0] = q[3] * dtn;
    k[1] = q[4] * dtn;
    k[2] = q[5] * dtn;
    k[3] = coeff * q[0] * dtn;
    k[4] = coeff * q[1] * dtn;
    k[5] = coeff * q[2] * dtn;
    coeff = G*planet.mp/(rp*rp*rp);
    k[3] -= coeff * q[0] * dtn;
    k[4] -= coeff * q[1] * dtn;
    k[5] -= coeff * q[2] * dtn;
    return;
}
void move_planet_step(double dtn) {
    int i;
    double k1[6] ={0,0,0,0,0,0};
    double k2[6] ={0,0,0,0,0,0};
    double k3[6] ={0,0,0,0,0,0};
    double k4[6] ={0,0,0,0,0,0};
    double k5[6] ={0,0,0,0,0,0};
    double k6[6] ={0,0,0,0,0,0};
    double q0[6] ={0,0,0,0,0,0};
    double q1[6] ={0,0,0,0,0,0};
    double cvals1[5] =  {0.2, 0.0, 0.0, 0.0, 0.0};
    double cvals2[5] = {0.075, 0.225, 0.0, 0.0, 0.0};
    double cvals3[5] = { 0.3, -0.9, 1.2, 0.0, 0.0};
    double cvals4[5] = {-11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0, 0.0 };
    double cvals5[5] = {1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0 };
    q0[0] = planet.x;
    q0[1] = planet.y;
    q0[2] = planet.z;
    q0[3] = planet.vx;
    q0[4] = planet.vy;
    q0[5] = planet.vz;

    get_accel(q0, k1, dtn);
    for(i=0;i<6;i++) {
        q1[i] = q0[i] + cvals1[0]*k1[i] + cvals1[1]*k2[i] + cvals1[2]*k3[i] + cvals1[3]*k4[i] + cvals1[4]*k5[i];
    }

    get_accel(q1, k2, dtn);
    for(i=0;i<6;i++) {
        q1[i] = q0[i] + cvals2[0]*k1[i] + cvals2[1]*k2[i] + cvals2[2]*k3[i] + cvals2[3]*k4[i] + cvals2[4]*k5[i];
    }


    get_accel(q1, k3, dtn);
    for(i=0;i<6;i++) {
        q1[i] = q0[i] + cvals3[0]*k1[i] + cvals3[1]*k2[i] + cvals3[2]*k3[i] + cvals3[3]*k4[i] + cvals3[4]*k5[i];
    }

    get_accel(q1, k4, dtn);
    for(i=0;i<6;i++) {
        q1[i] = q0[i] + cvals4[0]*k1[i] + cvals4[1]*k2[i] + cvals4[2]*k3[i] + cvals4[3]*k4[i] + cvals4[4]*k5[i];
    }

    get_accel(q1, k5, dtn);
    for(i=0;i<6;i++) {
        q1[i] = q0[i] + cvals5[0]*k1[i] + cvals5[1]*k2[i] + cvals5[2]*k3[i] + cvals5[3]*k4[i] + cvals5[4]*k5[i];
    }

    get_accel(q1, k6, dtn);

    for (i = 0; i < 6; i++) {
        q1[i] = (q0[i]+
                    37.0/378.0*k1[i]  + 
                    250.0/621.0*k3[i] + 
                    125.0/594.0*k4[i] + 
                    512.0/1771.0*k6[i]);
    }
    planet.x = q1[0];
    planet.y = q1[1];
    planet.z = q1[2];
    planet.vx = q1[3];
    planet.vy = q1[4];
    planet.vz = q1[5];
    return;
}
void rotate_sys(double angle) {

    double xt, yt,vxt,vyt;
    xt = planet.x;
    yt = planet.y;
    vxt = planet.vx;
    vyt = planet.vy;

    planet.x = xt * cos(angle)  + yt * sin(angle);
    planet.y = -xt * sin(angle)  + yt * cos(angle);

    planet.vx = vxt * cos(angle)  + vyt * sin(angle);
    planet.vy = -vxt * sin(angle)  + vyt * cos(angle);

    return;
}
void move_planet(void) {
    int i,j,k;
    int subcycles = 5;
    double dt_frac = 1./(double)subcycles;
    double omf_new;
    double xpl0,ypl0;
    double xpl1,ypl1;
    double d1, d2,cross,domega;

    xpl0 = planet.x;
    ypl0 = planet.y;

    for(i=0;i<subcycles;i++) {
        move_planet_step(dt*dt_frac);
    }
    xpl1 = planet.x;
    ypl1 = planet.y;
    
    d2 = sqrt(xpl1*xpl1 + ypl1*ypl1);
    d1 = sqrt(xpl0*xpl0+ypl0*ypl0);
    cross = xpl0*ypl1-xpl1*ypl0;
    omf_new = asin(cross/(d1*d2))/dt;
    domega = omf_new - omf;
    omf = omf_new;
    rotate_sys(omf*dt);

    i = j = k = 0;
    double res=0;
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            res = 0;
            for(i=0;i<size_x;i++) {
                vx[l] -= domega * ymed(j);
                res += vx[l];
            }
            res /=(double)nx;
            vxbar[j] = res;
        }
    }


    return;
}
