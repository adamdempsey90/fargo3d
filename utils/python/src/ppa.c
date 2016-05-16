#include "evolve.h"


void vanleer_ppa_a(double *q) {
    int i,j,k;
    i=j=k=0;

    double diff;
    double cord;
    double dqm, dqp,temp;


    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {

                dqm = (q[l]-q[lxm]);
                dqp = (q[lxp]-q[l]);
                if(dqp*dqm<=0.0)  {
                    slope[l] = 0.0;
                }
                else { // Monotonized centered slope limited
                    slope[l] = 0.5*(q[lxp]-q[lxm]);
                    temp = fabs(slope[l]);
                    if (2.0*fabs(dqm) < temp) {
                        temp = 2.0*fabs(dqm);
                    }
                    if (2.0*fabs(dqp) < temp) {
                        temp = 2.0*fabs(dqp);
                    }
                    if (slope[l] < 0) {
                        slope[l] = -temp;
                    }
                    else {
                        slope[l] = temp;
                    }
                }
            }
        }
    }


    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                temp = q[l]+0.5*(q[lxp]-q[l])-1.0/6.0*(slope[lxp]-slope[l]);
                qR[l] = temp;
                qL[lxp] = temp;
            }
        }
    }
    


    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                if ((qR[l]-q[l])*(q[l]-qL[l]) < 0.0) {
                  qL[l] = q[l];
                  qR[l] = q[l];
                }
                diff = qR[l] - qL[l];
                cord = q[l] - 0.5*(qL[l]+qR[l]);
                if (6.0*diff*cord > diff*diff) { 
                    qL[l] = 3.0*q[l]-2.0*qR[l];
                }
                if (-diff*diff > 6.0*diff*cord) { 
                    qR[l] = 3.0*q[l]-2.0*qL[l];
                }
            }
        }
    }

    return;
}
void vanleer_ppa_b(double dt, double *q, double *qs, double *vxt) {
    int i,j,k;
    i=j=k=0;
    double ksi;


    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                if (vxt[l] > 0.0) {
                  ksi = vxt[l]*dt/zone_size_x(j,k);
                  qs[l] = qR[lxm]+ksi*(q[lxm]-qR[lxm]);
                  qs[l]+= ksi*(1.0-ksi)*(2.0*q[lxm]-qR[lxm]-qL[lxm]);
                } else {
                  ksi = -vxt[l]*dt/zone_size_x(j,k);
                  qs[l] = qL[l]+ksi*(q[l]-qL[l]);
                  qs[l]+= ksi*(1.0-ksi)*(2.0*q[l]-qR[l]-qL[l]);
                }
            }
        }
    }


    return;
}
