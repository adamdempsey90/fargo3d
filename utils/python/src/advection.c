#include "evolve.h"

void advect_shift(double *q, double *nshift) {
    /* Assume NGHX = 0 */
    int i,j,k;
    i=j=k=0;
    
    int itarget,ltarget;

    for(k=0;k<size_x;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
            

	            itarget = i-nshift[((j)+((k)*(ny+2*NGHY)))]; // l2D_int
	            while (itarget <  0)  {
                    itarget += nx;
                }
	            while (itarget >= nx+0) {
                    itarget -= nx;
                }
	            ltarget = l-i+itarget;
	            Pres[l] = q[ltarget];
            }
        }
    }

    memcpy(q,Pres,sizeof(double)*size_x*size_y*size_z);
    return;

}
void compute_residual(double dt) {
    int i,j,k;
    i=j=k=0;
    double ntilde, nround;
    double res;

    for(k=0;k<size_x;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                res += vx_temp[l];
            }
            res /=(double)nx;
            for(i=0;i<size_x;i++) {
                ntilde = res*dt/zone_size_x(j,k);
                nround = floor(ntilde+0.5);
                if(i == 0) {
                    nshift[l2D_int] = (int)nround;
                }
                vx[l] = vx_temp[l] - res;
                vx_temp[l] = (ntilde-nround)*zone_size_x(j,k)/dt;
            }
        }
    }




    return;
}
