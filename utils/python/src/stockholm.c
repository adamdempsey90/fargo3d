#include "evolve.h"


void stockholm(void) {
    int i,j,k;
    i=j=k=0;
    double tau,taud;
    double vy_target, vx_target,dens_target;
    double wkzin = params.wkzin;
    double wkzout = params.wkzout;
    double Y_inf = params.ymin + (params.ymax-params.ymin)*wkzin;
    double Y_sup = params.ymax - (params.ymax-params.ymin)*wkzout;


    double ds = 0.03333;
    double rampy = 0;
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                rampy = 0;

                if (ymed(j) > Y_sup) {
                    rampy = (ymed(j)-Y_sup)/(params.ymax-Y_sup);
                    vy_target = vy0[j];
                    vx_target = vx0[j];
                    dens_target= dens0[j];
                }
                if (ymed(j) < Y_inf) {
                    rampy = (Y_inf-ymed(j))/(Y_inf-params.ymin);
                    vy_target = vy0[j];
                    vx_target = vx0[j];
                    dens_target= dens0[j];
    
                }
                rampy *= rampy;
                tau = ds*pow(ymed(j),1.5);
                if (rampy > 0.0) {
                    taud = tau/rampy;
	                vx_target -= (omf-omf0)*ymed(j);
                    dens[l] = (dens[l]*taud + dens_target*dt)/(dt+taud);
                    vx[l] = (vx[l]*taud + vx_target*dt)/(dt+taud);
                    vy[l] = (vy[l]*taud + vy_target*dt)/(dt+taud);


                }

            }
        }
    }



}
void read_stockholm(char *directory) {
    char filename[512];
    FILE *fd,*fx,*fy;


    printf("Inside.\n");
    sprintf(filename,"%sdensity0_2d.dat",directory);
    fd = fopen(filename,"r");
    if (fd == NULL) printf("Error loading %s\n",filename); 

    sprintf(filename,"%svx0_2d.dat",directory);
    fx = fopen(filename,"r");
    if (fx == NULL) printf("Error loading %s\n",filename); 

    sprintf(filename,"%svy0_2d.dat",directory);
    fy = fopen(filename,"r");
    if (fy == NULL) printf("Error loading %s\n",filename); 

    fread(&dens0[0],sizeof(double),size_y*size_z,fd);
    fread(&vy0[0],sizeof(double),size_y*size_z,fy);
    fread(&vx0[0],sizeof(double),size_y*size_z,fx);
    fclose(fd);
    fclose(fx);
    fclose(fy);

    printf("%lg\t%lg\t%lg\t%lg\n",vy0[size_y-4],vy0[size_y-3],vy0[size_y-2],vy0[size_y-1]);
    
    return;

}
