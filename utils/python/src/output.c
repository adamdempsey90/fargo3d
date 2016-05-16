#include "evolve.h"
void output_init(char *directory) {
    FILE *f;
    char fname[256];
    sprintf(fname,"%soutput_init.dat",directory);
    f = fopen(fname,"w");
    fwrite(&dens[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vy[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vx[0],sizeof(double),size_x*size_y*size_z,f);
    fclose(f);
    return;
}
void output(char *directory) {
    FILE *f;
    char fname[256];
    sprintf(fname,"%soutput.dat",directory);
    f = fopen(fname,"w");
    fwrite(&Ymin[0],sizeof(double),size_y+1,f);
    fwrite(&Xmin[0],sizeof(double),size_x+1,f);
    fwrite(&dens[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vy[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vx[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&Pot[0],sizeof(double),size_x*size_y*size_z,f);

    fclose(f);
    return;
}
void output_torque(char *directory) {
    FILE *f;
    char fname[256];
    sprintf(fname,"%storque.dat",directory);
    f = fopen(fname,"w");
    fwrite(&Ymed[NGHY],sizeof(double),ny,f);
    fwrite(&dbart[NGHY],sizeof(double),ny,f);
    fwrite(&Lt[NGHY],sizeof(double),ny,f);
    fwrite(&Lt[NGHY+size_y],sizeof(double),ny,f);
    fwrite(&Ld[NGHY],sizeof(double),ny,f);
    fwrite(&Ld[NGHY+size_y],sizeof(double),ny,f);
    fwrite(&Lw[NGHY],sizeof(double),ny,f);
    fwrite(&Lw[NGHY+size_y],sizeof(double),ny,f);
    fwrite(&drFt[NGHY],sizeof(double),ny,f);
    fwrite(&drFd[NGHY],sizeof(double),ny,f);
    fwrite(&drFw[NGHY],sizeof(double),ny,f);
    fwrite(&Lamex[NGHY],sizeof(double),ny,f);
    fwrite(&Lamex[NGHY+size_y],sizeof(double),ny,f);
    fwrite(&Lamdep[NGHY],sizeof(double),ny,f);
    fwrite(&dtLt[NGHY],sizeof(double),ny,f);
    fwrite(&dtLd[NGHY],sizeof(double),ny,f);
    fwrite(&dtLw[NGHY],sizeof(double),ny,f);
    fwrite(&mdotavg[NGHY],sizeof(double),ny,f);

    fclose(f);
    return;
}
