#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif


#define DIMSLENGTH 14
#define DIMSNXINDEX 6
#define NGHY 3
#define l ( (i) + nx*(j) + nx*ny*(k) )
#define l2d ( (i) + nx*(j) )

double zmin = .2;
double h = .05;
double *zm, *ym, *ymed, *zmed;
int nx,ny,nz;


double scaleH(double r);
void create_grid_file(char *directory);
void create_z_grid(char *directory);
void expand_field(char *name,char *directory, int n);
void read_domain(char *directory);
void read_dims(char *directory);

int main(int argc, char *argv[]) {
    
    int n = atoi(argv[1]);
    nz = atoi(argv[2]);
    h = atof(argv[3]);
    printf("%d\t%s\n",n,argv[4]);

    read_dims(argv[4]);
    printf("Detected grid with nx = %d, ny = %d, nz = %d\n",nx,ny,nz);

    zm = (double *)malloc(sizeof(double)*(nz+1));
    ym = (double *)malloc(sizeof(double)*(ny+1));
    zmed = (double *)malloc(sizeof(double)*(nz));
    ymed = (double *)malloc(sizeof(double)*(ny));

    read_domain(argv[4]);

    create_z_grid(argv[4]);
    printf("Converting density\n");
    expand_field("dens",argv[4],n);
    printf("Converting energy\n");
    expand_field("energy",argv[4],n);

    printf("Converting grid\n");
    create_grid_file(argv[4]);

    free(ym); free(ymed);
    free(zm); free(zmed);
    return 1;
}

double scaleH(double r) {
    return h*r;
}

void create_grid_file(char *directory) {
    FILE *f;

    char fname[512];
    int i,j,k;

    double *xc = (double *)malloc(sizeof(double)*nx*(ny)*(nz));
    double *yc = (double *)malloc(sizeof(double)*nx*(ny)*(nz));
    double *zc = (double *)malloc(sizeof(double)*nx*(ny)*(nz));

    sprintf(fname,"%scart_grid.dat",directory);


#ifdef _OPENMP
#pragma omp parallel for collapse(3) private(i,j,k)
#endif
    for(k=0;k<nz;k++) {
        for(j=0;j<ny;j++) {
            for(i=0;i<nx;i++) {
               xc[l] = 
                   ymed[j]*cos((i+.5)*2*M_PI/nx);
               yc[l] = 
                   ymed[j]*sin((i+.5)*2*M_PI/nx);
               zc[l] = 
                   ymed[j]/tan(zmed[k]);
            }
        }
    }


    f = fopen(fname,"w");
    fwrite(xc,sizeof(double),nx*(ny)*(nz),f);
    fwrite(yc,sizeof(double),nx*(ny)*(nz),f);
    fwrite(zc,sizeof(double),nx*(ny)*(nz),f);
    fclose(f);

    free(xc); free(yc); free(zc);
    return;
}

void create_z_grid(char *directory) {
    FILE *f;
    char fname[512];
    int k;
    sprintf(fname,"%sdomain_z.dat",directory);


    f = fopen(fname,"r");
    if (f == NULL) printf("Error finding grid file %s\n",fname);

    for(k=0;k<nz+1;k++) {
        fprintf(f,"%lg\n",&zm[k]);
    }
    fclose(f);

    return;
}



void expand_field(char *name,char *directory, int n){
    FILE *f;
    char fname[512];
    int i,j,k;
    double zval,norm;
    double *fld, *fld2d;

    sprintf(fname,"%sgas%s%d.dat",directory,name,n);

    fld2d = (double *)malloc(sizeof(double)*nx*ny);
    fld = (double *)malloc(sizeof(double)*nx*ny*nz);


    f = fopen(fname,"r");

    fread(fld2d,sizeof(double),nx*ny,f);
    fclose(f);

#ifdef _OPENMP
#pragma omp parallel for collapse(3) private(i,j,k,zval,norm)
#endif
    for(k=0;k<nz;k++) {
        for(j=0;j<ny;j++) {
            for(i=0;i<nx;i++) {
                zval = 1./tan(zmed[k]);
                zval *= zval/2;
                norm = scaleH(ymed[j])*sqrt(2*M_PI);
                fld[l] = fld2d[l2d] * exp(-zval)/norm;
            }
        }
    }

    f = fopen(fname,"w");

    fwrite(fld,sizeof(double),nx*ny*nz,f);
    fclose(f);

    free(fld);
    free(fld2d);
    return;
}   

void read_domain(char *directory) {
    char fname[512];
    FILE *f;
    char tempstr[512];
    double temp;
    int i;

    sprintf(fname,"%sdomain_y.dat",directory);

    f = fopen(fname,"r");

    for(i=0;i<ny+NGHY+1;i++) {
        if (i >= NGHY) {
            fscanf(f,"%lg\n",&ym[i-NGHY]);
        }

    }

    for(i=0;i<ny;i++) {
        ymed[i] = .5*(ym[i] + ym[i+1]);
    }


    for(i=0;i<nz+1;i++) {
        if (i <= nz/2) {
            zm[i] = M_PI/2 + -zmin + i*zmin/(nz/2);
        }
        else {
            zm[i] = M_PI/2 + i*zmin/(nz/2);
        }
    }
    for(i=0;i<nz;i++) {
        zmed[i] = .5*(zm[i]+zm[i+1]);
    }

    fclose(f);
    return;
}

void read_dims(char *directory) {
    char fname[512];
    FILE *f;
    char tempstr[512];
    double temp;
    int i;

    sprintf(fname,"%sdimensions.dat",directory);

    f = fopen(fname,"r");

    for(i=0;i<DIMSLENGTH;i++) {
        fscanf(f,"%s\n",tempstr);

    }

    for(i=0;i<DIMSLENGTH;i++) {
        fscanf(f,"%lg\t",&temp);
        if (i == DIMSNXINDEX) {
            nx = (int)temp;
        }
        if (i == DIMSNXINDEX+1) {
            ny = (int)temp;
        }
    }
    fclose(f);
    return;
}
