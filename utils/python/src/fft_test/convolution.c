#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <rfftw.h>


#define MMAX 50
#define NGHY 0
#define ymed(j) Ymed[(j)]
#define xmed(i) Xmed[(i)]
#define ymin(j) Ymin[(j)]
#define xmin(i) Xmed[(i)]
#define zone_size_x(j,k) (dx*ymed(j))
#define zone_size_y(j,k) (ymin(j+1)-ymin(j))
#define SurfY(j,k) ymin(j)*dx
#define SurfX(j,k) (ymin(j+1)-ymin(j))
#define InvVol(j,k) 2/(dx*(ymin(j+1)*ymin(j+1) - ymin(j)*ymin(j)))
#define Vol(j,k) 0.5*(dx*(ymin(j+1)*ymin(j+1) - ymin(j)*ymin(j)))
#define XC ymed(j)*cos(xmed(i))
#define YC ymed(j)*sin(xmed(i))

#define l   ((i)+(j)*(nx)+((k)*stride))
#define l_f(ii,jj,kk)   ((ii)+(jj)*(nx)+((kk)*stride))
#define lact   ((i)+(jact)*(nx)+((k)*stride))
#define lxp (((i)<(nx-1)) ? ((l)+1) : ((l)-(nx-1)))
#define lxm (((i)>0) ? ((l)-1) : ((l)+(nx-1)))
#define l2D ((j)+((k)*pitch2d))
#define l2D_int ((j)+((k)*pitch2d))

#define ixm ((i)>0 ? ((i)-1) : nx-1)
#define ixp ((i)<nx-1 ? ((i)+1) : 0)

#define lyp ((l)+nx)
#define lym ((l)-nx)

#define lzp (l+stride)
#define lzm (l-stride)


void read_files(void);
void read_domain(void);
void allocate(void);

int nx,ny,size_x,size_y,stride,pitch,pitch2d,pitch2d_int;
int size_f;
int nz,size_z;
double dx;
double *Xmin, *Ymin, *Xmed, *Ymed;
fftw_real *in1, *in2;
double *dens, *vy;
double *res;
fftw_real *out1,*out2;
rfftw_plan p1;

int main(int argc, char *argv[]) {

    ny = 512;
    nx = 700;
    nz = 1; size_z = nz;
    size_x = nx;
    size_y = ny + 2*NGHY;
    size_f = (size_x+1)/2;
    stride = size_x*size_y;
    pitch = size_x;
    pitch2d = 0;
    pitch2d_int = 0 ;//pitch*(int)double_to_int;
    dx = 2*M_PI/nx;
    
    allocate();
    read_domain();
    read_files();

    printf("Read domain\n");

    p1 = rfftw_create_plan(size_x,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE);
    printf("Allocated plan\n");

    int i,j,k;
    i = j = k =0;
    double fac = (double)(size_x * size_x);
    for(j=0;j<size_y;j++) {
        i = 0;
        memcpy(in1,&dens[l],sizeof(double)*size_x);
        memcpy(in2,&vy[l],sizeof(double)*size_x);
        rfftw_one(p1,in1,out1);
        rfftw_one(p1,in2,out2);

        res[j*size_f] = out1[0]*out2[0]/fac;
        for(i=1;i< size_f ;i++) {
            res[i + j*size_f] = 2 *( out1[i]*out2[i] + out1[size_x-i]*out2[size_x-i] )/fac;
        }
    }

    

    printf("%d\t%d\t%d\n",size_y,(size_x+1)/2,size_y*(size_x+1)/2);
    FILE *f = fopen("output.dat","w");
    fwrite(res,sizeof(double),size_y*MMAX,f);
    fclose(f);

    rfftw_destroy_plan(p1);
    fftw_free(out1); fftw_free(out2);

    return 0;
}

void read_domain(void) {
    FILE *fx, *fy;
    char filename[512];
    char filename2[512];
    int i,j;

    fx = fopen("domain_x.dat","r");

    for(i=0;i<size_x+1;i++) {
        fscanf(fx,"%lg\n",&Xmin[i]);
    }
    fclose(fx);
    fx = fopen("domain_y.dat","r");

    for(i=0;i<size_y+1;i++) {
        fscanf(fx,"%lg\n",&Ymin[i]);
    }
    fclose(fx);

    for(i=0;i<size_x;i++) {
        Xmed[i] = .5*(Xmin[i] + Xmin[i+1]);
    }
    for(j=0;j<size_y;j++) {
        Ymed[j] = .5*(Ymin[j] + Ymin[j+1]);
    }

    return;
}

void read_files(void) {
    FILE *fd,*fy;
    int i,j,k;


    fd = fopen("gasdens30.dat","r");
    fy = fopen("gasvy30.dat","r");

    i=j=k=0;
    for(k=0;k<size_z;k++) {
        for(j =NGHY; j<size_y-NGHY;j++) {
            for(i=0;i<size_x;i++) {
                fread(&dens[l],sizeof(double),1,fd);
                fread(&vy[l],sizeof(double),1,fy);
            }
        }
    }
    fclose(fd);
    fclose(fy);

    return;

}
void allocate(void) {
    Ymin =(double *) malloc(sizeof(double)*(size_y+1));
    Xmin =(double *) malloc(sizeof(double)*(size_x+1));
    Ymed =(double *) malloc(sizeof(double)*size_y);
    Xmed =(double *) malloc(sizeof(double)*size_x);
    dens =(double *) malloc(sizeof(double)*size_y*size_x);
    vy =  (double *) malloc(sizeof(double)*size_y*size_x);
    in1 =  (fftw_real *) malloc(sizeof(fftw_real)*size_x);
    out1 = (fftw_real *) fftw_malloc(sizeof(fftw_real)*size_x);
    in2 =  (fftw_real *) malloc(sizeof(fftw_real)*size_x);
    out2 = (fftw_real *) fftw_malloc(sizeof(fftw_real)*size_x);
    res = (double *)malloc(sizeof(double)*size_y*size_f);
    return;


}
