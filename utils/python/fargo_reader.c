#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ymed(j) Ymed[(j)]
#define xmed(i) Xmed[(i)]
#define ymin(j) Ymin[(j)]
#define xmin(i) Xmed[(i)]
#define zone_size_x(j,k) (dx*ymed(j))
#define zone_size_y(j,k) (ymin(j+1)-ymin(j))
#define SurfY(j,k) ymin(j)*dx
#define SurfX(j,k) (ymin(j+1)-ymin(j))
#define InvVol(j,k) 2/(dx*(ymin(j+1)*ymin(j+1) - ymin(j)*ymin(j)))

#define l   ((i)+(j)*(nx)+((k)*stride))

#define lxp (((i)<(nx-1)) ? ((l)+1) : ((l)-(nx-1)))
#define lxm (((i)>0) ? ((l)-1) : ((l)+(nx-1)))

#define ixm ((i)>0 ? ((i)-1) : nx-1)
#define ixp ((i)<nx-1 ? ((i)+1) : 0)

#define lyp ((l)+nx)
#define lym ((l)-nx)

#define lzp (l+stride)
#define lzm (l-stride)

typedef struct Parameters {
    double alpha;
    double mp;
    double a;
    double omf;
    double h;
    double flaringindex;
    double nuindex;
    double vrindex;
    double mdot;
    double soft;

} Parameters;


int nx, ny, nz;
int size_x, size_y, size_z;
int stride, pitch;
double *Ymin, *Xmin, *Ymed, *Xmed;
double *rho, *vx, *vy;
double *rhos, *mdot;
double *tauxx, *tauyy, *tauxy;
double *work, *work1;
double dx;

Parameters params;

void read_files(int n, char *directory); 
void read_domain(char *directory);
void read_param_file(char *filename); 
void add_boundary(void);
double Nu(double x);
double Cs(double x);
void set_tensor(void);
void set_mdot(void);
void set_slopes(double *q); 
void set_rhostar(void);




int main(int arc, char *argv[]) { 
  int n;
  int i,j,k;
  char *directory, *param_fname;

  n = atoi(argv[1]);
  strcpy(directory,argv[2]);
  strcpy(param_fname,argv[3]);

  read_param_file(param_fname);

  size_x = nx;
  size_y = ny+2;
  size_z = nz;
  stride = nx*ny;
  pitch = nx;
  dx = 2*M_PI/nx;

  Xmin  =     (double*)malloc(sizeof(double)*(size_x+1));
  Ymin  =     (double*)malloc(sizeof(double)*(size_y+1));
  Xmed  =     (double*)malloc(sizeof(double)*(size_x));
  Ymed  =     (double*)malloc(sizeof(double)*(size_y));

  vx =     (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  vy =     (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  rho = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  rhos = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  
  tauxx =     (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  tauyy =     (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  tauxy = (double*)malloc(sizeof(double)*size_x*size_y*size_z);

  work = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  work1 = (double*)malloc(sizeof(double)*size_x*size_y*size_z);

  for(k=0; k<size_z; k++) {
    for(j=0; j<size_y; j++) {
      for(i=0; i<size_x; i++) {
            vx[l] = 0;
            vy[l] = 0;
            rho[l] = 0;
            tauxx[l] = 0;
            tauyy[l] = 0;
            tauxy[l] = 0;
            work[l] = 0;
            work1[l] = 0;
      }
    }
  }

  read_domain(directory);
  read_files(n,directory);
  add_boundary();
  set_tensor();

}
double Nu(double x) {
    return params.alpha*params.h*params.h*pow(x,params.nuindex);
}
double Cs(double x) {
    return params.h*pow(x,params.flaringindex-0.5);
}

void set_tensor(void) {
    int i,j,k;
    double visc, viscm, div_v;
    k = 0;

    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {
            visc = Nu(ymed(j));
            viscm = Nu(ymin(j));

            div_v = (vx[lxp]-vx[l])*SurfX(j,k);
            div_v += (vy[lyp]*SurfY(j+1,k)-vy[l]*SurfY(j,k));
            div_v *= 2.0/3.0*InvVol(j,k);

            tauxx[l] = visc*rho[l]*(2.0*(vx[lxp]-vx[l])/zone_size_x(j,k) - div_v);
            tauxx[l] += visc*rho[l]*(vy[lyp]+vy[l])/ymed(j);
            tauyy[l] = visc*rho[l]*(2.0*(vy[lyp]-vy[l])/(ymin(j+1)-ymin(j)) - div_v);


            tauxy[l] = viscm*.25*(rho[l]+rho[lxm]+rho[lym]+rho[lxm-pitch])*((vy[l]-vy[lxm])/(dx*ymin(j)) + (vx[l]-vx[lym])/(ymed(j)-ymed(j-1))-.5*(vx[l]+vx[lym])/ymin(j)); //centered on left, inner vertical edge in z
            

        }
    }
    FILE *f= fopen("temp_files/tensor.dat","w");
    k=0;

    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {        
            fwrite(&tauxx[l],sizeof(double),1,f);
        }
    }
    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {        
            fwrite(&tauyy[l],sizeof(double),1,f);
        }
    }
    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {        
            fwrite(&tauxy[l],sizeof(double),1,f);
        }
    }
    fclose(f);

    return;
}

void set_slopes(double *q) {
    int i,j,k;
    k=0;
    double dqm, dqp;

    for (j=1; j<size_y-1; j++) {
      for (i=0; i<size_x; i++) {

        dqm = (q[l]-q[lym])/zone_size_y(j,k);
        dqp = (q[lyp]-q[l])/zone_size_y(j+1,k);
        if(dqp*dqm<=0) work[l] = 0;
        else  work[l] = 2.*dqp*dqm/(dqm+dqp);
      }
    }
    
    return;
}
void set_rhostar(void) {
    int i,j,k;
    FILE *f;
    k = 0;
    set_slopes(rho);
    f = fopen("temp_files/rhostar.dat","w");

    for (j=1; j<size_y-1; j++) {
      for (i=0; i<size_x; i++) {
    	if(vy[l]>0.)
	  rhos[l] = rho[lym] + 0.5 * (zone_size_y(j-1,k))*work[lym];
	else
	  rhos[l] = rho[l] - 0.5 * (zone_size_y(j,k))*work[l];
        
        fwrite(&rhos[l],sizeof(double),1,f);
      }

    }
    fclose(f);
    return;
}

void set_mdot(void) {
    int i,j,k;
    k = 0;
    FILE *f;


    f = fopen("temp_files/mdot.dat","w");

    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {
            
	        mdot[l] = -nx*vy[l]*rhos[l]*SurfY(j,k) ;
            fwrite(&mdot[l],sizeof(double),1,f);
        }
    }
    fclose(f);
    return;

}


void lam_ex(void) {
    int i,j,k;
    FILE *f;
    double smoothing = params.soft*params.h*pow(params.a,1 + params.flaringindex);
    smoothing *= smoothing;
    double rad,res;
    double pot;

    f = fopen("temp_files/lambda_ex.dat","w");

    k = 0;
    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {    
            rad = ymed(j)*ymed(j) + params.a*params.a -2*params.a*ymed(j)*cos(xmed(i))+smoothing;
            rad *= sqrt(rad);
            pot = params.mp * ymed(j)*params.a*sin(xmed(i))/rad;
            res = 2*M_PI*ymed(j)*rho[l]*pot;
            fwrite(&res,sizeof(double),1,f);

        }
    }
    fclose(f);
    return;

}

void set_Fd(void) {
    FILE *f;
    int i,j,k;
    k=0;

    double *mdbar = (double *)malloc(sizeof(double)*(ny+2+1));
    double *lbar = (double *)malloc(sizeof(double)*(ny+2));
    double res1, res2;
    for(j=1;j<size_y-1;j++) {
        res1=0;
        res2=0;

        



}



void add_boundary(void) {
    double fh1;
    double fac,rhom;
    double facm;
    int i,j,k;
    j  = 0;
    for(k=0; k<size_z; k++) {
        for(i=0;i<size_x;i++) {
            rho[l] = rho[lyp]*pow(ymed(j)/ymed(j+1),-params.nuindex);
            vx[l] = vx[lyp]*pow(ymed(j)/ymed(j+1),-.5);
            vy[l] = vy[lyp]*pow(ymin(j)/ymin(j+1),params.vrindex);

        }
    }
    j=size_y-1;
    for(k=0; k<size_z; k++) {
        for(i=0;i<size_x;i++) {
            fh1 = rho[lym]*3*M_PI*Nu(ymed(j-1))*sqrt(ymed(j-1));
            fac = 3*M_PI*Nu(ymed(j))*sqrt(ymed(j));
            facm = 3*M_PI*Nu(ymin(j))*sqrt(ymin(j));
            
            rho[l] = (fh1 + params.mdot*(sqrt(ymed(j))-sqrt(ymed(j-1)))/(3*M_PI))/fac;
            rhom = (fh1 + params.mdot*(sqrt(ymin(j))-sqrt(ymed(j-1)))/(3*M_PI))/facm;
            vx[l] = params.mdot/(-2*M_PI*ymin(j)*rhom);
            vy[l] = vy[lym]*pow(ymin(j)/ymin(j-1),params.vrindex);

        }
    }

    return;

}

void read_param_file(char *filename) {
    FILE *f;
    
    f = fopen(filename,"r");
    fscanf(f,"%d\n",&nx);
    fscanf(f,"%d\n",&ny);
    fscanf(f,"%d\n",&nz);
    fscanf(f,"%lg\n",&params.alpha);
    fscanf(f,"%lg\n",&params.mp);
    fscanf(f,"%lg\n",&params.a);
    fscanf(f,"%lg\n",&params.omf);
    fscanf(f,"%lg\n",&params.h);
    fscanf(f,"%lg\n",&params.flaringindex);
    fscanf(f,"%lg\n",&params.mdot);

    params.nuindex = 2*params.flaringindex + 0.5;
    params.vrindex = params.nuindex - 1.0;

    fclose(f);

    return;
}

void read_domain(char *directory) {
    FILE *fx, *fy;
    char filename[512];
    double temp;
    int i,j;
    sprintf(filename,"%sdomain_x.dat",directory);
    fx = fopen(filename,"r");
    sprintf(filename,"%sdomain_y.dat",directory);
    fy = fopen(filename,"r");

    for(i=0;i<size_x+1;i++) {
        fscanf(fx,"%lg\n",&Xmin[i]);
    }
    for(j=0;j<ny+6;j++) {
        if ( (j < 2) || j >= (ny+1+6-2) ) {
            fscanf(fy,"%lg\n",temp);
        }
        else {
            fscanf(fy,"%lg\n",&Ymin[j]);
        }
    }
    
    for(i=0;i<size_x;i++) {
        Xmed[i] = .5*(Xmin[i] + Xmin[i+1]);
    }

    for(j=0;j<size_y;j++) {
        Ymed[j] = .5*(Ymin[j] + Ymin[j+1]);
    }

    fclose(fx);
    fclose(fy);
    return;
}

void read_files(int n, char *directory) {
    char filename[512];
    FILE *fd,*fx,*fy;
    int i,j,k;
    sprintf(filename,"%sgasdensity%d.dat",directory,n);
    fd = fopen(filename,"r");
    sprintf(filename,"%sgasvx%d.dat",directory,n);
    fx = fopen(filename,"r");
    sprintf(filename,"%sgasvy%d.dat",directory,n);
    fy = fopen(filename,"r");

    for(k=0;k<size_z;k++) {
        for(j = 1; j<size_y-1;j++) {
            for(i=0;i<size_x;i++) {
                fread(&rho[l],sizeof(double),1,fd);
                fread(&vx[l],sizeof(double),1,fx);
                fread(&vy[l],sizeof(double),1,fy);
            }
        }
    }
    fclose(fd);
    fclose(fx);
    fclose(fy);
    return;

}




