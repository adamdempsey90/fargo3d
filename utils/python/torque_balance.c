#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NGHY 3
#define TRUE 1
#define FALSE 0
#define MALLOC_SAFE(ptr) if (ptr == NULL) printf("Malloc error at line %d!\n",__LINE__);

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

#define l   ((i)+(j)*(nx)+((k)*stride))
#define l_f(ii,jj,kk)   ((ii)+(jj)*(nx)+((kk)*stride))
#define lact   ((i)+(jact)*(nx)+((k)*stride))
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


double *dens, *vx, *vy, *Pres, *Pot;
double *dbar, *vxbar, *vybar, *dbarstar;
double *Ymed, *Xmed, *Ymin, *Xmin;
double *tauxx, *tauxy, *tauyy;
double *Pixp, *Pixm;
double *slope, *divrho, *denstar, *Qs;
double *Lstar, *Ldstar;
double *Lang,*drF,*Lambda;
double *Ld, *drFd, *Lamdep;

double dt,omf,dx;
int nx, ny, nz, size_x, size_y, size_z,stride,pitch;

Parameters params;

double Nu(double x);
double Cs(double x);
void allocate_all(void);
void free_all(void);
void viscosity(void);
void source_step(void);
void transport_step(void);
void set_momenta(void);
void transportY(void);
void DividebyRho(double *q);
void vanleer_y_a(double *q);
void vanleer_y_b(double *q, double *qs, double dt);
void updateY(double *q, double *qs);
void source_step_avg(void);
void transport_step_avg(void);
void set_momenta_avg(void);
void transportY_avg(void);
void DividebyRho_avg(double *q);
void vanleer_y_a_avg(double *q);
void vanleer_y_b_avg(double *q, double *qs, double dt);
void updateY_avg(double *q, double *qs);
void set_averages(void);
void set_bc(void);
void ymax_bound(void);
void ymin_bound(void);
void read_param_file(char *filename);
void read_domain(char *directory);
void read_files(int n, char *directory);
void output(char *directory);
void output_avg(char *directory);
void set_lamdep(void);

int main(int argc, char *argv[]) {

    int n,nsteps;
    char directory[256];
    char param_fname[100];

    n = atoi(argv[1]);
    nsteps = atoi(argv[2]);
    strcpy(directory,argv[3]);
    sprintf(param_fname,"%sparam_file.txt",directory);


    printf("Reading param file %s\n",param_fname);
    read_param_file(param_fname);

    size_x = nx;
    size_y = ny+2*NGHY;
    size_z = nz;
    stride = size_x*size_y;
    pitch = size_x;
    dx = 2*M_PI/nx;
    dt = M_PI*1e-5;

    allocate_all();
    read_domain(directory);
    read_files(n,directory);
    set_bc();
    set_averages(); 
    viscosity();
    source_step();
    transportY();

    source_step_avg();
    transportY_avg();
    set_lamdep();


    output(directory);
    output_avg(directory);
    free_all();
    return 0;
}
double Nu(double x) {
    return params.alpha*params.h*params.h*pow(x,params.nuindex);
}
double Cs(double x) {
    return params.h*pow(x,params.flaringindex-0.5);
}
void allocate_all(void) {
    
    MALLOC_SAFE((Ymin = (double *)malloc(sizeof(double)*(size_y+1))));
    MALLOC_SAFE((Xmin = (double *)malloc(sizeof(double)*(size_x+1))));

    MALLOC_SAFE((Ymed = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((Xmed = (double *)malloc(sizeof(double)*(size_x))));
    
    MALLOC_SAFE((dens = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((vx = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((vy = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Lang = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((drF = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Lambda = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Pres = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Pot = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Pixm = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Pixp = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((denstar = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Lstar = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Qs = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((slope = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((divrho = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((tauxx = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((tauxy = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((tauyy = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));

    MALLOC_SAFE((dbar = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((vxbar = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((vybar = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((Ld = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((drFd = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((Lamdep = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((dbarstar = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((Ldstar = (double *)malloc(sizeof(double)*(size_y))));

    return;
}
void free_all(void) {
    
    free(Ymin);
    free(Xmin);

    free(Ymed);
    free(Xmed);
    
    free(dens);
    free(vx);
    free(vy);
    free(Lang);
    free(drF);
    free(Lambda);
    free(Pres);
    free(Pot);
    free(Pixm);
    free(Pixp);
    free(denstar);
    free(Lstar);
    free(Qs);
    free(slope);
    free(divrho);
    free(tauxx);
    free(tauxy);
    free(tauyy);

    free(dbar);
    free(vxbar);
    free(vybar);
    free(Ld);
    free(drFd);
    free(Lamdep);
    free(dbarstar);
    free(Ldstar);
    return;
}
void compute_Pres(void) {
    int i,j,k;
    k=0;
    double cs2;
    for(j=0;j<size_y;j++) {
        cs2 = Cs(ymed(j));
        cs2 *= cs2;
        for(i=0;i<size_x;i++) {
            Pres[l] = cs2*dens[l];
        }
    }
    return;
}
void compute_potential(void) {
    int i,j,k;
    k=0;
    double xpl, ypl, xd, yd;
    double smoothing,distp,rad;
    xpl = params.a;
    ypl = 0;
    distp = sqrt(xpl*xpl + ypl*ypl);

	smoothing = params.h*pow(distp,params.flaringindex)*distp*params.soft;
    smoothing *= smoothing;
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            yd = ymed(j)*sin(xmed(i));
            xd = ymed(j)*cos(xmed(i));
            rad = (xd-xpl)*(xd-xpl) + (yd-ypl)*(yd-ypl);
            Pot[l] = -1./ymed(j);
            Pot[l] -= params.mp/sqrt(rad + smoothing);
        }
    }
    return;
}
void viscosity(void) {
    int i,j,k;
    k=0;
    double visc, viscm,div_v;
    double vf,vfm;
    for(j=1;j<size_y-1;j++) {
        visc = Nu(ymed(j));
        viscm = Nu(ymin(j));
        vf = 0;//omf*ymed(j);
        vfm = 0;//omf*ymed(j-1);
        for(i=0;i<size_x;i++) {

            div_v = 0.0;
            div_v += (vx[lxp]-vx[l])*SurfX(j,k);
            div_v += (vy[lyp]*SurfY(j+1,k)-vy[l]*SurfY(j,k));
            div_v *= 2.0/3.0*InvVol(j,k);

            tauxx[l] = visc*dens[l]*(2.0*(vx[lxp]-vx[l])/zone_size_x(j,k) - div_v);
            tauxx[l] += visc*dens[l]*(vy[lyp]+vy[l])/ymed(j);
            tauyy[l] = visc*dens[l]*(2.0*(vy[lyp]-vy[l])/(ymin(j+1)-ymin(j)) - div_v);
            tauxy[l] = viscm*.25*(dens[l]+dens[lxm]+dens[lym]+dens[lxm-pitch])*((vy[l]-vy[lxm])/(dx*ymin(j)) + (vx[l]+vf-(vx[lym]+vfm))/(ymed(j)-ymed(j-1))-.5*(vx[l]+vf+vx[lym]+vfm)/ymin(j)); //centered on left, inner vertical edge in z


        }
    }

    return;
}
void source_step(void) {
    int i,j,k;
    k=0;
    compute_Pres();
    compute_potential();

    for(j=1;j<size_y-1;j++) {

        for(i=0;i<size_x;i++) {
            // X
            drF[l] =  ymed(j)*(Pres[lxp]-Pres[lxm])/(2*zone_size_x(j,k));
            Lambda[l] = -.5*ymed(j)*(dens[l]+dens[lxm ])*.5*(Pot[l] - Pot[lxm])/(zone_size_x(j,k));
            Lambda[l] += -.5*ymed(j)*(dens[l]+dens[lxp])*.5*(Pot[lxp] - Pot[l])/(zone_size_x(j,k));
  
            drF[l] += -.5*ymed(j)*(tauxx[lxp]-tauxx[lxm])/zone_size_x(j,k);
            drF[l] += -ymed(j)*(ymin(j+1)*ymin(j+1)*.5*(tauxy[lyp]+tauxy[lxp+pitch]) - ymin(j)*ymin(j)*.5*(tauxy[l]+tauxy[lxp]))/(ymed(j)*ymed(j)*(ymin(j+1)-ymin(j)));
       }
    }


    return;
}
void source_step_avg(void) {
    int i,j,k;
    k=0;
    double res, resp;

    for(j=1;j<size_y-1;j++) {
        resp = 0;
        res = 0;
        for(i=0;i<size_x;i++) {
            resp += tauxy[lyp] + tauxy[lxp+pitch];
            res += tauxy[l] + tauxy[lxp];
        }
        resp  /= (double)nx;
        res /= (double)nx;
        drFd[j] = -ymed(j)*(ymin(j+1)*ymin(j+1)*.5*resp - ymin(j)*ymin(j)*.5*res)/(ymed(j)*ymed(j)*(ymin(j+1)-ymin(j)));
    }
    return;
}
void set_lamdep(void) {
    int i,j,k;
    k=0;
    double res, resp;

    for(j=1;j<size_y;j++) {
        res1 = 0; // dr(r * sig*vy)
        res2 = 0; // dr(vyp * l)-lp*dr*vyp 
        res3 = 0; // <sig/dens>
        res4 = 0; // dr (tau)
        for(i=0;i<size_x;i++) {
            res1 -=InvVol(j,k) *( SurfY(j+1,k)*(vy[lyp]-vybar[j+1])*(denstar[lyp] - dbarstar[j+1]) - Surfy(j,k)*(vy[l]-vybar[j])*(denstar[l]-dbarstar[j]));
            res2 -= InvVol(j,k) * ( SurfY(j+1,k)*(vy[lyp]-vybar[j+1])*(Lstar[lyp]-Ldstar[j+1]) - SurfY(j,k)*(vy[l]-vybar[j])*(Lstar[l]-Ldstar[j]));
            res2 += (Lang[l]-Ld[j])*InvVol(j,k) *( SurfY(j+1,k)*(vy[lyp]-vybar[j+1])-Surfy(j,k)*(vy[l]-vybar[j]));

            res3 += (dens[l] - dbar[j])/dens[l];
            res4 += 

        }
        res1 *= Ld[j]/(double)nx;
        res2 *= dbar[j]/(double)nx;
        res3  /= (double)nx;
        res4 *= dbar[j]/(double)nx;



        
    }


    return;
}
void set_momenta_avg(void) {
    int i,j,k;
    k=0;

    for(j=1;j<size_y-1;j++) {
        Pixm[j] = ymed(j)*(vxbar[j] + omf*ymed(j))*dens[j];
        Lang[j] = Pixm[j];
    }
    return;
}
void set_momenta(void) {
    int i,j,k;
    k=0;

    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {
            Pixm[l] = ymed(j)*(vx[l]+omf*ymed(j))*dens[l];
            Pixp[l] = ymed(j)*(vx[lxp]+omf*ymed(j))*dens[l];
            Lang[l] = .5*(Pixm[l] + Pixp[l]);
        }
    }
    return;
}
void DividebyRho_avg(double *q) {
    int j;
    for(j=0;j<size_y;j++) {
        divrho[j] = q[j]/dbar[j];
    }
    return;
}
void DividebyRho(double *q) {
    int i,j,k;
    k=0;
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            divrho[l] = q[l]/dens[l];
        }
    }
    return;
}
void transportY(void) {
    // Y direction
    int i,j,k;
    k=0;
    dt = 0;
    set_momenta();
    vanleer_y_a(dens);
    vanleer_y_b(dens,denstar,dt);
    
    DividebyRho(Pixm);
    vanleer_y_a(divrho);
    vanleer_y_b(divrho,Qs,dt);
    
    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {
            Lstar[l] = .5*Qs[l];
        }
    }


    updateY(Pixm,Qs);

    DividebyRho(Pixp);
    vanleer_y_a(divrho);
    vanleer_y_b(divrho,Qs,dt);
    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {
            Lstar[l] += .5*Qs[l];
        }
    }
    updateY(Pixp,Qs);

    return;

}
void transportY_avg(void) {
    // Y direction
    int j;
    dt = 0;
    set_momenta_avg();
    vanleer_y_a_avg(dbar);
    vanleer_y_b_avg(dbar,dbarstar,dt);
    
    DividebyRho_avg(Pixm);
    vanleer_y_a_avg(divrho);
    vanleer_y_b_avg(divrho,Qs,dt);
    for(j=1;j<size_j-1;j++) {
        Ldstar[j] = Qs[j];
    }
    updateY_avg(Pixm,Qs);


    return;

}
void updateY_avg(double *q, double *qs) {
    int i,j,k;
    k = 0;
    for(j=0;j<size_y-1;j++) {
        drFd[j] -= .5*((vybar[j]*qs[j]*dbarstar[j]*SurfY(j,k)-vybar[j+1]*qs[j+1]*dbarstar[j+1]*SurfY(j+1,k))*InvVol(j,k));
    }
    return;
}
void updateY(double *q, double *qs) {
    int i,j,k;
    k = 0;
    for(j=0;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {

            drF[l] -= .5*((vy[l]*qs[l]*denstar[l]*SurfY(j,k)-vy[lyp]*qs[lyp]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));

        }
    }
    return;
}
void vanleer_y_a_avg(double *q) {
    int i,j,k;
    k=0;
    double dqm, dqp;
    for(j=1;j<size_y-1;j++) {
        
        dqm = (q[j] - q[j-1])/zone_size_y(j,k);
        dqp = (q[j+1]-q[j])/zone_size_y(j,k);

        if (dqp*dqm <= 0) {
            slope[j] = 0;
        }
        else {
            slope[j] = 2*dqp*dqm/(dqp+dqm);
        }
    }
    
    return;
}
void vanleer_y_a(double *q) {
    int i,j,k;
    k=0;
    double dqm, dqp;
    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {

            dqm = (q[l] - q[lym])/zone_size_y(j,k);
            dqp = (q[lyp]-q[l])/zone_size_y(j,k);

            if (dqp*dqm <= 0) {
                slope[l] = 0;
            }
            else {
                slope[l] = 2*dqp*dqm/(dqp+dqm);
            }
        }
    }
    
    return;
}
void vanleer_y_b_avg(double *q, double *qs, double dt) {
    int i,j,k;
    k=0;

    for(j=1;j<size_y-1;j++) {
            
        if (vybar[j] > 0.) {
            qs[j] = q[j-1] + .5*(zone_size_y(j-1,k)-vybar[j]*dt)*slope[j-1];
        }
        else {
            qs[j] = q[j] - .5*(zone_size_y(j,k)+vybar[j]*dt)*slope[j];
        }

    }
    return;
}
void vanleer_y_b(double *q, double *qs, double dt) {
    int i,j,k;
    k=0;

    for(j=1;j<size_y-1;j++) {

        for(i=0;i<size_x;i++) {
            
            if (vy[l] > 0.) {
                qs[l] = q[lym] + .5*(zone_size_y(j-1,k)-vy[l]*dt)*slope[lym];
            }
            else {
                qs[l] = q[l] - .5*(zone_size_y(j,k)+vy[l]*dt)*slope[l];
            }

        }
    }
    return;
}
void set_bc(void) {
    
    ymin_bound();
    ymax_bound();
    return;

}

void set_averages(void) {
    int i,j,k;
    k=0;
    double resd, resx, resy;
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            resd = 0;
            resx = 0;
            resy = 0;
            for(i=0;i<size_y;i++) {
                resd += dens[l];
                resx += vx[l];
                resy += vy[l];
            }
            dbar[j] = resd/(double)nx;
            vxbar[j] = resx/(double)nx;
            vybar[j] = resy/(double)nx;
        }
    }
    return;
}
void ymax_bound(void) {
    int i,j,k;

  int jact;
  int jgh;
  int lgh;
  int lghs;
  int lactb;
  int lactbs;
  int lactbs_null;
  for(k=0; k<size_z; k++) {
    for(j=0; j<NGHY; j++) {
      for(i=0; i<size_x; i++) {

	lgh = i + (ny+NGHY+j)*pitch + k*stride;
	lghs = i + (ny+NGHY+1+j)*pitch + k*stride;
	lactb = i + (ny+NGHY-1-j)*pitch + k*stride;
	lactbs = i + (ny+NGHY-1-j)*pitch + k*stride;
	lactbs_null = i + (ny+NGHY)*pitch + k*stride;
	jgh = (ny+NGHY+j);
	jact = (ny+NGHY-1-j);

	dens[lgh] = (dens[i+(ny+NGHY-1)*pitch]-dens[i+(ny+NGHY-2)*pitch])/(ymed(ny+NGHY-1)-ymed(ny+NGHY-2))*(ymed(jgh)-ymed(ny+NGHY-1))+dens[i+(ny+NGHY-1)*pitch];
	vx[lgh] = (vx[lactb]+ymed(jact)*omf)*sqrt(ymed(jact)/ymed(jgh)) - ymed(jgh)*omf;
	if (j<size_y-1)
		vy[lghs] = -vy[lactbs];
	vy[lactbs_null] = 0;
      }
    }
  }

    return;
}

void ymin_bound(void) {

  int i;
  int j;
  int k;
  int jact;
  int jgh;
  int lgh;
  int lghs;
  int lactb;
  int lactbs;
  int lactbs_null;
for(k=0; k<size_z; k++) {
    for(j=0; j<NGHY; j++) {
      for(i=0; i<size_x; i++) {

	lgh = l;
	lghs = l;
	lactb = i + (2*NGHY-j-1)*pitch + k*stride;
	lactbs = i + (2*NGHY-j)*pitch + k*stride;
	lactbs_null = i + NGHY*pitch + k*stride;
	jgh = j;
	jact = (2*NGHY-j-1);

	dens[lgh] = (dens[i+(NGHY+1)*pitch]-dens[i+NGHY*pitch])/(ymed(NGHY+1)-ymed(NGHY))*(ymed(jgh)-ymed(NGHY))+dens[i+NGHY*pitch];
	vx[lgh] = (vx[lactb]+ymed(jact)*omf)*sqrt(ymed(jact)/ymed(jgh)) - ymed(jgh)*omf;
	vy[lghs] = -vy[lactbs];
	vy[lactbs_null] = 0;
//<\#>
      }
    }
  }
    return;
}
void read_param_file(char *filename) {
    FILE *f;
    
    printf("Reading %s\n",filename);
    f = fopen(filename,"r");
    if (f == NULL) {
        printf("Can't find parameter file, %s\n",filename);
        exit(0);
    }
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
    fscanf(f,"%lg\n",&params.soft);

    params.nuindex = 2*params.flaringindex + 0.5;
    params.vrindex = params.nuindex - 1.0;
    omf = params.omf;
    fclose(f);
    printf("nx=%d\tny=%d\tnz=%d\n",nx,ny,nz);
    printf("alpha=%.1e\tmp=%.1e\ta=%lg\n",params.alpha,params.mp,params.a);
    printf("omf=%lg\th=%lg\tflaring=%lg\n",params.omf,params.h,params.flaringindex);
    printf("mdot=%.2e\tsoft=%lg\n",params.mdot,params.soft);
    return;
}

void read_domain(char *directory) {
    FILE *fx, *fy;
    char filename[512];
    char filename2[512];
    int i,j;

    sprintf(filename,"%sdomain_x.dat",directory);
    printf("Reading %s\n",filename);
    fx = fopen(filename,"r");
    if (fx == NULL) printf("Error reading %s\n",filename);

    for(i=0;i<size_x+1;i++) {
        fscanf(fx,"%lg\n",&Xmin[i]);
    }
    fclose(fx);
    sprintf(filename2,"%sdomain_y.dat",directory);
    printf("Reading %s\n",filename2);
    fy = fopen(filename2,"r");
    if (fy == NULL) printf("Error reading %s\n",filename);
    for(j=0;j<size_y+1;j++) {
            fscanf(fy,"%lg\n",&Ymin[j]);
    }
    fclose(fy);
   

    for(i=0;i<size_x;i++) {
        Xmed[i] = .5*(Xmin[i] + Xmin[i+1]);
    }
    for(j=0;j<size_y;j++) {
        Ymed[j] = .5*(Ymin[j] + Ymin[j+1]);
    }

    return;
}

void read_files(int n, char *directory) {
    char filename[512];
    FILE *fd,*fx,*fy;
    int i,j,k;
    sprintf(filename,"%sgasdens%d.dat",directory,n);
    printf("Reading %s\n",filename);
    fd = fopen(filename,"r");
    if (fd == NULL) printf("Error loading %s\n",filename); 
    sprintf(filename,"%sgasvx%d.dat",directory,n);
    printf("Reading %s\n",filename);
    fx = fopen(filename,"r");
    if (fx == NULL) printf("Error loading %s\n",filename); 
    sprintf(filename,"%sgasvy%d.dat",directory,n);
    printf("Reading %s\n",filename);
    fy = fopen(filename,"r");
    if (fy == NULL) printf("Error loading %s\n",filename); 
    for(k=0;k<size_z;k++) {
        for(j =NGHY; j<size_y-NGHY;j++) {
            for(i=0;i<size_x;i++) {
                fread(&dens[l],sizeof(double),1,fd);
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
void output_avg(char *directory) {
    int i,j,k;
    k=0;
    double res;
    FILE *f;
    char fname[256];
    sprintf(fname,"%stemp_files/torques.dat",directory);
    f = fopen(fname,"w");
    fwrite(&Ymin[0],sizeof(double),size_y+1,f);

    for(j=0;j<size_y;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
            res += Lang[l];
        }
        res /= (double)nx;
        fwrite(&res,sizeof(double),1,f);
    }
    for(j=0;j<size_y;j++) {
        res = Ld[j];
        fwrite(&res,sizeof(double),1,f);
    }
    for(j=0;j<size_y;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
            res += Lang[l];
        }
        res /= (double)nx;
        res -= Ld[j];
        fwrite(&res,sizeof(double),1,f);
    }

    for(j=0;j<size_y;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
            res += drF[l];
        }
        res /= (double)nx;
        fwrite(&res,sizeof(double),1,f);
    }
    for(j=0;j<size_y;j++) {
        res = drFd[j];
        fwrite(&res,sizeof(double),1,f);
    }
    for(j=0;j<size_y;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
            res += drF[l];
        }
        res /= (double)nx;
        res -= drFd[j];
        fwrite(&res,sizeof(double),1,f);
    }

    for(j=0;j<size_y;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
            res += Lambda[l];
        }
        res /= (double)nx;
        fwrite(&res,sizeof(double),1,f);
    }
    for(j=0;j<size_y;j++) {
        res = Lamdep[j];
        fwrite(&res,sizeof(double),1,f);
    }

    fclose(f);
    return;
}
void output(char *directory) {
    FILE *f;
    char fname[256];
    sprintf(fname,"%stemp_files/output2.dat",directory);
    f = fopen(fname,"w");
    fwrite(&Ymin[0],sizeof(double),size_y+1,f);
    fwrite(&Xmin[0],sizeof(double),size_x+1,f);
    fwrite(&dens[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vy[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vx[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&tauxy[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&Lang[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&drF[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&Lambda[0],sizeof(double),size_x*size_y*size_z,f);
    fclose(f);
    return;
}
