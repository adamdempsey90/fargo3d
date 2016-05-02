#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NGHY 3
#define TRUE 1
#define FALSE 0
#define G 1.0
#define MSTAR 1.0

#define MALLOC_SAFE(ptr) if (ptr == NULL) printf("Malloc error at line %d!\n",__LINE__);
#define FREE_SAFE(ptr) free(ptr); ptr=NULL;

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

double *drFt, *drFd, *drFw, *Ld, *Lt, *Lw, *Lamex, *Lamdep;
double *dbarstar;
double *dens, *vx, *vy, *Pres, *Pot, *energy;
double *Pixp, *Pixm;
double *slope, *divrho, *denstar, *Qs;
double *Ymed, *Xmed, *Ymin, *Xmin;
double *tauxx, *tauxy, *tauyy, *tauxyavg;

double dt,omf,dx;
int nx, ny, nz, size_x, size_y, size_z,stride,pitch;

Parameters params;

double Nu(double x);
double Cs(double X);


void set_ang(void);
void allocate_all(void);
void free_all(void);
void compute_Pres(void);
void compute_energy(void);
void viscosity(void);
void potential(void);
void temp_to_vel(void);
void vel_to_temp(void);
void source_step(void);
void transport_step(void);
void set_momenta(void);
void transportY(void);
void DividebyRho(double *q);
void transportX(void);
void set_vel(void);
void vanleer_y_a(double *q);
void vanleer_x_a(double *q);
void vanleer_y_b(double *q, double *qs, double dt);
void vanleer_x_b(double *q, double *qs, double dt);
void updateX(double *q, double *qs,double dt);
void updateY(double *q, double *qs,double dt);
void update_density_X(double dt);
void update_density_Y(double dt);
void set_bc(void);
void ymax_bound(void);
void ymin_bound(void);
void read_param_file(char *filename);
void read_domain(char *directory);
void read_single_file(int n,int i, char *directory);
void read_files(int n, char *directory);
void output(char *directory);
void output_init(char *directory);


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
    dt = 0; //dt = .0000314159265359;

    allocate_all();
    read_domain(directory);
    read_single_file(n,0,directory);
    //set_bc();
    printf("%lg\n",omf);
    output_init(directory);
    
        printf("potential\n"); 
        potential();
        printf("viscl\n"); 
        viscosity();
        printf("sourcel\n"); 
        source_step();
        transport_step();
        set_Lamdep();
        set_wave();
    output(directory);
    free_all();
    return 0;
}
double Nu(double x) {
    return params.alpha*params.h*params.h*pow(x,2*params.flaringindex+.5);
}
double Cs(double x) {
    return params.h*pow(x,params.flaringindex-0.5);
}
void set_wave(void) {
    int i,j,k;
    k=0;
    for(j=NGHY;j<size_y-NGHY;j++) {
        drFw[j] = drFt[j]- drFd[j];
    }
    return;
}
void allocate_all(void) {
    
    MALLOC_SAFE((Ymin = (double *)malloc(sizeof(double)*(size_y+1))));
    MALLOC_SAFE((Xmin = (double *)malloc(sizeof(double)*(size_x+1))));

    MALLOC_SAFE((Ymed = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((Xmed = (double *)malloc(sizeof(double)*(size_x))));
    
    MALLOC_SAFE((dens = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((vx = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((vy = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((vx_temp = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((vy_temp = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Lang = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Pres = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((energy = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Pot = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Pixm = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Pixp = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Piym = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Piyp = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((denstar = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Qs = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((slope = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((divrho = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((tauxx = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((tauxy = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((tauyy = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));

    MALLOC_SAFE((tauxyavg = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((Lt = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((Lw = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((Ld = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((drFt = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((drFw = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((drFd = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((Lamex = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((Lamdep = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((dbarstar = (double *)malloc(sizeof(double)*(size_y))));

    int i,j,k;
    i = j = k = 0;

    for(j=0;j<size_y;j++) {
        Ymed[j] = 0;
        Ymin[j] = 0;
    }
    Ymin[size_y] = 0;
    for(i=0;i<size_x;i++) {
        Xmed[i] = 0;
        Xmin[i] = 0;
    }
    Xmin[size_x] = 0;
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
        
                
                dens[l]=0;
                vx[l]=0;
                vy[l]=0;
                vx_temp[l]=0;
                vy_temp[l]=0;
                Lang[l]=0;
                Pres[l]=0;
                energy[l]=0;
                Pot[l]=0;
                Pixm[l]=0;
                Pixp[l]=0;
                Piym[l]=0;
                Piyp[l]=0;
                denstar[l]=0;
                Qs[l]=0;
                slope[l]=0;
                divrho[l]=0;
                tauxx[l]=0;
                tauxy[l]=0;
                tauyy[l]=0;
            }
        }
    }

    return;
}
void free_all(void) {
    printf("Freeing.\n");
    
    FREE_SAFE(Ymin);
    FREE_SAFE(Xmin);

    FREE_SAFE(Ymed);
    FREE_SAFE(Xmed);
    
    FREE_SAFE(dens);
    FREE_SAFE(vx);
    FREE_SAFE(vy);
    FREE_SAFE(vx_temp);
    FREE_SAFE(vy_temp);
    FREE_SAFE(Lang);
    FREE_SAFE(Pres);
    FREE_SAFE(energy);
    FREE_SAFE(Pot);
    FREE_SAFE(Pixm);
    FREE_SAFE(Pixp);
    FREE_SAFE(Piym);
    FREE_SAFE(Piyp);
    FREE_SAFE(denstar);
    FREE_SAFE(Qs);
    FREE_SAFE(slope);
    FREE_SAFE(divrho);
    FREE_SAFE(tauxx);
    FREE_SAFE(tauxy);
    FREE_SAFE(tauyy);
    FREE_SAFE(tauxyavg);
    FREE_SAFE(Lt);
    FREE_SAFE(Lw);
    FREE_SAFE(Ld);
    FREE_SAFE(drFt);
    FREE_SAFE(drFw);
    FREE_SAFE(drFd);
    FREE_SAFE(Lamex);
    FREE_SAFE(Lamdep);
    FREE_SAFE(dbarstar);
    return;
}
void viscosity(void) {
    int i,j,k;
    i=j=k=0;
    double visc, viscm,div_v;
    compute_energy();
    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {

        visc = params.alpha*energy[l]*energy[l]*sqrt(ymed(j)*ymed(j)*ymed(j));
        viscm = params.alpha*.5*(energy[l]*energy[l]+energy[lym]*energy[lym])*sqrt(ymin(j)*ymin(j)*ymin(j));
            div_v = 0.0;
            div_v += (vx[lxp]-vx[l])*SurfX(j,k);
            div_v += (vy[lyp]*SurfY(j+1,k)-vy[l]*SurfY(j,k));
            div_v *= 2.0/3.0*InvVol(j,k);

            tauxx[l] = visc*dens[l]*(2.0*(vx[lxp]-vx[l])/zone_size_x(j,k) - div_v);
            tauxx[l] += visc*dens[l]*(vy[lyp]+vy[l])/ymed(j);
            tauyy[l] = visc*dens[l]*(2.0*(vy[lyp]-vy[l])/(ymin(j+1)-ymin(j)) - div_v);
            tauxy[l] = viscm*.25*(dens[l]+dens[lxm]+dens[lym]+dens[lxm-pitch])*((vy[l]-vy[lxm])/(dx*ymin(j)) + (vx[l]-vx[lym])/(ymed(j)-ymed(j-1))-.5*(vx[l]+vx[lym])/ymin(j)); //centered on left, inner vertical edge in z


        }
            tauxyavg[j] = viscm*.5*(dens[j]+dens[j-1])*((vxbar[j]-vxbar[j-1])/(ymed(j)-ymed(j-1))-.5*(vxbar[j]+vxbar[j-1])/ymin(j)); //centered on left, inner vertical edge in z
    }

    return;
}
void compute_Pres(void) {
    int i,j,k;
    i=j=k=0;
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
void compute_energy(void) {
    int i,j,k;
    i=j=k=0;
    double cs2;
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            energy[l] = params.h * pow(ymed(j),params.flaringindex) * sqrt(1./ymed(j));
        }
    }
    return;
}
void potential(void) {
    int i,j,k;
    i=j=k=0;
    double xpl, ypl;
    double smoothing,distp,rad;
    xpl = params.a;
    ypl = 0;
    distp = sqrt(xpl*xpl + ypl*ypl);

	smoothing = params.h*pow(distp,params.flaringindex)*distp*params.soft;
    smoothing *= smoothing;
    printf("Smoothing = %.16f\n",smoothing);
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            rad = (XC-xpl)*(XC-xpl) + (YC-ypl)*(YC-ypl);
            Pot[l] = -1./ymed(j);
            Pot[l] -= params.mp/sqrt(rad + smoothing);
	        //Pot[l] += params.mp*(xd*xpl+yd*ypl)/(distp*distp*distp);
        }
    }
    return;
}
void source_step(void) {
    int i,j,k;
    i=j=k=0;
    compute_Pres();
    double resfm, resfp, reslamp, reslamm;
    for(j=NGHY;j<size_y-NGHY;j++) {
        reslamm = 0;
        reslamp = 0;
        resfm = 0;
        resfp = 0;
        for(i=0;i<size_x;i++) {
            // X
            reslamm -= .5*(dens[l]+dens[lxm])*(Pot[l]-Pot[lxm])/zone_size_x(j,k);
            reslamp -= .5*(dens[l]+dens[lxp])*(Pot[lxp]-Pot[l])/zone_size_x(j,k);
	        resfm -= (ymin(j+1)*ymin(j+1)*tauxy[lyp]-ymin(j)*ymin(j)*tauxy[l])/((ymin(j+1)-ymin(j))*ymed(j));
	        resfp -= (ymin(j+1)*ymin(j+1)*tauxy[lxp+pitch]-ymin(j)*ymin(j)*tauxy[lxp])/((ymin(j+1)-ymin(j))*ymed(j));
       }
        reslamm /= (double)nx;
        reslamp /= (double)nx;
        resfm /= (double)nx;
        resfp /= (double)nx;
        drFt[j] = .5*(resfm + resfp);
        Lamex[j] = .5*(reslamm + reslamp);
        drFd[j] = -(ymin(j+1)*ymin(j+1)*tauxyavg[j+1]-ymin(j)*ymin(j)*tauxyavg[j])/((ymin(j+1)-ymin(j))*ymed(j));
    }


    return;
}
void transport_step(void) {
    set_momenta();    
    transportY();
    return;
}

void set_momenta(void) {
    int i,j,k;
    i=j=k=0;
    double res;
    for(j=0;j<size_y-1;j++) {
        res = 0;
        resw = 0;
        Ld[j] = dbar[j]*ymed(j)*(vxbar[j] + omf*ymed(j));
        for(i=0;i<size_x;i++) {
            Pixm[l] = ymed(j)*(vx[l]+omf*ymed(j))*dens[l];
            Pixp[l] = ymed(j)*(vx[lxp]+omf*ymed(j))*dens[l];
            res += .5*(Pixm[l] + Pixp[l]);
            resw += res - Ld[j];
        }
        res /= (double)nx;
        Lt[j] = res;
        Lw[j] = Lt[j] - Ld[j];
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
    updateY(Pixm,Qs,dt);

    DividebyRho(Pixp);
    vanleer_y_a(divrho);
    vanleer_y_b(divrho,Qs,dt);
    updateY(Pixp,Qs,dt);


    vanleer_y_a_avg(dbar);
    vanleer_y_b_avg(dbar,dbarstar,dt);
   
    DividebyRhoavg(Ld);
    vanleer_y_a_avg(divrho);
    vanleer_y_b_avg(divrho,Qs,dt);
    updateYavg(Ld,Qs,dt);

 
//    update_density_Y(dt);

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
    int i,j,k;
    i=j=k=0;
    for(j=0;j<size_y;j++) {
        divrho[j] = q[j]/dbar[j];
    }
    return;
}

void vanleer_y_a(double *q) {
    int i,j,k;
    i=j=k=0;
    double dqm, dqp;
    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {

            dqm = (q[l] - q[lym])/zone_size_y(j,k);
            dqp = (q[lyp]-q[l])/zone_size_y(j+1,k);

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
void vanleer_y_b(double *q, double *qs, double dt) {
    int i,j,k;
    i=j=k=0;

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
void vanleer_y_a_avg(double *q) {
    int i,j,k;
    i=j=k=0;
    double dqm, dqp;
    for(j=1;j<size_y-1;j++) {

            dqm = (q[j] - q[j-1])/zone_size_y(j,k);
            dqp = (q[j+1]-q[j])/zone_size_y(j+1,k);

            if (dqp*dqm <= 0) {
                slope[j] = 0;
            }
            else {
                slope[j] = 2*dqp*dqm/(dqp+dqm);
            }
    }
    
    return;
}
void vanleer_y_b_avg(double *q, double *qs, double dt) {
    int i,j,k;
    i=j=k=0;

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
void updateYavg(double *q, double *qs,double dt) {
    int i,j,k;

    for(j=NGHY;j<size_y-NGHY;j++) {
        drFd[j] -= ((vybar[j]*qs[j]*dbarstar[j]*SurfY(j,k)-vybar[j+1]*qs[j+1]*dbarstar[j+1]*SurfY(j+1,k))*InvVol(j,k));
    }
    return;
}
void updateY(double *q, double *qs,double dt) {
    int i,j,k;

    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {
            resp += ((vy[l]*qs[l]*denstar[l]*SurfY(j,k)-vy[lyp]*qs[lyp]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));

        }
        resp /= (double) nx;
        drFt[j] -= resp*.5;
    }
    return;
}
void set_Lamdep(double dt) {
    int i,j,k;
    double fac1, fac2, fac3, fac4;
    for(j=NGHY;j<size_y-NGHY;j++) {
        fac1 = 0; // Mass flux term
        fac2 = 0; // shear term
        fac3 = 0; // visc term 1
        fac4 = 0; // visc term 2

        for(i=0;i<size_x;i++) {

            fac1 += (((vy[l]-vybar[j])*(denstar[l]-dbarstar[j])*SurfY(j,k)
                        -(vy[lyp]-vybar[j+1])*(denstar[lyp]-dbarstar[j+1])*SurfY(j+1,k))*InvVol(j,k));

            fac2 += .5*(vy[l]-vybar[j])*(ymed(j)*(vx[l]-vxbar[j])-ymed(j-1)*(vx[lym]-vxbar[j+1]))/(ymed(j)-ymed(j-1));
            fac2 += .5*(vy[lyp]-vybar[j+1])*(ymed(j+1)*(vx[lyp]-vxbar[j+1])-ymed(j)*(vx[l]-vxbar[j]))/(ymed(j+1)-ymed(j));

        
            fac3 += (dens[l]-dbar[j])/dens[l];

            fac4 -= ymed(j)*2.0*(Pres[l]-Pres[lxm])/(zone_size_x(j,k)*(dens[l]+dens[lxm]));
            fac4 += ymed(j)*2.0*(tauxx[l]-tauxx[lxm])/(zone_size_x(j,k)*(dens[l]+dens[lxm]));
	        fac4 += 2.0*(ymin(j+1)*ymin(j+1)*(tauxy[lyp]-tauxyavg[j+1])-ymin(j)*ymin(j)*(tauxy[l]-tauxyavg[j]))/((ymin(j+1)-ymin(j))*ymed(j)*(dens[lxm]+dens[l]));

        }
        fac1 /= (double)nx;
        fac1 *= Ld[j]/dbar[j];
        
        fac2 /= (double)nx;
        fac2 *= -dbar[j];
        
        fac3 /= (double)nx;
        fac3 *= -(ymin(j+1)*ymin(j+1)*tauxyavg[j+1]-ymin(j)*ymin(j)*tauxyavg[j])/((ymin(j+1)-ymin(j))*ymed(j));



        fac4 /= (double)nx;
        fac4 *= dbar[j];

        Lamdep[j] = fac1 + fac2 + fac3 + fac4;
    }
    return;
}
void set_bc(void) {
    
    ymin_bound();
    ymax_bound();
    return;

}

void ymax_bound(void) {
    int i,j,k;

  int jact;
  int jgh;
  int kact;
  int kgh;
  int lgh;
  int lghs;
  int lactb;
  int lactbs;
  int lactbs_null;
  int nghy = NGHY;
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
  int kact;
  int kgh;
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
    double temp;
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
void read_single_file(int n, int i,char *directory) {
    char filename[512];
    FILE *f;
    sprintf(filename,"%ssubstep_%d_%d.dat",directory,i,n);
    printf("Reading %s\n",filename);
    f = fopen(filename,"r");
    if (f == NULL) printf("Error loading %s\n",filename); 
    fread(dens,sizeof(double),size_x*size_y*size_z,f);
    fread(vy,sizeof(double),size_x*size_y*size_z,f);
    fread(vx,sizeof(double),size_x*size_y*size_z,f);
    fclose(f);
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
void output_init(char *directory) {
    int i,j,k;
    i=j=k=0;
    FILE *f;
    char fname[256];
    sprintf(fname,"%stemp_files/output_init.dat",directory);
    printf("Outputting initial conditions to %s\n",fname);
    f = fopen(fname,"w");
    fwrite(&dens[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vy[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vx[0],sizeof(double),size_x*size_y*size_z,f);
    fclose(f);
    return;
}
void output(char *directory) {
    int i,j,k;
    i=j=k=0;
    FILE *f;
    char fname[256];
    sprintf(fname,"%stemp_files/output.dat",directory);
    printf("Outputing final results to %s\n",fname);
    f = fopen(fname,"w");
    fwrite(&Ymin[0],sizeof(double),size_y+1,f);
    fwrite(&Xmin[0],sizeof(double),size_x+1,f);
    fwrite(&dens[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vy[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vx[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&Lang[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&Pot[0],sizeof(double),size_x*size_y*size_z,f);

    fclose(f);
    return;
}
