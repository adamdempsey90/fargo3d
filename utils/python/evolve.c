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
double *vx_temp, *vy_temp;
double *Pixp, *Pixm, *Piym, *Piyp;
double *slope, *divrho, *denstar, *Qs;
double *Ymed, *Xmed, *Ymin, *Xmin;
double *tauxx, *tauxy, *tauyy;
double *Lang;

double dt,omf,dx;
int nx, ny, nz, size_x, size_y, size_z,stride,pitch;

Parameters params;

double Nu(double x);
double Cs(double X);


void set_ang(void);
void allocate_all(void);
void free_all(void);
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
void read_single_file(int n, char *directory);
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
    dt = M_PI*1e-5;

    allocate_all();
    read_domain(directory);
    read_single_file(n,directory);
    set_bc();
    printf("%lg\n",omf);
    output_init(directory);
    int i;
    for(i=0; i <nsteps; i++) {
        
        potential();
        source_step();
        viscosity();
       temp_to_vel();
        set_bc();
        vel_to_temp();
        transport_step();
        set_bc();
        set_ang();
    }
    

    
   
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
void set_ang(void) {
    int i,j,k;
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            Lang[l] = .5*(Pixp[l] + Pixm[l]);
        }
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
    return;
}
void free_all(void) {
    
    free(Ymin );
    free(Xmin );

    free(Ymed );
    free(Xmed );
    
    free(dens );
    free(vx );
    free(vy );
    free(vx_temp );
    free(vy_temp );
    free(Lang);

    free(Pres);
    free(Pot);
    free(Pixm );
    free(Pixp );
    free(Piym );
    free(Piyp );
    free(denstar );
    free(Qs );
    free(slope );
    free(divrho );
    free(tauxx );
    free(tauxy );
    free(tauyy );
    return;
}
void viscosity(void) {
    int i,j,k;
    k=0;
    double visc, viscm,div_v;
    for(j=1;j<size_y-1;j++) {
        visc = Nu(ymed(j));
        viscm = Nu(ymin(j));
        for(i=0;i<size_x;i++) {

            div_v = 0.0;
            div_v += (vx[lxp]-vx[l])*SurfX(j,k);
            div_v += (vy[lyp]*SurfY(j+1,k)-vy[l]*SurfY(j,k));
            div_v *= 2.0/3.0*InvVol(j,k);

            tauxx[l] = visc*dens[l]*(2.0*(vx[lxp]-vx[l])/zone_size_x(j,k) - div_v);
            tauxx[l] += visc*dens[l]*(vy[lyp]+vy[l])/ymed(j);
            tauyy[l] = visc*dens[l]*(2.0*(vy[lyp]-vy[l])/(ymin(j+1)-ymin(j)) - div_v);
            tauxy[l] = viscm*.25*(dens[l]+dens[lxm]+dens[lym]+dens[lxm-pitch])*((vy[l]-vy[lxm])/(dx*ymin(j)) + (vx[l]-vx[lym])/(ymed(j)-ymed(j-1))-.5*(vx[l]+vx[lym])/ymin(j)); //centered on left, inner vertical edge in z


        }
    }

    for(j=1;j<size_y-2;j++) {

        for(i=0;i<size_x;i++) {
            // X
  
	        vx_temp[l] += 2.0*(tauxx[l]-tauxx[lxm])/(zone_size_x(j,k)*(dens[l]+dens[lxm]))*dt;
	        vx_temp[l] += 2.0*(ymin(j+1)*ymin(j+1)*tauxy[lyp]-ymin(j)*ymin(j)*tauxy[l])/((ymin(j+1)-ymin(j))*ymed(j)*ymed(j)*(dens[lxm]+dens[l]))*dt;
            // Y
            vy_temp[l] += 2.0*(ymed(j)*tauyy[l]-ymed(j-1)*tauyy[lym])/((ymed(j)-ymed(j-1))*(dens[l]+dens[lym])*ymin(j))*dt;
            vy_temp[l] += 2.0*(tauxy[lxp]-tauxy[l])/(dx*ymin(j)*(dens[l]+dens[lym]))*dt;
            vy_temp[l] -= (tauxx[l]+tauxx[lym])/(ymin(j)*(dens[l]+dens[lym]))*dt;
       }
    }
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
void potential(void) {
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
void source_step(void) {
    int i,j,k;
    k=0;
    double vxc;
    compute_Pres();

    for(j=1;j<size_y-1;j++) {

        for(i=0;i<size_x;i++) {
            // X
            vx_temp[l] =vx[l] -2*dt/(dens[l]+dens[lxm]) *(Pres[l]-Pres[lxm])/zone_size_x(j,k);
            vx_temp[l] -= dt*(Pot[l]-Pot[lxm])/zone_size_x(j,k);
  
            // Y
            vy_temp[l] = vy[l] -2*dt/(dens[l]+dens[lym])*(Pres[l]-Pres[lym])/(ymed(j)-ymed(j-1));
            vxc = .25*(vx[l]+vx[lxp]+vx[lym]+vx[lxp-pitch])+omf*ymin(j);
            vy_temp[l] += dt*vxc*vxc/ymin(j);
            vy_temp[l] -= dt*(Pot[l]-Pot[lym])/(ymed(j)-ymed(j-1));
       }
    }


    return;
}
void vel_to_temp(void) {
    int i;
    for(i=0;i<size_x*size_y*size_z;i++) {
        vy_temp[i] = vy[i];
        vx_temp[i] = vx[i];
    }
    return;
}
void temp_to_vel(void) {
    int i;
    for(i=0;i<size_x*size_y*size_z;i++) {
        vy[i] = vy_temp[i];
        vx[i] = vx_temp[i];
    }
    return;
}
void transport_step(void) {
    set_momenta();    
    transportY();
    transportX();
    set_vel();
    return;
}

void set_momenta(void) {
    int i,j,k;
    k=0;

    for(j=0;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {
            Pixm[l] = ymed(j)*(vx_temp[l]+omf*ymed(j))*dens[l];
            Pixp[l] = ymed(j)*(vx_temp[lxp]+omf*ymed(j))*dens[l];
            Piym[l] = vy_temp[l]*dens[l];
            Piyp[l] = vy_temp[lyp]*dens[l];
        }
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

    DividebyRho(Piym);
    vanleer_y_a(divrho);
    vanleer_y_b(divrho,Qs,dt);
    updateY(Piym,Qs,dt);

    DividebyRho(Piyp);
    vanleer_y_a(divrho);
    vanleer_y_b(divrho,Qs,dt);
    updateY(Piyp,Qs,dt);
 
    update_density_Y(dt);

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
void transportX(void) {
    // X direction

    vanleer_x_a(dens);
    vanleer_x_b(dens,denstar,dt);
  
    DividebyRho(Pixm);
    vanleer_x_a(divrho);
    vanleer_x_b(divrho,Qs,dt);
    updateX(Pixm,Qs,dt);

    DividebyRho(Pixp);
    vanleer_x_a(divrho);
    vanleer_x_b(divrho,Qs,dt);
    updateX(Pixp,Qs,dt);

    DividebyRho(Piym);
    vanleer_x_a(divrho);
    vanleer_x_b(divrho,Qs,dt);
    updateX(Piym,Qs,dt);

    DividebyRho(Piyp);
    vanleer_x_a(divrho);
    vanleer_x_b(divrho,Qs,dt);
    updateX(Piyp,Qs,dt);
  
    update_density_X(dt);

    return;

}
void set_vel(void) {

    int i,j,k;
    k=0;

    for(j=1;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            vy[l] = (Piym[l] + Piyp[lym])/(dens[l]+dens[lym]);
            vx[l] = (Pixm[l] + Pixp[lxm])/(ymed(j)*(dens[l]+dens[lxm])) - omf*ymed(j);
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
void vanleer_x_a(double *q) {
    int i,j,k;
    k=0;
    double dqm, dqp;
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {

            dqm = (q[l] - q[lxm]);
            dqp = (q[lxp]-q[l]);

            if (dqp*dqm <= 0) {
                slope[l] = 0;
            }
            else {
                slope[l] = 2*dqp*dqm/((dqp+dqm)*zone_size_x(j,k));
            }
        }
    }
    
    return;
}
void vanleer_y_b(double *q, double *qs, double dt) {
    int i,j,k;
    k=0;

    for(j=1;j<size_y-1;j++) {

        for(i=0;i<size_x;i++) {
            
            if (vy_temp[l] > 0.) {
                qs[l] = q[lym] + .5*(zone_size_y(j-1,k)-vy_temp[l]*dt)*slope[lym];
            }
            else {
                qs[l] = q[l] - .5*(zone_size_y(j,k)+vy_temp[l]*dt)*slope[l];
            }

        }
    }
    return;
}
void vanleer_x_b(double *q, double *qs, double dt) {
    int i,j,k;
    k=0;

    for(j=0;j<size_y;j++) {

        for(i=0;i<size_x;i++) {
            
            if (vx_temp[l] > 0.) {
                qs[l] = q[lxm] + .5*(zone_size_x(j,k)-vx_temp[l]*dt)*slope[lxm];
            }
            else {
                qs[l] = q[l] - .5*(zone_size_x(j,k)+vx_temp[l]*dt)*slope[l];
            }

        }
    }
    return;
}
void updateX(double *q, double *qs,double dt) {
    int i,j,k;

    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {

            q[l] += ((vx_temp[l]*qs[l]*denstar[l]-vx_temp[lxp]*qs[lxp]*denstar[lxp])*SurfX(j,k)*dt*InvVol(j,k));

        }
    }



    return;
}
void updateY(double *q, double *qs,double dt) {
    int i,j,k;

    for(j=0;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {

            q[l] += ((vy_temp[l]*qs[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*qs[lyp]*denstar[lyp]*SurfY(j+1,k))*dt*InvVol(j,k));

        }
    }
    return;
}
void update_density_Y(double dt) {
    int i,j,k;

    for(j=0;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {

            dens[l] += ((vy_temp[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*denstar[lyp]*SurfY(j+1,k))*dt*InvVol(j,k));

        }
    }
    return;
}
void update_density_X(double dt) {
    int i,j,k;

    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {

            dens[l] += ((vx_temp[l]*denstar[l]-vx_temp[lxp]*denstar[lxp])*SurfX(j,k)*dt*InvVol(j,k));

        }
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
void read_single_file(int n, char *directory) {
    char filename[512];
    FILE *f;
    sprintf(filename,"%ssubstep_0_%d.dat",directory,n);
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
    k=0;
    FILE *f;
    char fname[256];
    sprintf(fname,"%stemp_files/output_init.dat",directory);
    f = fopen(fname,"w");
    fwrite(&dens[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vy[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vx[0],sizeof(double),size_x*size_y*size_z,f);
    fclose(f);
    return;
}
void output(char *directory) {
    int i,j,k;
    k=0;
    FILE *f;
    char fname[256];
    sprintf(fname,"%stemp_files/output.dat",directory);
    f = fopen(fname,"w");
    fwrite(&Ymin[0],sizeof(double),size_y+1,f);
    fwrite(&Xmin[0],sizeof(double),size_x+1,f);
    fwrite(&dens[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vy[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&vx[0],sizeof(double),size_x*size_y*size_z,f);
    fwrite(&Lang[0],sizeof(double),size_x*size_y*size_z,f);
    fclose(f);
    return;
}
