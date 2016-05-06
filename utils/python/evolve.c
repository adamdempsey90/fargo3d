#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

#define NGHY 3
#define TRUE 1
#define FALSE 0
#define G 1.0
#define MSTAR 1.0
#define CFL 0.44
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define MALLOC_SAFE(ptr) if (ptr == NULL) printf("Malloc error at line %d!\n",__LINE__);
//#define FREE_SAFE(ptr) free(ptr); ptr=NULL; 

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

typedef struct Orbit {

    double a; // semi-major axis
    double e; // eccentricity
    double M; // mean anomaly
    double V; // true anomaly
    double psi; // arg of periastron from ascending ndoe
    double phi; // angle b/w actual and initial position of x-axis
    double i; // inclination
    double w; // longitude of ascending node w.r.t actual x-axis
    double alpha; // projection of perihelion w.r.t actual x-axis

} Orbit;

typedef struct Planet {

    double x,y,z;
    double vx,vy,vz;
    double mp, omf;
    double dist;
    double t;
    Orbit orbit;
} Planet;


double *dens, *vx, *vy, *Pres, *indPot,*Pot, *energy;
double *dbar, *vxbar, *vybar, *dbarstar;
double *vx_temp, *vy_temp;
double *Pixp, *Pixm, *Piym, *Piyp;
double *slope, *divrho, *denstar, *Qs;
double *Ymed, *Xmed, *Ymin, *Xmin;
double *tauxx, *tauxy, *tauyy, *tauxyavg;
double *Lt, *Ld, *Lw, *drFw, *drFd, *drFt, *Lamdep, *Lamex;
double *dtLt, *dtLd, *dtLw;

double dt,omf,dx;
int nx, ny, nz, size_x, size_y, size_z,stride,pitch;

Parameters params;
Planet planet;

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
void DividebyRhoavg(double *q);
void transportX(void);
void set_vel(void);
void vanleer_y_a(double *q);
void vanleer_y_a_avg(double *q);
void vanleer_x_a(double *q);
void vanleer_y_b(double *q, double *qs, double dt);
void vanleer_y_b_avg(double *q, double *qs, double dt);
void vanleer_x_b(double *q, double *qs, double dt);
void updateX(double *q, double *qs,double dt);
void updateY(double *q, double *qs,double dt);
void update_flux(double *qs);
void update_flux_avg(double *qs);
void update_density_X(double dt);
void update_density_Y(double dt);
void set_bc(void);
void ymax_bound(void);
void ymin_bound(void);
void ymin_bound_acc(void);
void ymax_bound_acc(void);
void read_param_file(char *filename);
void read_domain(char *directory);
void read_single_file(int n,int i, char *directory);
void read_files(int n, char *directory);
void output(char *directory);
void output_init(char *directory);
void set_Lamdep(void);
void set_waves(void);
void set_Lamex(void);
void output_torque(char *directory);
void set_avg(int p);
double cfl(void);
void read_planet_file(int n, char *directory);
void get_accel(double *q, double *k, double dtn);
void move_planet_step(double dtn);
void rotate_sys(double angle);
void move_planet(void);


int main(int argc, char *argv[]) {

    int n,nsteps;
    int status;
    double cfl_dt;
    char directory[256],outputdir[256];
    char param_fname[100];

    if (argc > 4) {
        printf("Too many arguments.Only using first 3.\n");
    }
    n = atoi(argv[1]);
    nsteps = atoi(argv[2]);
    strcpy(directory,argv[3]);

    sprintf(outputdir,"%stemp_files/",directory);
    status = mkdir(outputdir,S_IRWXU | S_IRWXG | S_IRWXO);
    if (status == -33)  {
        printf("Status is -33 for some reason.\n");
    }
    sprintf(param_fname,"%sparam_file.txt",directory);


    read_param_file(param_fname);

    size_x = nx;
    size_y = ny+2*NGHY;
    size_z = nz;
    stride = size_x*size_y;
    pitch = size_x;
    dx = 2*M_PI/nx;

    allocate_all();
    read_domain(directory);
    read_files(n,directory);
    //read_single_file(n,1,directory);
    cfl_dt = cfl();
    dt = cfl_dt;
    /*
    if (cfl_dt < dt) {
        printf("Using cfl limited timestep of %.4e instead of %.4e\n",cfl_dt,dt);
        dt = cfl_dt;
    }
    */
    set_bc();
    set_avg(0);
    output_init(outputdir);
    int i;
    for(i=0; i <nsteps; i++) {
    
        potential();
        move_planet();
    
       source_step();
       
        set_Lamex();
       viscosity();
    
        temp_to_vel();
        set_bc();
        
        vel_to_temp();
        
        transport_step();
      

    }
//    temp_to_vel();   
    set_avg(1);
    set_Lamdep();
    output(outputdir);
    output_torque(outputdir);
    //free_all();
    return 0;
}
double Nu(double x) {
    return params.alpha*params.h*params.h*pow(x,2*params.flaringindex+.5);
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
    MALLOC_SAFE((vx_temp = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((vy_temp = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Pres = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((energy = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((Pot = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
    MALLOC_SAFE((indPot = (double *)malloc(sizeof(double)*(size_y*size_x*size_z))));
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


    MALLOC_SAFE((Lt = (double *)malloc(sizeof(double)*(size_y*2))));
    MALLOC_SAFE((Ld = (double *)malloc(sizeof(double)*(size_y*2))));
    MALLOC_SAFE((Lw = (double *)malloc(sizeof(double)*(size_y*2))));
    MALLOC_SAFE((drFt = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((drFd = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((drFw = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((Lamdep = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((Lamex = (double *)malloc(sizeof(double)*(size_y*2))));
    MALLOC_SAFE((tauxyavg = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((dbar = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((vxbar = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((vybar = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((dbarstar = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((dtLt = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((dtLd = (double *)malloc(sizeof(double)*(size_y))));
    MALLOC_SAFE((dtLw = (double *)malloc(sizeof(double)*(size_y))));
    int i,j,k;
    i = j = k = 0;

    for(j=0;j<size_y;j++) {
        Ymed[j] = 0;
        Ymin[j] = 0;
        Lt[j] = 0;
        Ld[j] = 0;
        Lw[j] = 0;
        Lt[j+size_y] = 0;
        Ld[j+size_y] = 0;
        Lw[j+size_y] = 0;
        drFt[j] = 0;
        drFd[j] = 0;
        drFw[j] = 0;
        Lamex[j] = 0;
        Lamex[j+size_y] = 0;

        Lamdep[j] = 0;
        tauxyavg[j] = 0;
        vxbar[j] = 0;
        vybar[j] = 0;
        dbar[j] = 0;
        dbarstar[j] = 0;
        dtLt[j] = 0;
        dtLd[j] = 0;
        dtLw[j] = 0;
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
                Pres[l]=0;
                energy[l]=0;
                Pot[l]=0;
                indPot[l]=0;
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
    
    free(Ymin);
    free(Xmin);

    free(Ymed);
    free(Xmed);
    
    free(dens);
    free(vx);
    free(vy);
    free(vx_temp);
    free(vy_temp);
    free(Pres);
    free(energy);
    free(Pot);
    free(indPot);
    free(Pixm);
    free(Pixp);
    free(Piym);
    free(Piyp);
    free(denstar);
    free(Qs);
    free(slope);
    free(divrho);
    free(tauxx);
    free(tauxy);
    free(tauyy);

    free(Lt);
    free(Ld);
    free(Lw);
    free(drFt);
    free(drFd);
    free(drFw);
    free(Lamdep);
    free(Lamex);
    free(tauxyavg);
    free(dbar);
    free(vxbar);
    free(vybar);
    free(dbarstar);
    free(dtLt);
    free(dtLd);
    free(dtLw);
    return;
}
void viscosity(void) {
    int i,j,k;
    i=j=k=0;
    double visc, viscm,div_v;
    double res,fac,facp;
    visc = 0;
    viscm = 0;
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
        tauxyavg[j] = viscm*.5*(dbar[j]+dbar[j-1])*( (vxbar[j]-vxbar[j-1])/(ymed(j)-ymed(j-1))-.5*(vxbar[j]+vxbar[j-1])/ymin(j)); //centered on left, inner vertical edge in z
    }

    for(j=1;j<size_y-2;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
            // X

            vx_temp[l] += 2.0*(tauxx[l]-tauxx[lxm])/(zone_size_x(j,k)*(dens[l]+dens[lxm]))*dt;
            fac =  (ymin(j+1)*ymin(j+1)*tauxy[lyp]-ymin(j)*ymin(j)*tauxy[l])/((ymin(j+1)-ymin(j))*ymed(j));
            facp =  (ymin(j+1)*ymin(j+1)*tauxy[lxp+pitch]-ymin(j)*ymin(j)*tauxy[lxp])/((ymin(j+1)-ymin(j))*ymed(j));
            vx_temp[l] += dt*fac*2.0/(ymed(j)*(dens[l]+dens[lxm]));
            // Y
            vy_temp[l] += 2.0*(ymed(j)*tauyy[l]-ymed(j-1)*tauyy[lym])/((ymed(j)-ymed(j-1))*(dens[l]+dens[lym])*ymin(j))*dt;
            vy_temp[l] += 2.0*(tauxy[lxp]-tauxy[l])/(dx*ymin(j)*(dens[l]+dens[lym]))*dt;
            vy_temp[l] -= (tauxx[l]+tauxx[lym])/(ymin(j)*(dens[l]+dens[lym]))*dt;
            res += .5*(fac+facp);
        }
         drFt[j] = -res/(double)nx;
        //drFd[j] = -(ymin(j+1)*ymin(j+1)*tauxyavg[j+1]*SurfY(j+1,k) - ymin(j)*ymin(j)*tauxyavg[j]*SurfY(j,k))*InvVol(j,k);
        drFd[j]  =  -(ymin(j+1)*ymin(j+1)*tauxyavg[j+1]-ymin(j)*ymin(j)*tauxyavg[j])/((ymin(j+1)-ymin(j))*ymed(j));

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
    xpl = planet.x;
    ypl = planet.y;
    double mp = planet.mp;
    distp = sqrt(xpl*xpl + ypl*ypl);

    smoothing = params.h*pow(distp,params.flaringindex)*distp*params.soft;
    smoothing *= smoothing;

    double resx = 0;
    double resy = 0;
    double cellmass;
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++)  {
            cellmass = Vol(j,k)*dens[l];
            rad = XC*XC + YC*YC;
            rad = pow(rad+smoothing,-1.5);
            resx += G * cellmass * XC * rad;
            resy += G * cellmass * YC * rad;
        }
    }
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            rad = (XC-xpl)*(XC-xpl) + (YC-ypl)*(YC-ypl);
            Pot[l] = -G*MSTAR/ymed(j);
            Pot[l] += -G*mp/sqrt(rad + smoothing);
	        indPot[l] = G*planet.mp*(XC*xpl+YC*ypl)/(distp*distp*distp); 
	        indPot[l]  -= resx*XC + resy*YC ;
        }
    }
    return;
}
void source_step(void) {
    int i,j,k;
    i=j=k=0;
    double vxc;
    compute_Pres();
    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {
            // X
            vx_temp[l] = vx[l] -2*dt/(dens[l]+dens[lxm]) *(Pres[l]-Pres[lxm])/zone_size_x(j,k);
            vx_temp[l] -= dt*(Pot[l]-Pot[lxm])/zone_size_x(j,k);
            vx_temp[l] -= dt*(indPot[l]-indPot[lxm])/zone_size_x(j,k);

            // Y
            vy_temp[l] = vy[l] -2*dt/(dens[l]+dens[lym])*(Pres[l]-Pres[lym])/(ymed(j)-ymed(j-1));
            vxc = .25*(vx[l]+vx[lxp]+vx[lym]+vx[lxp-pitch])+omf*ymin(j);
            vy_temp[l] += dt*vxc*vxc/ymin(j);
            vy_temp[l] -= dt*(Pot[l]-Pot[lym])/(ymed(j)-ymed(j-1));
            vy_temp[l] -= dt*(indPot[l]-indPot[lym])/(ymed(j)-ymed(j-1));
        }
    }


    return;
}
void vel_to_temp(void) {
    int i,j,k;
    i = j = k =0;
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                vx_temp[l] = vx[l];
                vy_temp[l] = vy[l];
            }
        }
    }
    return;
}
void temp_to_vel(void) {
    int i,j,k;
    i = j = k = 0;
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                vx[l] = vx_temp[l];
                vy[l] = vy_temp[l];
            }
        }
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
    i=j=k=0;
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
    update_flux(Qs);
    updateY(Pixm,Qs,dt);

    DividebyRho(Pixp);
    vanleer_y_a(divrho);
    vanleer_y_b(divrho,Qs,dt);
    update_flux(Qs);
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


    vanleer_y_a_avg(dbar);
    vanleer_y_b_avg(dbar,dbarstar,dt);

    DividebyRhoavg(Ld);
    vanleer_y_a_avg(divrho);
    vanleer_y_b_avg(divrho,Qs,dt);
    update_flux_avg(Qs);



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
    int j;
    for(j=0;j<size_y;j++) {
        divrho[j] = q[j]/dbar[j];
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
    i=j=k=0;
    for(j=1;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            vy[l] = (Piym[l] + Piyp[lym])/(dens[l]+dens[lym]);
            vx[l] = (Pixm[l] + Pixp[lxm])/(ymed(j)*(dens[l]+dens[lxm])) - omf*ymed(j);
        }
    }
/*
    for(j=1;j<size_y;j++) {
        resv = 0;
        for(i=0;i<size_x;i++) {
            resv += .5*(Pixm[l] + Pixp[l]);
            //resv += dens[l]*(.5*(vx[lxp]+vx[l])+omf*ymed(j))*ymed(j);
        }
        resv /= (double)nx;
        dtLt[j] = -resv;
        dtLd[j] = -dbar[j]*(vx[j] + omf*ymed(j))*ymed(j);
    }
*/
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
void vanleer_y_a_avg(double *q) {
    int j;
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
void vanleer_x_a(double *q) {
    int i,j,k;
    i=j=k=0;
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
    i=j=k=0;

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
void vanleer_y_b_avg(double *q, double *qs, double dt) {
    int i,j,k;
    i=j=k=0;
    double res;
    for(j=1;j<size_y-1;j++) {
        res =0 ;
        for(i=0;i<size_x;i++) {
            res += vy_temp[l];
        }
        res /=(double)nx;
        if (res > 0.) {
            qs[j] = q[j-1] + .5*(zone_size_y(j-1,k)-res*dt)*slope[j-1];
        }
        else {
            qs[j] = q[j] - .5*(zone_size_y(j,k)+res*dt)*slope[j];
        }

    }
    return;
}
void vanleer_x_b(double *q, double *qs, double dt) {
    int i,j,k;
    i=j=k=0;

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
    i=j=k=0;
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {

            q[l] += ((vx_temp[l]*qs[l]*denstar[l]-vx_temp[lxp]*qs[lxp]*denstar[lxp])*SurfX(j,k)*dt*InvVol(j,k));

        }
    }

    return;
}
void update_flux_avg(double *qs) {
    int i,j,k;
    i=j=k=0;
    double res;
    for(j=0;j<size_y-1;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {    
            res += ((vy_temp[l]*qs[j]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*qs[j+1]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));
        }
        res /=(double)nx;

        drFd[j] -= res;
    }
    return;
}
void update_flux(double *qs) {
    int i,j,k;
    i=j=k=0;
    double res;
    for(j=0;j<size_y-1;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {    
            res += ((vy_temp[l]*qs[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*qs[lyp]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));
        }
        res /=(double)nx;
        drFt[j] -= res*.5;
    }
    return;
}
void updateY(double *q, double *qs,double dt) {
    int i,j,k;
    i=j=k=0;
    for(j=0;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {
            
            q[l] += dt*((vy_temp[l]*qs[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*qs[lyp]*denstar[lyp]*SurfY(j+1,k))*InvVol(j,k));
        }
    }
    return;
}
void update_density_Y(double dt) {
    int i,j,k;
    i=j=k=0;

    for(j=0;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {

            dens[l] += ((vy_temp[l]*denstar[l]*SurfY(j,k)-vy_temp[lyp]*denstar[lyp]*SurfY(j+1,k))*dt*InvVol(j,k));

        }
    }
    return;
}
void update_density_X(double dt) {
    int i,j,k;
    i = j= k =0;

    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {

            dens[l] += ((vx_temp[l]*denstar[l]-vx_temp[lxp]*denstar[lxp])*SurfX(j,k)*dt*InvVol(j,k));

        }
    }
    return;
}
void set_Lamdep(void) {
    int j;
    for(j=NGHY;j<size_y - NGHY;j++) {
        dtLt[j] = (Lt[j+size_y]-Lt[j])/dt;
        dtLd[j] = (Ld[j+size_y]-Ld[j])/dt;
        dtLw[j] = dtLt[j] - dtLd[j];
        Lamdep[j] = dtLd[j] + drFd[j];
        drFw[j] = drFt[j] - drFd[j];

    }

    return;
}
void set_Lamex(void) {
    int i,j,k;
    i=j=k=0;
    
    double res,resi;
    for(j=NGHY;j<size_y-NGHY;j++) {
        res = 0;
        resi = 0;
        for(i=0;i<size_x;i++) {
            res -= dens[l]*(Pot[lxp]-Pot[lxm])/(2*dx);
            resi -= dens[l]*(indPot[lxp]-indPot[lxm])/(2*dx);
        }
        res /=(double)nx;
        resi /=(double)nx;
        Lamex[j] = res;
        Lamex[j + size_y] = resi;

    }

    return;

}
void set_bc(void) {
    
    ymin_bound_acc();
    ymax_bound_acc();
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
void ymin_bound_acc(void) {

  int i;
  int j;
  int k;
  int jact;
  int jgh;
  int lgh;
  int lghs;
  int lactb;
  int lactbs_null;
  int nghy = NGHY;
  double sig1;
 double vr1;
 double ri1;
 double rm1;
  double omegaframe = omf;
  double nu_index = 0.5 + 2*params.flaringindex;
  double vr_index = -0.5 + 2*params.flaringindex;


  i = j = k = 0;

  for(k=0; k<size_z; k++) {
    for(j=0; j<NGHY; j++) {
      for(i=0; i<size_x; i++) {

	lgh = l;
	lghs = l;
	lactb = i + (2*nghy-j-1)*pitch + k*stride;
	lactbs_null = i + nghy*pitch + k*stride;
	jgh = j;
	jact = (2*nghy-j-1);

	sig1 = dens[ i + (nghy)*pitch + k*stride];
	ri1 = ymed(nghy);
	rm1 = ymin(nghy);
	vr1 = vy[ i + (nghy)*pitch + k*stride];
	dens[lgh] = sig1*pow(ri1/ymed(jgh),nu_index);
	vx[lgh] = (vx[lactb]+ymed(jact)*omegaframe)*sqrt(ymed(jact)/ymed(jgh))-ymed(jgh)*omegaframe;
	vy[lghs] = vr1*pow(ymed(jgh)/ri1,vr_index);
	vy[lactbs_null] = vr1*pow(rm1/ri1,vr_index);
      }
    }
  }
  return;
}

void ymax_bound_acc(void) {
  int i;
  int j;
  int k;
  int jact;
  int jgh;
  int lgh;
  int lghs;
  int lactb;
  int lactbs_null;
  double sig1;
 double ri1;
 double fh1;
  int nghy = NGHY;
  double mdot = params.mdot;
  double omegaframe = omf;
  double nu_0 = params.alpha*params.h*params.h;
  double nu_index = 0.5 + 2*params.flaringindex;
  double vnorm = -1.5*params.alpha*params.h*params.h;
  double vr_index = -0.5 + 2*params.flaringindex;
  double pi = M_PI;

  i = j = k = 0;

  for(k=0; k<size_z; k++) {
    for(j=0; j<NGHY; j++) {
      for(i=0; i<size_x; i++) {

	lgh = i + (ny+nghy+j)*pitch + k*stride;
	lghs = i + (ny+nghy+1+j)*pitch + k*stride;
	lactb = i + (ny+nghy-1-j)*pitch + k*stride;
	lactbs_null = i + (ny+nghy)*pitch + k*stride;
	jgh = (ny+nghy+j);
	jact = (ny+nghy-1-j);

	sig1 = dens[ i + (ny+nghy-1)*pitch + k*stride];
	ri1 = ymed(ny+nghy-1);
	fh1 = 3*pi*nu_0*pow(ri1,nu_index+0.5)*sig1;
	dens[lgh] = (fh1+mdot*(sqrt(ymed(jgh))-sqrt(ri1)))/(3*pi*nu_0*pow(ymed(jgh),nu_index+0.5));
	vx[lgh] = (vx[lactb]+ymed(jact)*omegaframe)*sqrt(ymed(jact)/ymed(jgh))-ymed(jgh)*omegaframe;
	if (j<size_y-1)
		vy[lghs] = vnorm*pow(ymin(jgh),vr_index)*mdot*sqrt(ymin(jgh))/(fh1+mdot*(sqrt(ymin(jgh))-sqrt(ri1)));
	vy[lactbs_null] = vnorm*pow(ymin(jgh),vr_index)*mdot*sqrt(ymin(jgh))/(fh1+mdot*(sqrt(ymin(jgh))-sqrt(ri1)));
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
    //fscanf(f,"%lg\n",&params.omf);
    fscanf(f,"%lg\n",&params.h);
    fscanf(f,"%lg\n",&params.flaringindex);
    fscanf(f,"%lg\n",&params.mdot);
    fscanf(f,"%lg\n",&params.soft);
    fscanf(f,"%lg\n",&dt);

    params.nuindex = 2*params.flaringindex + 0.5;
    params.vrindex = params.nuindex - 1.0;
    omf = 1.0;
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
    fx = fopen(filename,"r");
    if (fx == NULL) printf("Error reading %s\n",filename);

    for(i=0;i<size_x+1;i++) {
        fscanf(fx,"%lg\n",&Xmin[i]);
    }
    fclose(fx);
    sprintf(filename2,"%sdomain_y.dat",directory);
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
void read_planet_file(int n, char *directory) {
    char filename[512],filename2[512];
    FILE *f,*f1;
    int i;
    double xpl,ypl,zpl,vxpl,vypl,vzpl,mpl,tpl,omfpl;
    double epl,apl,vpl,psipl,phipl,ipl,wpl,alphapl;
    int scanned = 0;

    sprintf(filename,"%splanet0.dat",directory);
    f = fopen(filename,"r");
    while ( (scanned = fscanf(f,"%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                &i,&xpl,&ypl,&zpl,&vxpl,&vypl,&vzpl,&mpl,&tpl,&omfpl)) != EOF) {

            if (i == n) {
                planet.x = xpl; 
                planet.y = ypl; 
                planet.z = zpl; 
                planet.vx = vxpl; 
                planet.vy = vypl; 
                planet.vz = vzpl; 
                planet.mp = mpl; 
                planet.omf = omfpl; 
                planet.t = tpl;
                break;
            }

    }


    fclose(f);
    sprintf(filename2,"%sorbit0.dat",directory);
    f1 = fopen(filename2,"r");

    while ( (scanned = fscanf(f1,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                &tpl,&epl,&apl,&mpl,&vpl,&psipl,&phipl,&ipl,&wpl,&alphapl)) != EOF) {

            if (fabs(tpl-planet.t)<1e-8) {
                planet.orbit.e = epl; 
                planet.orbit.a = apl; 
                planet.orbit.M = mpl; 
                planet.orbit.V = vpl; 
                planet.orbit.psi = psipl; 
                planet.orbit.phi = phipl; 
                planet.orbit.i = ipl; 
                planet.orbit.w = wpl; 
                planet.orbit.alpha = alphapl; 
                break;
            }

    }
    fclose(f1);
    omf = planet.omf;
    return;

}
void read_single_file(int n, int i,char *directory) {
    char filename[512];
    FILE *f;
    sprintf(filename,"%ssubstep_%d_%d.dat",directory,i,n);
    f = fopen(filename,"r");
    if (f == NULL) printf("Error loading %s\n",filename); 
    fread(dens,sizeof(double),size_x*size_y*size_z,f);
    fread(vy,sizeof(double),size_x*size_y*size_z,f);
    fread(vx,sizeof(double),size_x*size_y*size_z,f);
    fclose(f);
    read_planet_file(n,directory);
    return;

}

void read_files(int n, char *directory) {
    char filename[512];
    FILE *fd,*fx,*fy;
    int i,j,k;
    sprintf(filename,"%sgasdens%d.dat",directory,n);
    fd = fopen(filename,"r");
    if (fd == NULL) printf("Error loading %s\n",filename); 
    sprintf(filename,"%sgasvx%d.dat",directory,n);
    fx = fopen(filename,"r");
    if (fx == NULL) printf("Error loading %s\n",filename); 
    sprintf(filename,"%sgasvy%d.dat",directory,n);
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

    read_planet_file(n,directory);
    return;

}
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
    fwrite(&dbar[NGHY],sizeof(double),ny,f);
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

    fclose(f);
    return;
}
void set_avg(int p) {
    int i,j,k;
    i=j=k=0;
    double resx,resy,resd;
    double resL;
    for(j=0;j<size_y;j++) {
        resx = 0;
        resy = 0;
        resd = 0;
        resL = 0;
        for(i=0;i<size_x;i++) {
            resx += vx[l];
            resy += vy[l];
            resd += dens[l];
            resL += ymed(j)*dens[l]*(.5*(vx[l] + vx[lxp]) + omf*ymed(j));
        }
        vxbar[j] = resx/(double)nx;
        vybar[j] = resy/(double)nx;
        dbar[j] = resd/(double)nx;
        Lt[j + p*size_y] = resL/(double)nx;

        Ld[j + p*size_y] = ymed(j)*( vxbar[j] + omf*ymed(j))*dbar[j];
        Lw[j + p*size_y] = Lt[j + p*size_y] - Ld[j + p*size_y];
    }
    return;
}
double cfl(void) {
    int i,j,k;
    i=j=k=0;
    double soundspeed2,soundspeed,visc;
    double cfl1_a, cfl1_b, cfl1;
    double cfl7_a, cfl7_b, cfl7;
    double cfl2, cfl3;
    double res,fac;

    res = 1e99;
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {
            soundspeed2 = energy[l]*energy[l];
            visc = params.alpha*energy[l]*energy[l]*sqrt(ymed(j)*ymed(j)*ymed(j));
            soundspeed = sqrt(soundspeed2);
            cfl1_a = soundspeed/zone_size_x(j,k);
            cfl1_b = soundspeed/zone_size_y(j,k);
            cfl1 = fmax(cfl1_a,cfl1_b);

	        cfl2 =fmax(fabs(vx[l]),fabs(vx[lxp]))/zone_size_x(j,k);
	        cfl3 = fmax(fabs(vy[l]),fabs(vy[lyp]))/zone_size_y(j,k);

	        cfl7_a = 1.0/zone_size_x(j,k);	
        	cfl7_b = 1.0/zone_size_y(j,k);
	        cfl7 = 4.0*visc*pow(fmax(cfl7_a,cfl7_b),2);


	        fac = CFL/sqrt(cfl1*cfl1 + cfl2*cfl2 + 
			     cfl3*cfl3 + cfl7*cfl7);

            if (fac < res) {
                res = fac;
            }
        }
    }

    return res;

}
void get_accel(double *q, double *k, double dtn) {
    double rp,coeff;
    rp = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);
    coeff = -G*MSTAR*pow(rp,-3.0);
    k[0] = q[3] * dtn;
    k[1] = q[4] * dtn;
    k[2] = q[5] * dtn;
    k[3] = coeff * q[0] * dtn;
    k[4] = coeff * q[1] * dtn;
    k[5] = coeff * q[2] * dtn;
    coeff = G*planet.mp/(rp*rp*rp);
    k[3] -= coeff * q[0] * dtn;
    k[4] -= coeff * q[1] * dtn;
    k[5] -= coeff * q[2] * dtn;
    return;
}
void move_planet_step(double dtn) {
    int i;
    double k1[6] ={0,0,0,0,0,0};
    double k2[6] ={0,0,0,0,0,0};
    double k3[6] ={0,0,0,0,0,0};
    double k4[6] ={0,0,0,0,0,0};
    double k5[6] ={0,0,0,0,0,0};
    double k6[6] ={0,0,0,0,0,0};
    double q0[6] ={0,0,0,0,0,0};
    double q1[6] ={0,0,0,0,0,0};
    double cvals1[5] =  {0.2, 0.0, 0.0, 0.0, 0.0};
    double cvals2[5] = {0.075, 0.225, 0.0, 0.0, 0.0};
    double cvals3[5] = { 0.3, -0.9, 1.2, 0.0, 0.0};
    double cvals4[5] = {-11.0/54.0, 2.5, -70.0/27.0, 35.0/27.0, 0.0 };
    double cvals5[5] = {1631.0/55296.0, 175.0/512.0, 575.0/13824.0, 44275.0/110592.0, 253.0/4096.0 };
    q0[0] = planet.x;
    q0[1] = planet.y;
    q0[2] = planet.z;
    q0[3] = planet.vx;
    q0[4] = planet.vy;
    q0[5] = planet.vz;

    get_accel(q0, k1, dtn);
    for(i=0;i<6;i++) {
        q1[i] = q0[i] + cvals1[0]*k1[i] + cvals1[1]*k2[i] + cvals1[2]*k3[i] + cvals1[3]*k4[i] + cvals1[4]*k5[i];
    }

    get_accel(q1, k2, dtn);
    for(i=0;i<6;i++) {
        q1[i] = q0[i] + cvals2[0]*k1[i] + cvals2[1]*k2[i] + cvals2[2]*k3[i] + cvals2[3]*k4[i] + cvals2[4]*k5[i];
    }


    get_accel(q1, k3, dtn);
    for(i=0;i<6;i++) {
        q1[i] = q0[i] + cvals3[0]*k1[i] + cvals3[1]*k2[i] + cvals3[2]*k3[i] + cvals3[3]*k4[i] + cvals3[4]*k5[i];
    }

    get_accel(q1, k4, dtn);
    for(i=0;i<6;i++) {
        q1[i] = q0[i] + cvals4[0]*k1[i] + cvals4[1]*k2[i] + cvals4[2]*k3[i] + cvals4[3]*k4[i] + cvals4[4]*k5[i];
    }

    get_accel(q1, k5, dtn);
    for(i=0;i<6;i++) {
        q1[i] = q0[i] + cvals5[0]*k1[i] + cvals5[1]*k2[i] + cvals5[2]*k3[i] + cvals5[3]*k4[i] + cvals5[4]*k5[i];
    }

    get_accel(q1, k6, dtn);

    for (i = 0; i < 6; i++) {
        q1[i] = (q0[i]+
                    37.0/378.0*k1[i]  + 
                    250.0/621.0*k3[i] + 
                    125.0/594.0*k4[i] + 
                    512.0/1771.0*k6[i]);
    }
    planet.x = q1[0];
    planet.y = q1[1];
    planet.z = q1[2];
    planet.vx = q1[3];
    planet.vy = q1[4];
    planet.vz = q1[5];
    return;
}
void rotate_sys(double angle) {

    double xt, yt,vxt,vyt;
    xt = planet.x;
    yt = planet.y;
    vxt = planet.vx;
    vyt = planet.vy;

    planet.x = xt * cos(angle)  + yt * sin(angle);
    planet.y = -xt * sin(angle)  + yt * cos(angle);

    planet.vx = vxt * cos(angle)  + vyt * sin(angle);
    planet.vy = -vxt * sin(angle)  + vyt * cos(angle);

    return;
}
void move_planet(void) {
    int i,j,k;
    int subcycles = 5;
    double dt_frac = 1./(double)subcycles;
    double omf_new;
    double xpl0,ypl0;
    double xpl1,ypl1;
    double d1, d2,cross,domega;

    xpl0 = planet.x;
    ypl0 = planet.y;

    for(i=0;i<subcycles;i++) {
        move_planet_step(dt*dt_frac);
    }
    xpl1 = planet.x;
    ypl1 = planet.y;
    
    d2 = sqrt(xpl1*xpl1 + ypl1*ypl1);
    d1 = sqrt(xpl0*xpl0+ypl0*ypl0);
    cross = xpl0*ypl1-xpl1*ypl0;
    omf_new = asin(cross/(d1*d2))/dt;
    domega = omf_new - omf;
    omf = omf_new;
    rotate_sys(omf*dt);

    i = j = k = 0;
    double res=0;
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            res = 0;
            for(i=0;i<size_x;i++) {
                vx[l] -= domega * ymed(j);
                res += vx[l];
            }
            res /=(double)nx;
            vxbar[j] = res;
        }
    }


    return;
}

