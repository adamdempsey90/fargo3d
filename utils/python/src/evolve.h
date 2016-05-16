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
#define MAXSTEPS 1000000
#define CFL 0.44
#define CVNR 1.41
#define CVNL 0.05
#define MINDT 1e-10
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define ARTIFICIALVISCOSITY


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
double *dbar,*dbart, *vxbar, *vybar, *dbarstar;
double *vx_temp, *vy_temp;
double *Pixp, *Pixm, *Piym, *Piyp;
double *slope, *divrho, *denstar, *Qs;
double *Ymed, *Xmed, *Ymin, *Xmin;
double *tauxx, *tauxy, *tauyy, *tauxyavg;
double *Lt, *Ld, *Lw, *drFw, *drFd, *drFt, *Lamdep, *Lamex;
double *dtLt, *dtLd, *dtLw;

double dt,omf,dx,time_step;
int nx, ny, nz, size_x, size_y, size_z,stride,pitch,nsteps;
int nb;
int IndirectTerm;

Parameters params;
Planet psys[2];

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
void time_avg(void);
void init_rk5(void);
void free_rk5(void);
void artificial_visc(void);
void move_to_com(void);
