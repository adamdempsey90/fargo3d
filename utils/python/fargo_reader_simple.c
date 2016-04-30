#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NGHY 3
#define SYMMETRIC
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


int nx, ny, nz;
int size_x, size_y, size_z;
int stride, pitch;
double dx;
double *Ymin, *Xmin, *Ymed, *Xmed;
double *rho, *vx, *vy;
double *pot;
double *sig, *vxp, *vyp;
double *tauxx, *tauxy, *tauyy;
double *rhobar, *vxbar, *vybar;
double *tauxyavg, *divpavg;

double *Ltotstar, *rhostar;
double *Lstar;
double *Mdottot, *Mdotd;
double *vrbar,*Ldstar;
double *shear;
double *work, *work1;

Parameters params;


double *Ltot, *Ld, *Lw;
double *drFtot, *drFd, *drFw;
double *lamex, *lamind, *lamdep;


double Nu(double x);
double Cs(double x);
void allocate_all(void);
void free_all(void);
void set_tensor(void);
void set_div(void);
void set_Ltot(void);
void set_Ld(void);
void set_Mdotd(void);
void set_Mdottot(void);
void set_shear(void);
void set_drFtot(void);
void set_drFd(void);
void set_waves(void);
void set_lam_ex(void);
void set_pot(void);
void set_lam_dep(void);
void set_lamind(void);
void add_boundary(void);
void read_param_file(char *filename);
void read_domain(char *directory);
void read_files(int n, char *directory);
void set_averages(void);
void set_slopesbar(double *q);
void transportbar(double *qs,double *vel,double *res);
void set_slopes(double *q);
void transport(double *qs,double *vel,double *res);
void clear_work(void);
void output(char *directory);




int main(int argc, char *argv[]) {

    int n;
    int i,j,k;
    char directory[256];
    char param_fname[100];

    n = atoi(argv[1]);
    strcpy(directory,argv[2]);
    sprintf(param_fname,"%sparam_file.txt",directory);

    printf("Reading param file %s\n",param_fname);
    read_param_file(param_fname);

    size_x = nx;
    size_y = ny+2*NGHY;
    size_z = nz;
    stride = size_x*size_y;
    pitch = size_x;
    dx = 2*M_PI/nx;

    printf("Allocating.\n");
    allocate_all();


    printf( "Reading domain.\n");
    read_domain(directory); 
    read_files(n,directory);
    printf( "Add bc.\n");
    add_boundary();
    printf( "Computing averages.\n");
    set_averages();

    printf( "Computing Mdotd.\n");
    set_Mdotd();
    printf( "Computing tensor.\n");
    set_tensor();
    printf( "Computing Ltot.\n");
    set_Ltot();
    printf( "Computing Ld.\n");
    set_Ld();
    printf( "Computing Mdottot.\n");
    set_Mdottot();
    printf( "Computing divpavg.\n");
    set_div();
    printf( "Computing shear.\n");
    set_shear();
    printf( "Computing drFtot.\n");
    set_drFtot();
    printf( "Computing drFd.\n");
    set_drFd();
    printf( "Computing waves.\n");
    set_waves();
    printf( "Computing lamex.\n");
    set_lam_ex();
    set_lam_dep();

    printf("Output.\n");
    output(directory);
    printf("Exiting.\n");
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
    MALLOC_SAFE((rho = (double *)malloc(sizeof(double)*size_x*size_y*size_z))); 
    MALLOC_SAFE((pot = (double *)malloc(sizeof(double)*size_x*size_y*size_z))); 
    MALLOC_SAFE((vx = (double *)malloc(sizeof(double)*size_x*size_y*size_z))); 
    MALLOC_SAFE((vy = (double *)malloc(sizeof(double)*size_x*size_y*size_z))); 
    MALLOC_SAFE((sig = (double *)malloc(sizeof(double)*size_x*size_y*size_z))); 
    MALLOC_SAFE((vxp = (double *)malloc(sizeof(double)*size_x*size_y*size_z))); 
    MALLOC_SAFE((vyp = (double *)malloc(sizeof(double)*size_x*size_y*size_z))); 
    MALLOC_SAFE((tauxx = (double *)malloc(sizeof(double)*size_x*size_y*size_z))); 
    MALLOC_SAFE((tauxy = (double *)malloc(sizeof(double)*size_x*size_y*size_z))); 
    MALLOC_SAFE((tauyy = (double *)malloc(sizeof(double)*size_x*size_y*size_z))); 
    MALLOC_SAFE((Ltotstar = (double *)malloc(sizeof(double)*size_x*size_y*size_z))); 
    MALLOC_SAFE((rhostar = (double *)malloc(sizeof(double)*size_x*size_y*size_z))); 
    MALLOC_SAFE((Lstar = (double *)malloc(sizeof(double)*size_x*size_y*size_z))); 

    MALLOC_SAFE((tauxyavg = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((divpavg = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((rhobar = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((vxbar = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((vybar = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((Ltot = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((Ld = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((Lw = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((drFd = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((drFw = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((drFtot = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((Mdottot = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((Mdotd = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((vrbar = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((Ldstar = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((shear = (double *)malloc(sizeof(double)*size_y*size_z))); 

    MALLOC_SAFE((lamex = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((lamdep = (double *)malloc(sizeof(double)*size_y*size_z))); 
    MALLOC_SAFE((lamind = (double *)malloc(sizeof(double)*size_y*size_z))); 

    MALLOC_SAFE((Xmed = (double *)malloc(sizeof(double)*size_x))); 
    MALLOC_SAFE((Xmin = (double *)malloc(sizeof(double)*size_x))); 
    MALLOC_SAFE((Ymed = (double *)malloc(sizeof(double)*size_y))); 
    MALLOC_SAFE((Ymin = (double *)malloc(sizeof(double)*(size_y+1)))); 

    MALLOC_SAFE((work = (double *)malloc(sizeof(double)*size_x*size_y*size_z))); 
    MALLOC_SAFE((work1 = (double *)malloc(sizeof(double)*size_x*size_y*size_z))); 

    int i,j,k;
    k=0;

    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            tauxyavg[j]= 0;
            rhobar[j]= 0;
            vxbar[j]= 0;
            vybar[j]= 0;
            Ltot[j]= 0;
            Ld[j]= 0;
            Lw[j]= 0;
            drFd[j]= 0;
            drFw[j]= 0;
            drFtot[j]= 0;
            Mdottot[j]= 0;
            Mdotd[j]= 0;
            vrbar[j]= 0;
            Ldstar[j]= 0;
            shear[j]= 0;

            lamex[j]= 0;
            lamdep[j]= 0;
            lamind[j]= 0;

            for(i=0;i<size_x;i++) {

                rho[l]= 0;
                pot[l]= 0;
                vx[l]= 0;
                vy[l]= 0;
                sig[l]= 0;
                vxp[l]= 0;
                vyp[l]= 0;
                tauxx[l]= 0;
                tauxy[l]= 0;
                tauyy[l]= 0;
                Ltotstar[l]= 0;
                rhostar[l]= 0;
                Lstar[l]= 0;
                work[l]= 0;
                work1[l]= 0;

            }
        }
    }

    return;
}
void free_all(void) {
    free(rho);
    free(pot);
    free(vx);
    free(vy);
    free(sig);
    free(vxp);
    free(vyp);
    free(tauxx);
    free(tauxy);
    free(tauyy);
    free(Ltotstar);
    free(rhostar);
    free(Lstar);

    free(tauxyavg);
    free(divpavg);
    free(rhobar);
    free(vxbar);
    free(vybar);
    free(Ltot);
    free(Ld);
    free(Lw);
    free(drFd);
    free(drFw);
    free(drFtot);
    free(Mdottot);
    free(Mdotd);
    free(vrbar);
    free(Ldstar);
    free(shear);

    free(lamex);
    free(lamdep);
    free(lamind);

    free(Xmed);
    free(Xmin);
    free(Ymed);
    free(Ymin);

    free(work);
    free(work1);

    return;
}

void set_tensor(void) {
    int i,j,k;
    double visc, viscm, div_v;
    double cs, csm;
    double rhoc, rhocb;
    k = 0;


    for(j=1;j<size_y-1;j++) {
        cs = Cs(ymed(j));
        csm = Cs(ymin(j));
        cs *= cs; csm *= csm;
        visc = Nu(ymed(j));
        viscm = Nu(ymin(j));
        tauxyavg[j] = viscm*.5*(rhobar[j]+rhobar[j-1])*( 
                (vxbar[j]-vxbar[j-1])/(ymed(j)-ymed(j-1))
                -.5*(vxbar[j]+vxbar[j-1])/ymin(j));

        for(i=0;i<size_x;i++) {
            rhoc=.25*(rho[l]+rho[lxm]+rho[lym]+rho[lxm-pitch]);
            div_v = (vx[lxp]-vx[l])*SurfX(j,k);
            div_v += (vy[lyp]*SurfY(j+1,k)-vy[l]*SurfY(j,k));
            div_v *= 2.0/3.0*InvVol(j,k);

            tauxx[l] = visc*rho[l]*(2.0*(vx[lxp]-vx[l])/zone_size_x(j,k) - div_v);
            tauxx[l] += visc*rho[l]*(vy[lyp]+vy[l])/ymed(j);
            tauxx[l] += -cs*rho[l]; // Isothermal Pressure
            tauyy[l] = visc*rho[l]*(2.0*(vy[lyp]-vy[l])/(ymin(j+1)-ymin(j)) - div_v);
            tauyy[l] += -cs*rho[l]; // Isothermal Pressure
            tauxy[l] = viscm*rhoc*((vy[l]-vy[lxm])/(dx*ymin(j)) + ymin(j)*(vx[l]/ymed(j)-vx[lym]/ymed(j-1))/(ymed(j)-ymed(j-1)));//-.5*(vx[l]+vx[lym])/ymin(j)); //centered on left, inner vertical edge in z
            

        }
    }
    return;
}

void set_Ltot(void) {
    int i,j,k;
    k=0;
    double res;

    for(j=0;j<size_y;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
            Ltotstar[l] = 2*M_PI*ymed(j)*ymed(j)*rho[l]*.5*(vx[l] + vx[lxp]);
            res += Ltotstar[l];
        }
        res /= (double) nx;
        Ltot[j] = res;
    }

    return;
}

void set_Ld(void) {
    int i,j,k;
    k=0;

    for(j=0;j<size_y;j++) {
        Ld[j] = 2*M_PI*ymed(j)*ymed(j)*rhobar[j]*vxbar[j];
    }

    return;
}

void set_Mdotd(void) {
    int i,j,k;
    k=0;
    clear_work();
    transport(rho,vy,rhostar);

    double res,res1;

    for(j=0;j<size_y;j++) {
        res = 0;
        res1 = 0;
        for(i=0;i<size_x;i++) {
            res +=  rhostar[l]*vy[l];
            res1 += rhostar[l];
        }
        res /= (double) nx;
        res1 /= (double)nx;
        Mdotd[j] = -2*M_PI*ymin(j)*res;
        vrbar[j] = -Mdotd[j]/(2*M_PI*ymin(j)*res1);
    }


    return;
}

void set_Mdottot(void) {
    int i,j,k;
    k=0;

    clear_work();
    transport(Ltotstar,vy,NULL);

    double res;

    for(j=NGHY;j<size_y-NGHY+1;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
            res +=  Ltotstar[l]*vy[l];
        }
        res /= (double) nx;
        Mdottot[j] = res;
    }
    return;
}
void set_drFtot(void) {
    int i,j,k;
    k=0;

    double res;
    for(j=NGHY;j<size_y-NGHY;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
	        res -= 2*M_PI*(ymin(j+1)*ymin(j+1)*tauxy[lyp]-ymin(j)*ymin(j)*tauxy[l])/((ymin(j+1)-ymin(j)));
        }
        res /= (double) nx;
        res += (Mdottot[j+1]-Mdottot[j])/(ymin(j+1)-ymin(j));
        drFtot[j] = res;
    }

    return;
}
void set_drFd(void) {
    int i,j,k;
    k=0;
    
    clear_work();
    transportbar(Ld,vrbar,Ldstar);

    for(j=NGHY;j<size_y-NGHY;j++) {
        drFd[j] = -(Mdotd[j+1]*Ldstar[j+1]-Mdotd[j]*Ldstar[j])/(ymin(j+1)-ymin(j));
        drFd[j] -= 2*M_PI*(ymin(j+1)*ymin(j+1)*tauxyavg[j+1]-ymin(j)*ymin(j)*tauxyavg[j])/(ymin(j+1)-ymin(j));
    }

    return;
}
void set_waves(void) {
    int j;


    for(j=NGHY;j<size_y;j++) {
        Lw[j] = Ltot[j] - Ld[j];
        drFw[j] = drFtot[j] - drFd[j];
    }
    return;
}
void set_pot(void) {
    int i,j,k;
    k = 0;

    double smoothing;// = params.soft*params.h*pow(params.a,1 + params.flaringindex);
    //smoothing *= smoothing;
    double rad,dist;
    double potp,potm;
    double xpl, ypl, distp;
    double xd, yd;
    xpl = params.a;
    ypl = 0;
    distp = sqrt(xpl*xpl + ypl*ypl);
    
    smoothing = params.h*pow(distp,params.flaringindex)*params.soft*distp;
    smoothing *= smoothing;
    printf("%lg\n",params.mp);
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {    
            xd = ymed(j)*cos(xmin(i));
            yd = ymed(j)*sin(xmin(i));
            //pot[l] = -1/ymed(j);
            rad = (xd-xpl)*(xd-xpl) + (yd-ypl)*(yd-ypl);
            pot[l] = -params.mp / sqrt(rad + smoothing);
            //pot[l] += params.mp*(yd*ypl + xd*xpl)*pow(distp,-3.0); 
        }
    }
    return;

}

void set_lam_ex(void) {
    int i,j,k;
    k=0;
    double res;
    set_pot();
    for(j=NGHY;j<size_y-NGHY;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
            res += rho[l]*(pot[lxp]-pot[l])/(dx);
        }
        res /=(double)nx;
        res *= -2*M_PI*ymed(j);
        lamex[j] = res;
    }

    return;

}
/*
void set_lam_ex(void) {

    int i,j,k;
    k = 0;

    double smoothing = params.soft*params.h*pow(params.a,1 + params.flaringindex);
    smoothing *= smoothing;
    double radp,radm,res;
    double potp,potm;
    for(j=NGHY;j<size_y-NGHY;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {    
            radm = ymed(j)*ymed(j) + params.a*params.a -2*params.a*ymed(j)*cos(xmin(ixp))+smoothing;
            radp = ymed(j)*ymed(j) + params.a*params.a -2*params.a*ymed(j)*cos(xmin(i))+smoothing;
            potp = pow(radp,-1.5);
            potm = pow(radm,-1.5);
            res += 2*M_PI*ymed(j)*rho[l]*params.mp*(potp-potm)/dx;
            //pot = -params.mp * ymed(j)*params.a*sin(xmed(i))*pow(rad,-1.5);
            //res += -2*M_PI*ymed(j)*rho[l]*pot;
        }
        res /= (double)nx;
        lamex[j] = res;
    }

    return;
}
*/
void set_div(void) {
    int i,j,k;
    k=0;
    double res,res1;
    for(j=NGHY;j<size_y-NGHY;j++) {
        res = 0;
        res1=0;
       for(i=0;i<size_x;i++) {
            res1 += (rho[l]-rhobar[j])/rho[l];
	        res += 2*M_PI*ymed(j)*2.0*(tauxx[l]-tauxx[lxm])/(dx*(rho[l]+rho[lxm]));
	        res += 2*M_PI*2.0*(ymin(j+1)*ymin(j+1)*(tauxy[lyp]-tauxyavg[j+1])-ymin(j)*ymin(j)*(tauxy[l]-tauxyavg[j]))/((ymin(j+1)-ymin(j))*(rho[lxm]+rho[l]));
        }
       res /= (double)nx;
       res1 /= (double)nx;
       res *= rhobar[j]; //*2*M_PI*ymed(j)*ymed(j);
       res -= res1*2*M_PI*(ymin(j+1)*ymin(j+1)*tauxyavg[j+1]-ymin(j)*ymin(j)*tauxyavg[j])/(ymin(j+1)-ymin(j));
       divpavg[j] = res;
    }


    return;


}

void set_shear(void) {
    int i,j,k;
    k=0;
    double res;
    double dy, res1,res2,fac,facm,facmp,lspp, lsp, lp,lmed;
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            Lstar[l] = Ltotstar[l]/(2*M_PI*ymin(j)*rhostar[l]);
        }
    }
    for(j=NGHY;j<size_y-NGHY;j++) {
        res1 = 0;
        res2=0;
        for(i=0;i<size_x;i++) {
            res2 += Lstar[lyp];
            res1 += Lstar[l];
        }
        res1 /= (double)nx;
        res2 /=(double)nx;
        res = 0;

        dy = ymin(j+1)-ymin(j);
        fac = 2*M_PI*ymed(j);
        facm = 2*M_PI*ymin(j);
        facmp = 2*M_PI*ymin(j+1);
        lmed = Ld[j]/(fac*rhobar[j]);
        
        for(i=0;i<size_x;i++) {
            lspp = Lstar[lyp] - res2;
            lsp = Lstar[l] - res1;
            lp = Ltot[l]/(fac*rho[l]) - lmed;
            res += ( vyp[lyp] * facmp*sig[lyp] * lspp - vyp[l]*facm*sig[l]*lsp )/dy - lmed*(vyp[lyp] * facmp*sig[lyp] - vyp[l] * facm*sig[l])/dy;
            res += -fac*rhobar[j]*( (vyp[lyp]*lspp-vyp[l]*lsp)/dy - lp*(vyp[lyp]-vyp[l])/dy);

//res += ((vy[lyp]-vybar[j+1])*Lstar[lyp] - vy[l]*Lstar[l])/(ymin(j+1)-ymin(j));
//            res -= .5*(vx[l]+vx[lxp])*ymed(j)*(vy[lyp]-vy[l])/(ymin(j+1)-ymin(j));
        }
        res /= (double)nx;
        shear[j] = res;
    }

    return;
}
void set_lam_dep(void) {
    int i,j,k;
    k=0;
    double fac,res,sfac;

    for(j=NGHY;j<size_y-NGHY;j++) {
        /*
        fac = 2*M_PI*ymed(j);
        res=0;
        for(i=0;i<size_x;i++) {
          
            res -= (Mdotd[j+1]*Ldstar[j+1]-Mdotd[j]*Ldstar[j])/(ymin(j+1)-ymin(j));
            res += (Mdotd[j+1]-Mdotd[j])/(ymin(j+1)-ymin(j)) * Ld[j];
            

        }
        res /= (double)nx;
        res += -fac*rhobar[j]*shear[j]+divpavg[j];
    */
        lamdep[j] = shear[j] + divpavg[j];
    }

    return;
}
void set_lamind(void) {
    int i,j,k;
    k=0;

    double res,resj;

    double resx, resy;
    double dm;
    resx=0;
    resy=0;
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {
           
            dm= Vol(j,k)*rho[l];
            resx += dm*cos(xmed(i))/(ymed(j)*ymed(j));
            resy += dm*sin(xmed(i))/(ymed(j)*ymed(j));
            
        }

    }

    for(j=NGHY;j<size_y-NGHY;j++) {
        res = 0;
        for(i=0;i<size_x;i++) {
            
            resj = params.mp*ymed(j)*sin(xmed(i))/(params.a*params.a); 
            resj += ymed(j)*(cos(xmed(i))*resy - sin(xmed(i))*resx); 
            resj *= -2*M_PI*ymed(j)*rho[l];
            res += resj;
            
        }
        res /= (double)nx;
        lamind[j] = res;

    }
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

	rho[lgh] = (rho[i+(ny+NGHY-1)*pitch]-rho[i+(ny+NGHY-2)*pitch])/(ymed(ny+NGHY-1)-ymed(ny+NGHY-2))*(ymed(jgh)-ymed(ny+NGHY-1))+rho[i+(ny+NGHY-1)*pitch];
	vx[lgh] = (vx[lactb])*sqrt(ymed(jact)/ymed(jgh));
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

	rho[lgh] = (rho[i+(NGHY+1)*pitch]-rho[i+NGHY*pitch])/(ymed(NGHY+1)-ymed(NGHY))*(ymed(jgh)-ymed(NGHY))+rho[i+NGHY*pitch];
	vx[lgh] = (vx[lactb])*sqrt(ymed(jact)/ymed(jgh));
	vy[lghs] = -vy[lactbs];
	vy[lactbs_null] = 0;
//<\#>
      }
    }
  }
    return;
}

void add_boundary(void) {
    int i,j,k;
    k = 0;
#ifdef SYMMETRIC
    ymin_bound();
    ymax_bound();
#else
    double fh1;
    double fac,rhom;
    double facm;
    int jact;
    jact = NGHY;
    for(j=0;j<NGHY;j++) {
            for(i=0;i<size_x;i++) {
                vy[l] = vy[lact]*pow(ymin(j)/ymin(jact),params.vrindex);
                rho[l] = rho[lact]*pow(ymed(j)/ymed(jact),-params.nuindex);
                vx[l] = vx[lact]*pow(ymed(j)/ymed(jact),-.5);

            }
    }
    jact=size_y-NGHY-1;
    for(j=size_y-NGHY;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                fh1 = rho[lact]*3*M_PI*Nu(ymed(jact))*sqrt(ymed(jact));
                fac = 3*M_PI*Nu(ymed(j))*sqrt(ymed(j));
                facm = 3*M_PI*Nu(ymin(j))*sqrt(ymin(j));
                vx[l] = vx[lact]*pow(ymed(j)/ymed(jact),-.5);
                rho[l] = (fh1 + params.mdot*(sqrt(ymed(j))-sqrt(ymed(jact)))/(3*M_PI))/fac;
                rhom = (fh1 + params.mdot*(sqrt(ymin(j))-sqrt(ymed(jact)))/(3*M_PI))/facm;
                vy[l] = params.mdot/(-2*M_PI*ymin(j)*rhom);
                //vy[l] = vy[lact]*pow(ymin(j)/ymin(jact),params.vrindex);

            }
    }
#endif
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
                fread(&rho[l],sizeof(double),1,fd);
                fread(&vx[l],sizeof(double),1,fx);
                fread(&vy[l],sizeof(double),1,fy);
                vx[l] += params.omf*ymed(j);
            }
        }
    }



    fclose(fd);
    fclose(fx);
    fclose(fy);
    return;

}
void set_averages(void) {
    int i,j,k;
    k = 0;

    double res,res1,res2;


    for(j=0;j<size_y;j++) {
        res = 0;
        res1 = 0;
        res2 = 0;

        for(i=0;i<size_x;i++) {
            res += rho[l];
            res1 += vx[l];
            res2 += vy[l];
        }
        rhobar[j] = res/(double)nx;
        vxbar[j] = res1/(double)nx;
        vybar[j] = res2/(double)nx;
    
        for(i=0;i<size_x;i++) {
            sig[l] = rho[l] - rhobar[j];
            vxp[l] = vx[l] - vxbar[j];
            vyp[l] = vy[l] - vybar[j];
        }
    }
    return;
}

void set_slopesbar(double *q) {
    int i,j,k;
    k=0;
    double dqm, dqp;

    for (j=1; j<size_y-1; j++) {
        dqm = (q[j]-q[j-1])/zone_size_y(j,k);
        dqp = (q[j+1]-q[j])/zone_size_y(j+1,k);
        if(dqp*dqm<=0) work[j] = 0;
        else  work[j] = 2.*dqp*dqm/(dqm+dqp);
      }
    
    work[0] = 0; work[size_y-1]= 0;
    
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
    j=0;
    for(i=0;i<size_x;i++) {
        work[l] = 0;
    }
    j=size_y-1;
    for(i=0;i<size_x;i++) {
        work[l] = 0;
    }

    
    return;
}
void transportbar(double *qs,double *vel,double *res) {
    int i,j,k;
    k = 0;
   
    if (res == NULL) {
        res = qs;
    }

    for(j=0;j<size_y;j++) {
        work1[j] = qs[j];
    }
    set_slopesbar(qs);   // Slope is in work array

    for (j=1; j<size_y-1; j++) {

    	if (vel[j]>0.) {
            res[j] = work1[j-1] + 0.5 * (zone_size_y(j-1,k))*work[j-1];
        }
	    else {
	        res[j] = work1[j] - 0.5 * (zone_size_y(j,k))*work[j];
        }


    }
    return;
}
void transport(double *qs,double *vel,double *res) {
    int i,j,k;
    k = 0;
   
    if (res == NULL) {
        res = qs;
    }

    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            work1[l] = qs[l];
        }
    }
    set_slopes(qs);   // Slope is in work array

    for (j=1; j<size_y-1; j++) {
      for (i=0; i<size_x; i++) {

    	if (vel[l]>0.) {
            res[l] = work1[lym] + 0.5 * (zone_size_y(j-1,k))*work[lym];
        }
	    else {
	        res[l] = work1[l] - 0.5 * (zone_size_y(j,k))*work[l];
        }

      }

    }
    return;
}
void clear_work(void){
    int i,j,k;
    k=0;
    for(k=0;k<size_z;k++) {
        for(j=0;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                work[l] = 0;
                work1[l] =0;
            }
        }
    }
    return;
}
void output(char *directory) {
    int i,j,k;
    k=0;
    FILE *f;
    char fname[256];
    sprintf(fname,"%stemp_files/output.dat",directory);
    f = fopen(fname,"w");
    fwrite(&Ymin[NGHY],sizeof(double),ny,f);
    fwrite(&Ymed[NGHY],sizeof(double),ny,f);
    fwrite(&Ltot[NGHY],sizeof(double),ny,f);
    fwrite(&Ld[NGHY],sizeof(double),ny,f);
    fwrite(&Lw[NGHY],sizeof(double),ny,f);
    fwrite(&drFtot[NGHY],sizeof(double),ny,f);
    fwrite(&drFd[NGHY],sizeof(double),ny,f);
    fwrite(&drFw[NGHY],sizeof(double),ny,f);
    fwrite(&lamex[NGHY],sizeof(double),ny,f);
    fwrite(&lamind[NGHY],sizeof(double),ny,f);
    fwrite(&lamdep[NGHY],sizeof(double),ny,f);
    fwrite(&Mdottot[NGHY],sizeof(double),ny,f);
    fwrite(&Mdotd[NGHY],sizeof(double),ny,f);
    fwrite(&Ldstar[NGHY],sizeof(double),ny,f);
    fwrite(&vrbar[NGHY],sizeof(double),ny,f);
    fclose(f);
    return;

}
