#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NGHY 3

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
double *Ymin, *Xmin, *Ymed, *Xmed;
double *rho, *vx, *vy;
double *rhos, *mdot;
double *momp, *momm, *lstar;
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
void set_lam_ex(char *directory);
void set_tensor(char *directory);
void set_mdot(char *directory);
void set_slopes(double *q); 
void set_rhostar(char *directory);
void clear_work(void);
void set_Fd(char *directory);
void set_Twd(char *directory);
void set_lstar(char *directory);
void transport(double *q);
void set_Ld(char *directory);

int main(int arc, char *argv[]) { 
  int n;
  int i,j,k;
  char directory[256];
  char param_fname[100];

  n = atoi(argv[1]);
  strcpy(directory,argv[2]);
  sprintf(param_fname,"%sparam_file.txt",directory);
//strcpy(param_fname,argv[2]);

  read_param_file(param_fname);

  

  size_x = nx;
  size_y = ny+2*NGHY;
  size_z = nz;
  stride = size_x*size_y;
  pitch = size_x;
  dx = 2*M_PI/nx;



  Xmin  =     (double*)malloc(sizeof(double)*(size_x+1));
  Ymin  =     (double*)malloc(sizeof(double)*(size_y+1));
  Xmed  =     (double*)malloc(sizeof(double)*(size_x));
  Ymed  =     (double*)malloc(sizeof(double)*(size_y));

  vx =     (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  vy =     (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  rho = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  rhos = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  momp = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  momm = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  lstar = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  mdot = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  
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
            momp[l] = 0;
            momm[l] = 0;
            lstar[l] = 0;
      }
    }
  }

  read_domain(directory); 
  read_files(n,directory);
  add_boundary();
  set_rhostar(directory);
  set_mdot(directory);
  set_lstar(directory);
  clear_work();
  set_tensor(directory);
  set_lam_ex(directory);
  clear_work();
  set_Fd(directory);
  set_Twd(directory);
  set_Ld(directory);

  free(Xmin);
  free(Ymin);
  free(Xmed);
  free(Ymed);
  free(vx);
  free(vy);
  free(rho);
  free(rhos);
  free(tauxx);
  free(tauyy);
  free(tauxy);
  free(work);
  free(work1);
  free(lstar);
  free(momm);
  free(momp);
  return 1;

}
double Nu(double x) {
    return params.alpha*params.h*params.h*pow(x,params.nuindex);
}
double Cs(double x) {
    return params.h*pow(x,params.flaringindex-0.5);
}

void set_Ld(char *directory) {
    int i,j,k;
    k=0;

    double Ld,Lw,Ltot;
    double rhoavg, lavg;
    double fac;
    char filename[256];
    sprintf(filename,"%stemp_files/Ld.dat",directory);
    printf("Writing %s\n",filename);
    FILE *f= fopen(filename,"w");

/*
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {
            fwrite(&work[l],sizeof(double),1,f);
        }
    }
*/

    char filename2[256]; 
    sprintf(filename2,"%stemp_files/Lw.dat",directory);
    printf("Writing %s\n",filename2);
    FILE *f1= fopen(filename2,"w");

    for(j=NGHY;j<size_y-NGHY;j++) {
        rhoavg = 0;
        lavg = 0;
        fac  = 2*M_PI*ymed(j);
        for(i=0;i<size_x;i++) {
            rhoavg += rho[l];
            lavg += .5*(vx[l]+vx[lxm])*ymed(j);
        }
        rhoavg /=(double)nx;
        lavg /=(double)nx;

        for(i=0;i<size_x;i++) {
            
            Ld = fac*rhoavg*lavg;
            Lw = fac*(rho[l]-rhoavg)*( (vx[l]+vx[lxm])*.5*ymed(j)-lavg);


            fwrite(&Ld,sizeof(double),1,f);
            fwrite(&Lw,sizeof(double),1,f1);
        }

        }
    fclose(f);
    fclose(f1);

    return;
}
void set_tensor(char *directory) {
    int i,j,k;
    double visc, viscm, div_v;
    double cs, csm;
    k = 0;

    for(j=1;j<size_y-1;j++) {
        cs = Cs(ymed(j));
        csm = Cs(ymin(j));
        cs *= cs; csm *= csm;
        for(i=0;i<size_x;i++) {
            visc = Nu(ymed(j));
            viscm = Nu(ymin(j));

            div_v = (vx[lxp]-vx[l])*SurfX(j,k);
            div_v += (vy[lyp]*SurfY(j+1,k)-vy[l]*SurfY(j,k));
            div_v *= 2.0/3.0*InvVol(j,k);

            tauxx[l] = visc*rho[l]*(2.0*(vx[lxp]-vx[l])/zone_size_x(j,k) - div_v);
            tauxx[l] += visc*rho[l]*(vy[lyp]+vy[l])/ymed(j);
            tauxx[l] += -cs*rho[l]; // Isothermal Pressure
            tauyy[l] = visc*rho[l]*(2.0*(vy[lyp]-vy[l])/(ymin(j+1)-ymin(j)) - div_v);
            tauyy[l] += -cs*rho[l]; // Isothermal Pressure


           // tauxy[l] = viscm*.5*(rhos[l]+rhos[lxm])*((vy[l]-vy[lxm])/(dx*ymin(j)) + (vx[l]-vx[lym])/(ymed(j)-ymed(j-1))-.5*(vx[l]+vx[lym])/ymin(j)); //centered on left, inner vertical edge in z
            tauxy[l] = viscm*.25*(rho[l]+rho[lxm]+rho[lym]+rho[lxm-pitch])*((vy[l]-vy[lxm])/(dx*ymin(j)) + (vx[l]-vx[lym])/(ymed(j)-ymed(j-1))-.5*(vx[l]+vx[lym])/ymin(j)); //centered on left, inner vertical edge in z
            

        }
    }

    char filename[256];
    sprintf(filename,"%stemp_files/tensor.dat",directory);
    printf("Writing %s\n",filename);
    FILE *f= fopen(filename,"w");
    if (f==NULL) printf("Error loading temp_files/tensor.dat\n");
    k=0;

    fwrite(&tauxx[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
    fwrite(&tauyy[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
    fwrite(&tauxy[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
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
void transport(double *qs) {
    int i,j,k;
    FILE *f;
    k = 0;
    
    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            work1[l] = qs[l];
        }
    }
    set_slopes(qs);   // Slope is in work array

    for (j=1; j<size_y-1; j++) {
      for (i=0; i<size_x; i++) {

    	if (vy[l]>0.) {
            qs[l] = work1[lym] + 0.5 * (zone_size_y(j-1,k))*work[lym];
        }
	    else {
	        qs[l] = work1[l] - 0.5 * (zone_size_y(j,k))*work[l];
        }

      }

    }
    return;
}
    
void set_rhostar(char *directory) {
    int i,j,k;
    k = 0;
    set_slopes(rho);

    for (j=1; j<size_y-1; j++) {
      for (i=0; i<size_x; i++) {

    	if(vy[l]>0.)
	  rhos[l] = rho[lym] + 0.5 * (zone_size_y(j-1,k))*work[lym];
	else
	  rhos[l] = rho[l] - 0.5 * (zone_size_y(j,k))*work[l];

//        rhos[l] = (rho[lym]+rho[lyp])*.5;
      }

    }
    char filename[256];
    sprintf(filename,"%stemp_files/rhostar.dat",directory);
    printf("Writing %s\n",filename);
    FILE *f= fopen(filename,"w");
    fwrite(&rhos[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
    fclose(f);
    return;
}

void set_lstar(char *directory) {
    int i,j,k;
    k=0;
    transport(momp);
    transport(momm);
    // Now momenta are on the edge.

    for(j=0;j<size_y;j++) {
        for(i=0;i<size_x;i++) {
            lstar[l] = ymin(j)*.5*(momp[l] + momm[l])/rhos[l];
        }
    }

    char filename[256];
    sprintf(filename,"%stemp_files/lstar.dat",directory);
    printf("Writing %s\n",filename);
    FILE *f= fopen(filename,"w");
    fwrite(&lstar[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
    fclose(f);

    return;

}

void set_mdot(char *directory) {
    int i,j,k;
    k = 0;

    

    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {    
	        mdot[l] = -nx*vy[l]*rhos[l]*SurfY(j,k) ;
        }
    }
    char filename[256];
    sprintf(filename,"%stemp_files/mdot.dat",directory);
    printf("Writing %s\n",filename);
    FILE *f= fopen(filename,"w");

    fwrite(&mdot[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
    fclose(f);
    return;

}


void set_lam_ex(char *directory) {
    int i,j,k;

    double smoothing = params.soft*params.h*pow(params.a,1 + params.flaringindex);
    smoothing *= smoothing;
    double rad,res;
    double pot;

    char filename[256];
    sprintf(filename,"%stemp_files/lambda_ex.dat",directory);
    printf("Writing %s\n",filename);
    FILE *f= fopen(filename,"w");

    k = 0;
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {    
            rad = ymin(j)*ymin(j) + params.a*params.a -2*params.a*ymin(j)*cos(xmed(i))+smoothing;
            rad *= sqrt(rad);
            pot = params.mp * ymin(j)*params.a*sin(xmed(i))/rad;
            res = -2*M_PI*ymin(j)*rhos[l]*pot;
            fwrite(&res,sizeof(double),1,f);

        }
    }
    fclose(f);
    return;

}

void set_Twd(char *directory) {
    int i,j,k;
    double lc;
    k=0;

    double *mdbar = (double *)malloc(sizeof(double)*(size_y));
    double *rhobar = (double *)malloc(sizeof(double)*(size_y));
    double *lbar = (double *)malloc(sizeof(double)*(size_y));
    double *dlbar = (double *)malloc(sizeof(double)*(size_y));
    double *vybar = (double *)malloc(sizeof(double)*(size_y));
    double *dbar = (double *)malloc(sizeof(double)*(size_y));
    double res1, res2,res3, res4,res5;
    double dlbarc, lbarc;
    for(j=0;j<size_y;j++) {
        res1=0;
        res2=0;
        res3=0;
        res4=0;
        res5=0;
        for(i=0;i<size_x;i++) {
            res1 += mdot[l];
            res2 += lstar[l];
            res3 += rhos[l];
            res4 += vy[l];
            res5 += rho[l];
        }
        mdbar[j] = res1/(double)nx;
        lbar[j] = res2/(double)nx;
        rhobar[j] = res3/(double)nx;
        vybar[j] = res4/(double)nx;
        dbar[j] = res5/(double)nx;
    }
    for(j=1;j<size_y-1;j++) {
        dlbar[j] = (lbar[j+1]-lbar[j-1])/(ymin(j+1)-ymin(j-1));
    }
    dlbar[0] = (lbar[1]-lbar[0])/(ymin(1)-ymin(0));
    dlbar[size_y-1] = (lbar[size_y-1]-lbar[size_y-2])/(ymin(size_y-1)-ymin(size_y-2));

    double fac,fac1;
    double rho_bar, rho_edge;
    for(j=NGHY;j<size_y-NGHY;j++) {
        rho_bar = rhobar[j];
        for(i=0;i<size_x;i++) {
            rho_edge = rhos[l]; //(rho[lym]+rho[l])*.5;
            fac1 = M_PI*(ymin(j+1)*ymin(j+1)*(tauxy[lyp]+tauxy[lxp+pitch])-ymin(j-1)*ymin(j-1)*(tauxy[lym]+tauxy[lxp-pitch]))/(ymin(j+1)-ymin(j-1));
            work[l] = ( .5*M_PI*ymin(j)*(tauxx[lxp]+tauxx[lxp-pitch]-(tauxx[lxm]+tauxx[lxm-pitch]))/(dx) + fac1)*rho_bar/rho_edge;//(dbar[j]+dbar[j-1])*.5;
            work[l] -= fac1;

            work[l] += 2*M_PI*ymin(j)*(rhos[l]-rhobar[j])*(vy[l]-vybar[j])*(lbar[j+1]-lbar[j-1])/(ymin(j+1)-ymin(j-1));
            work[l] -= 2*M_PI*ymin(j)*rhobar[j]* (vy[l]-vybar[j])*(lstar[lyp]-lbar[j+1] - (lstar[lym]-lbar[j-1]))/(ymin(j+1)-ymin(j-1));
            //work[l] += (-mdot[l]-rho_bar*vybar[j]*2*M_PI*ymin(j))*dlbar[j];
           // work[l] -= 2*M_PI*ymin(j)*rho_bar*(vy[l]*(lstar[lyp]-lstar[lym])/(ymin(j+1)-ymin(j)) - vybar[j]*dlbar[j]);
            //work[l] -= 2*M_PI*ymin(j)*rho_bar*(vy[l]*.5*(ymed(j)*vx[l]-ymed(j-1)*vx[lym]+ymed(j)*vx[lxp]-ymed(j-1)*vx[lxp-pitch])/(ymed(j)-ymed(j-1)) - vybar[j]*dlbar[j]);

//	        work1[l] = 2.0*(tauxx[l]-tauxx[lxm])/(zone_size_x(j,k)*(rho[l]+rho[lxm]));
//            work1[l] += 2.0*(ymin(j+1)*ymin(j+1)*tauxy[lyp]-ymin(j)*ymin(j)*tauxy[l])/((ymin(j+1)-ymin(j))*ymed(j)*ymed(j)*(rho[lxm]+rho[l]));
        }
    }

    char filename[256];
    sprintf(filename,"%stemp_files/lambda_dep.dat",directory);
    printf("Writing %s\n",filename);
    FILE *f= fopen(filename,"w");

    fwrite(&work[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
/*
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {
            fwrite(&work[l],sizeof(double),1,f);
        }
    }
*/
    fclose(f);
    free(mdbar);
    free(lbar);
    free(dlbar);
    free(rhobar);
    free(vybar);
    free(dbar);
    return;

}


void set_Fd(char *directory) {

    int i,j,k;
    double lc;
    k=0;

    double *mdbar = (double *)malloc(sizeof(double)*(size_y));
    double *lbar = (double *)malloc(sizeof(double)*(size_y));
    double res1, res2;
    for(j=0;j<size_y;j++) {
        res1=0;
        res2=0;
        
        for(i=0;i<size_x;i++) {
            res1 += mdot[l];
            res2 += lstar[l];
        }
        res1 /= (double)nx;
        res2 /= (double)nx;
        mdbar[j] = res1;//res1/(double)nx;
        lbar[j] = res2;//res2/(double)nx;
    }
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {
            work[l] = -mdbar[j]*lbar[j] - M_PI*ymin(j)*ymin(j)*(tauxy[l]+tauxy[lxp]);
            //work[l] = -mdbar[j]*.5*(lbar[j]+lbar[j+1]) - M_PI*ymin(j)*ymin(j)*(tauxy[l]+tauxy[lxp]);

            work1[l] = -( mdot[l]-mdbar[j])*(lstar[l]-lbar[j]);
        }
    }

    char filename[256];
    sprintf(filename,"%stemp_files/fd.dat",directory);
    printf("Writing %s\n",filename);
    FILE *f= fopen(filename,"w");

    fwrite(&work[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);


/*
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {
            fwrite(&work[l],sizeof(double),1,f);
        }
    }
*/
    fclose(f);

    char filename2[256]; 
    sprintf(filename2,"%stemp_files/fw.dat",directory);
    printf("Writing %s\n",filename);
    f= fopen(filename2,"w");

    fwrite(&work1[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
/*
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {
            fwrite(&work1[l],sizeof(double),1,f);
        }
    }
*/
    fclose(f);


    free(mdbar);
    free(lbar);
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


void add_boundary(void) {
    double fh1;
    double fac,rhom;
    double facm;
    int i,j,k,jact;
    k = 0;
    
    jact = NGHY;
    for(j=0;j<NGHY;j++) {
            for(i=0;i<size_x;i++) {
                rho[l] = rho[lact]*pow(ymed(j)/ymed(jact),-params.nuindex);
                vx[l] = vx[lact]*pow(ymed(j)/ymed(jact),-.5);
                vy[l] = vy[lact]*pow(ymin(j)/ymin(jact),params.vrindex);

            }
    }
    jact=size_y-NGHY-1;
    for(j=size_y-NGHY;j<size_y;j++) {
            for(i=0;i<size_x;i++) {
                fh1 = rho[lact]*3*M_PI*Nu(ymed(jact))*sqrt(ymed(jact));
                fac = 3*M_PI*Nu(ymed(j))*sqrt(ymed(j));
                facm = 3*M_PI*Nu(ymin(j))*sqrt(ymin(j));
                
                rho[l] = (fh1 + params.mdot*(sqrt(ymed(j))-sqrt(ymed(jact)))/(3*M_PI))/fac;
                rhom = (fh1 + params.mdot*(sqrt(ymin(j))-sqrt(ymed(jact)))/(3*M_PI))/facm;
                vy[l] = params.mdot/(-2*M_PI*ymin(j)*rhom);
                vx[l] = vx[lact]*pow(ymed(j)/ymed(jact),-.5);
                //vy[l] = vy[lact]*pow(ymin(j)/ymin(jact),params.vrindex);

            }
    }
    for(j=0;j<size_y;j++) {
       for(i=0;i<size_x;i++) {
        momp[l] = vx[lxp]*rho[l]; 
        momm[l] = vx[l]*rho[l];
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
/*
    for(j=0;j<ny+2*NGHY;j++) {
        if ( (j < NGHY) || j >= (ny+1+2*NGHY-NGHY) ) {
            fscanf(fy,"%lg\n",&temp);
        }
        else {
            fscanf(fy,"%lg\n",&Ymin[j-2]);
        }
    }
*/
   

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




