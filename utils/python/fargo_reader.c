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
void set_lam_ex(void);
void set_tensor(void);
void set_mdot(void);
void set_slopes(double *q); 
void set_rhostar(void);
void clear_work(void);
void set_Fd(void);
void set_Twd(void);



int main(int arc, char *argv[]) { 
  int n;
  int i,j,k;
  char directory[50];
  char param_fname[100];

  printf("reading arguments\n");
  n = atoi(argv[1]);
  printf("Using %d file\n",n);
  strcpy(param_fname,argv[2]);
  printf("File name %s\n",param_fname);

  //sprintf(directory,"%s","");
  printf("reading param file\n");
  read_param_file(param_fname);

  

  size_x = nx;
  size_y = ny+2*NGHY;
  size_z = nz;
  stride = size_x*size_y;
  pitch = size_x;
  dx = 2*M_PI/nx;



  printf("allocating\n");
  Xmin  =     (double*)malloc(sizeof(double)*(size_x+1));
  Ymin  =     (double*)malloc(sizeof(double)*(size_y+1));
  Xmed  =     (double*)malloc(sizeof(double)*(size_x));
  Ymed  =     (double*)malloc(sizeof(double)*(size_y));

  vx =     (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  vy =     (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  rho = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  rhos = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  mdot = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  
  tauxx =     (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  tauyy =     (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  tauxy = (double*)malloc(sizeof(double)*size_x*size_y*size_z);

  work = (double*)malloc(sizeof(double)*size_x*size_y*size_z);
  work1 = (double*)malloc(sizeof(double)*size_x*size_y*size_z);

  printf("setting i.c\n");
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

  read_domain(NULL); 
  printf("Read domain\n");
  read_files(n,NULL);
  printf("Loaded data\n");
  add_boundary();
  printf("Added b.c\n");
  set_tensor();
  printf("Added tensor\n");
  set_mdot();
  printf("Added mdot\n");
  set_lam_ex();
  clear_work();
  printf("Added dTr\n");
  set_Fd();
  printf("Added Fd\n");
  set_Twd();
  printf("Added Twd\n");
  
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
  return 1;

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
    printf("Finished 1.\n");
    j=0;
    for(i=0;i<size_x;i++) {
        tauxx[l] = tauxx[lyp];
        tauyy[l] = tauyy[lyp];
        tauxy[l] = tauxy[lyp];
    }
    printf("Finished 2.\n");
    j=size_y-1;
    for(i=0;i<size_x;i++) {
        tauxx[l] = tauxx[lym];
        tauyy[l] = tauyy[lym];
        tauxy[l] = tauxy[lym];
    }
    printf("Finished 3.\n");

    FILE *f= fopen("temp_files/tensor.dat","w");
    if (f==NULL) printf("Error loading temp_files/tensor.dat\n");
    k=0;

    fwrite(&tauxx[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
    fwrite(&tauyy[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
    fwrite(&tauxy[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
/*
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {        
            fwrite(&tauxx[l],sizeof(double),1,f);
        }
    }
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {        
            fwrite(&tauyy[l],sizeof(double),1,f);
        }
    }
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {        
            fwrite(&tauxy[l],sizeof(double),1,f);
        }
    }
*/
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
void set_rhostar(void) {
    int i,j,k;
    FILE *f;
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
    f = fopen("temp_files/rhostar.dat","w");
    fwrite(&rhos[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
/*
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {
            fwrite(&rhos[l],sizeof(double),1,f);
        }
    }
*/
    fclose(f);
    return;
}

void set_mdot(void) {
    int i,j,k;
    k = 0;
    FILE *f;

    printf("set slopes.\n");
    set_rhostar();
    printf("set rhostar\n");
    

    for(j=1;j<size_y-1;j++) {
        for(i=0;i<size_x;i++) {    
	        mdot[l] = -nx*vy[l]*rhos[l]*SurfY(j,k) ;
        }
    }
    f = fopen("temp_files/mdot.dat","w");

    fwrite(&mdot[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
/*
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) { 
            fwrite(&mdot[l],sizeof(double),1,f);
        }
    }
*/
    fclose(f);
    return;

}


void set_lam_ex(void) {
    int i,j,k;
    FILE *f;
    double smoothing = params.soft*params.h*pow(params.a,1 + params.flaringindex);
    smoothing *= smoothing;
    double rad,res;
    double pot;

    f = fopen("temp_files/lambda_ex.dat","w");

    k = 0;
    for(j=NGHY;j<size_y-NGHY;j++) {
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

void set_Twd(void) {
    FILE *f;
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
            res2 += vx[l]*ymed(j);
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
        dlbar[j] = (lbar[j]-lbar[j-1])/(ymed(j)-ymed(j-1));
    }
    dlbar[0] = (lbar[1]-lbar[0])/(ymed(1)-ymed(0));
    dlbar[size_y-1] = (lbar[size_y-1]-lbar[size_y-2])/(ymed(size_y-1)-ymed(size_y-2));

    double fac,fac1;
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {
            fac = (rho[lym]+rho[l])*.5;
            fac1 = .5*(ymin(j+1)*ymin(j+1)*(tauxy[lyp]+tauxy[lxp+pitch])-ymin(j-1)*ymin(j-1)*(tauxy[lym]+tauxy[lxp-pitch]))/(ymin(j+1)-ymin(j-1));
            work[l] = ( .25*(tauxx[lxp]-tauxx[lxm]+tauxx[lxp-pitch]-tauxx[lxm-pitch])/(2*dx*fac) - fac1/(fac*ymin(j)*ymin(j)))*ymin(j)*dbar[j];
            work[l] -= 2*M_PI*fac1;
            work[l] += (-mdot[l]-rhobar[j]*vybar[j])*dlbar[j];
            work[l] -= rhobar[j]*(vy[l]*.5*(ymed(j)*vx[l]-ymed(j-1)*vx[lym]+ymed(j)*vx[lxp]-ymed(j-1)*vx[lxp-pitch])/(ymed(j)-ymed(j-1)) - vybar[j]*dlbar[j]);

//	        work1[l] = 2.0*(tauxx[l]-tauxx[lxm])/(zone_size_x(j,k)*(rho[l]+rho[lxm]));
//            work1[l] += 2.0*(ymin(j+1)*ymin(j+1)*tauxy[lyp]-ymin(j)*ymin(j)*tauxy[l])/((ymin(j+1)-ymin(j))*ymed(j)*ymed(j)*(rho[lxm]+rho[l]));
        }
    }

    f = fopen("temp_files/lambda_dep.dat","w");

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


void set_Fd(void) {
    FILE *f;
    int i,j,k;
    double lc;
    k=0;

    double *mdbar = (double *)malloc(sizeof(double)*(size_y));
    double *lbar = (double *)malloc(sizeof(double)*(size_y));
    double res1, res2;
    printf("Allocated\n");
    for(j=0;j<size_y;j++) {
        res1=0;
        res2=0;
        
        for(i=0;i<size_x;i++) {
            res1 += mdot[l];
            res2 += vx[l]*ymed(j);
        }
        res1 /= (double)nx;
        res2 /= (double)nx;
        mdbar[j] = res1;//res1/(double)nx;
        lbar[j] = res2;//res2/(double)nx;
    }
    printf("Average\n"); 
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {
            work[l] = -mdbar[j]*.5*(lbar[j]+lbar[j+1]) - M_PI*ymin(j)*ymin(j)*(tauxy[l]+tauxy[lxp]);
            work1[l] = -mdot[l] + -mdbar[j]*.5*(lbar[j]+lbar[j+1]);
        }
    }

    printf("opening file\n");
    f = fopen("temp_files/fd.dat","w");

    fwrite(&work[l_f(0,NGHY,0)],sizeof(double),nx*ny,f);
/*
    for(j=NGHY;j<size_y-NGHY;j++) {
        for(i=0;i<size_x;i++) {
            fwrite(&work[l],sizeof(double),1,f);
        }
    }
*/
    fclose(f);

    f = fopen("temp_files/fw.dat","w");

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

    jact = 0;
    for(j=0;j<NGHY;j++) {
        for(k=0; k<size_z; k++) {
            for(i=0;i<size_x;i++) {
                rho[l] = rho[lact]*pow(ymed(j)/ymed(j+1),-params.nuindex);
                vx[l] = vx[lact]*pow(ymed(j)/ymed(j+1),-.5);
                vy[l] = vy[lact]*pow(ymin(j)/ymin(j+1),params.vrindex);

            }
        }
    }
    jact=size_y-NGHY-1;
    for(j=size_y-NGHY;j<size_y;j++) {
        for(k=0; k<size_z; k++) {
            for(i=0;i<size_x;i++) {
                fh1 = rho[lact]*3*M_PI*Nu(ymed(j-1))*sqrt(ymed(j-1));
                fac = 3*M_PI*Nu(ymed(j))*sqrt(ymed(j));
                facm = 3*M_PI*Nu(ymin(j))*sqrt(ymin(j));
                
                rho[l] = (fh1 + params.mdot*(sqrt(ymed(j))-sqrt(ymed(j-1)))/(3*M_PI))/fac;
                rhom = (fh1 + params.mdot*(sqrt(ymin(j))-sqrt(ymed(j-1)))/(3*M_PI))/facm;
                vx[l] = params.mdot/(-2*M_PI*ymin(j)*rhom);
                vy[l] = vy[lact]*pow(ymin(j)/ymin(j-1),params.vrindex);

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
    double temp;
    int i,j;
    printf("Reading domain\n");
    fx = fopen("domain_x.dat","r");
    if (fx == NULL) printf("Error reading %s\n",filename);

    for(i=0;i<size_x+1;i++) {
        fscanf(fx,"%lg\n",&Xmin[i]);
    }
    printf("Read xmin\n");
    fy = fopen("domain_y.dat","r");
    if (fy == NULL) printf("Error reading %s\n",filename);
    for(j=0;j<size_y+1;j++) {
            fscanf(fy,"%lg\n",&Ymin[j]);
    }
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
    printf("Read ymin\n");
   

    for(i=0;i<size_x;i++) {
        Xmed[i] = .5*(Xmin[i] + Xmin[i+1]);
    }
    printf("Read xmed\n");
    for(j=0;j<size_y;j++) {
        Ymed[j] = .5*(Ymin[j] + Ymin[j+1]);
    }
    printf("Read ymed\n");

    fclose(fx);
    fclose(fy);
    return;
}

void read_files(int n, char *directory) {
    char filename[512];
    FILE *fd,*fx,*fy;
    int i,j,k;
    if (directory == NULL) {
    sprintf(filename,"gasdens%d.dat",n);
    fd = fopen(filename,"r");
    if (fd == NULL) printf("Error loading %s\n",filename); 
    sprintf(filename,"gasvx%d.dat",n);
    fx = fopen(filename,"r");
    if (fx == NULL) printf("Error loading %s\n",filename); 
    sprintf(filename,"gasvy%d.dat",n);
    fy = fopen(filename,"r");
    if (fy == NULL) printf("Error loading %s\n",filename); 
    }
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




