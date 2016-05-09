#include "evolve.h"
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
    fscanf(f,"%lg\n",&params.h);
    fscanf(f,"%lg\n",&params.flaringindex);
    fscanf(f,"%lg\n",&params.mdot);
    fscanf(f,"%lg\n",&params.soft);
    fscanf(f,"%d\n",&IndirectTerm);

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
    int i,j;
    double xpl,ypl,zpl,vxpl,vypl,vzpl,mpl,tpl,omfpl;
    double epl,apl,vpl,psipl,phipl,ipl,wpl,alphapl;
    int scanned = 0;

    for(j=0;j<nb;j++) {
        sprintf(filename,"%splanet%d.dat",directory,j);
        f = fopen(filename,"r");
        while ( (scanned = fscanf(f,"%d\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    &i,&xpl,&ypl,&zpl,&vxpl,&vypl,&vzpl,&mpl,&tpl,&omfpl)) != EOF) {

                if (i == n) {
                    psys[j].x = xpl; 
                    psys[j].y = ypl; 
                    psys[j].z = zpl; 
                    psys[j].vx = vxpl; 
                    psys[j].vy = vypl; 
                    psys[j].vz = vzpl; 
                    psys[j].mp = mpl; 
                    psys[j].omf = omfpl; 
                    psys[j].t = tpl;
                    break;
                }

        }


        fclose(f);
        sprintf(filename2,"%sorbit%d.dat",directory,j);
        f1 = fopen(filename2,"r");

        while ( (scanned = fscanf(f1,"%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\t%lg\n",
                    &tpl,&epl,&apl,&mpl,&vpl,&psipl,&phipl,&ipl,&wpl,&alphapl)) != EOF) {

                if (fabs(tpl-psys[j].t)<1e-8) {
                    psys[j].orbit.e = epl; 
                    psys[j].orbit.a = apl; 
                    psys[j].orbit.M = mpl; 
                    psys[j].orbit.V = vpl; 
                    psys[j].orbit.psi = psipl; 
                    psys[j].orbit.phi = phipl; 
                    psys[j].orbit.i = ipl; 
                    psys[j].orbit.w = wpl; 
                    psys[j].orbit.alpha = alphapl; 
                    break;
                }

        }
        fclose(f1);
    }
    omf = psys[0].omf;


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
