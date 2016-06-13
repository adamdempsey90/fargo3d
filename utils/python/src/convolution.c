#include "evolve.h"
#include "rfftw.h"


double *in1, *in2;
fftw_real *out1, *out2;
int size_f;
double fft_norm_fac;
rfftw_plan fplan;


void fft1d(const double *fld1, double *res) {
    memcpy(&in1[0],&fld1[0],sizeof(double)*size_x);
    rfftw_one(fplan,in1,out1);
    memcpy(&res[0],&out1[0],sizeof(double)*size_x);
    return;
}


void convolution(const double *fld1, const double *fld2, double *res, double fac, int jres, int ncols) {
    int mi;
    double temp;
    
    memcpy(&in1[0],&fld1[0],sizeof(double)*size_x);
    memcpy(&in2[0],&fld2[0],sizeof(double)*size_x);


    rfftw_one(fplan,in1,out1);
    rfftw_one(fplan,in2,out2);

    res[jres + ncols] += fac*out1[0]*out2[0]/fft_norm_fac; 
    
    for(mi=1;mi<MMAX;mi++) {
        res[jres + (mi+1)*ncols] += fac*2*(out1[mi]*out2[mi] + out1[size_x-mi]*out2[size_x-mi])/fft_norm_fac;
    }
    temp = 0;
    for(mi=MMAX; mi<size_f;mi++) {
        temp += 2*(out1[mi]*out2[mi] + out1[size_x-mi]*out2[size_x-mi]);
    }
    if ( NEVEN ) {
        temp += out1[size_x/2]*out2[size_x/2];
    }
    res[jres + (MMAX+1)*ncols] += fac*temp/fft_norm_fac;
    

    return;
}
void convolution_deriv(const double *fld1, const double *fld2, double *res, double fac, int jres,int ncols) {
    int mi;
    memcpy(&in1[0],&fld1[0],sizeof(double)*size_x);
    memcpy(&in2[0],&fld2[0],sizeof(double)*size_x);
    rfftw_one(fplan,in1,out1);
    rfftw_one(fplan,in2,out2);

    for(mi=1;mi<MMAX;mi++) {
        res[jres + (mi+1)*ncols] += mi*fac*2*(out1[mi]*out2[size_x-mi] - out1[size_x-mi]*out2[mi])/fft_norm_fac;
    }
    double temp = 0;
    for(mi=MMAX;mi< size_f;mi++) {
        temp += mi*2*(out1[mi]*out2[size_x-mi] - out1[size_x-mi]*out2[mi]);
    }
    res[jres + (MMAX+1)*ncols] += fac*temp/fft_norm_fac;


    return;
}

void allocate_conv(void) {
    size_f = (size_x + 1)/2;
    in1 = (double *)malloc(sizeof(double)*size_x);
    in2 = (double *)malloc(sizeof(double)*size_x);
    out1 = (fftw_real *)malloc(sizeof(double)*size_x);
    out2 = (fftw_real *)malloc(sizeof(double)*size_x);
    
    fplan = rfftw_create_plan(size_x,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE);

    fft_norm_fac = (double)size_x*size_x;
    int i;
    for(i=0;i<size_x;i++) {
        in1[i] = 0;
        in2[i] = 0;
        out1[i] = 0;
        out2[i] = 0;
    }

    return;
}
void free_conv(void) {
    free(in1);
    free(in2);
    fftw_free(out1);
    fftw_free(out2);
    rfftw_destroy_plan(fplan);
    return;
}
