#include "evolve.h"
#include "rfftw.h"


double *in1, *in2;
fftw_real *out1, *out2;
int size_f;
double fft_norm_fac;
rfftw_plan fplan;



void convolution(double *fld1, double *fld2, double *res, double fac, int jres, int j, int k) {
    int i = 0;
    int mi;
    memcpy(in1,&fld1[l],sizeof(double)*size_x);
    memcpy(in2,&fld2[l],sizeof(double)*size_x);
    rfftw_one(fplan,in1,out1);
    rfftw_one(fplan,in2,out2);

    res[jres + MMAX+2] += fac*out1[0]*out2[0]/fft_norm_fac; 
    
    for(mi=1;mi<MMAX;mi++) {
        res[jres + (mi+1)*(MMAX+2)] += fac*2*(out1[mi]*out2[mi] + out1[size_x-mi]*out2[size_x-mi])/fft_norm_fac;
    }
    double temp = 0;
    for(mi=MMAX;mi< size_f;mi++) {
        temp += 2*(out1[mi]*out2[mi] + out1[size_x-mi]*out2[size_x-mi]);
    }
    res[jres + (MMAX+1)*(MMAX+2)] += fac*temp/fft_norm_fac;

    return;
}
void convolution_deriv(double *fld1, double *fld2, double *res, double fac, int jres, int j, int k) {
    int i = 0;
    int mi;
    memcpy(in1,&fld1[l],sizeof(double)*size_x);
    memcpy(in2,&fld2[l],sizeof(double)*size_x);
    rfftw_one(fplan,in1,out1);
    rfftw_one(fplan,in2,out2);

    
    for(mi=1;mi<MMAX;mi++) {
        res[jres + (mi+1)*(MMAX+2)] += mi*fac*2*(out1[mi]*out2[size_x-mi] - out1[size_x-mi]*out2[mi])/fft_norm_fac;
    }
    double temp = 0;
    for(mi=MMAX;mi< size_f;mi++) {
        temp += mi*2*(out1[mi]*out2[size_x-mi] - out1[size_x-mi]*out2[mi]);
    }
    res[jres + (MMAX+1)*(MMAX+2)] += fac*temp/fft_norm_fac;

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
