#include "lib_poisson1D.h"

void GB_to_full(double **R, double *A, int la);
void print_full(double **R, int la);
void my_GB_tri_facto_lu(double * A, int la);
void my_GB_tri_up(int lab, int la, double *A, double *x, double *b);
void my_GB_tri_down(int lab, int la, double *A, double *x, double *b);
void my_GB_tri_LU(int lab, int la, double *A, double *x, double *b);