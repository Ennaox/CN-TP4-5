/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"

int main(int argc,char *argv[])
/* ** argc: Number of arguments */
/* ** argv: Values of arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  double *RHS, *EX_SOL, *X, *RHS_cpy;
  double *AB;

  double temp, relres;

  NRHS=1;
  nbpoints=102;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);
  RHS_cpy = malloc(la * sizeof(double));

  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_dense_RHS_DBC_1D(RHS_cpy,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);

  info=1;

  /* working array for pivot used by LU Factorization */
  ipiv = (int *) calloc(la, sizeof(int));

  int row = 0; //
  double alpha = 1;
  double beta = -1;
  int incx = 1;
  int incy = 1;

  if (row == 1){ // LAPACK_ROW_MAJOR
    //Calcul de la solution de A * x = B pour une matrice poisson 1D en format GB Row Major
    set_GB_operator_rowMajor_poisson1D(AB, &lab, &la,&kv);
    write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB_row.dat");
    
    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR,la, kl, ku, NRHS, AB, la, ipiv, RHS, NRHS);
    
    //Calcul de y = alpha * A + beta * y pour une matrice poisson 1D en format GB Row Major
    free(AB);
    kv = 0;
    lab=kv+kl+ku+1;

    AB = (double *) malloc(sizeof(double)*lab*la);
    set_GB_operator_rowMajor_poisson1D(AB, &lab, &la,&kv);

    cblas_dgbmv(CblasRowMajor,CblasNoTrans,la,la,kl,ku,alpha,AB,lab,EX_SOL,incx,beta,RHS_cpy,incy);
    temp = cblas_ddot(la, RHS_cpy, 1, RHS_cpy,1);
    temp = sqrt(temp);
    printf("Error for dgbmv: %e\n",temp);


  }
  else { // LAPACK_COL_MAJOR
    //Calcul de la solution de A * x = B pour une matrice poisson 1D en format GB Col Major
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_col.dat");

    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR,la, kl, ku, NRHS, AB, lab, ipiv, RHS, la);

    //Calcul de y = alpha * A + beta * y pour une matrice poisson 1D en format GB Col Major
    free(AB);
    kv = 0;
    lab=kv+kl+ku+1;

    AB = (double *) malloc(sizeof(double)*lab*la);
    set_GB_operator_rowMajor_poisson1D(AB, &lab, &la,&kv);

    //Faire appelle a DGBMV
    free(AB); //Réécriture de AB sans la ligne de 0
    kv = 0;
    lab=kv+kl+ku+1;
    AB = (double *) malloc(sizeof(double)*lab*la);

    set_GB_operator_colMajor_poisson1D(AB, &lab, &la,&kv);

    cblas_dgbmv(CblasColMajor,CblasNoTrans,la,la,kl,ku,alpha,AB,lab,EX_SOL,incx,beta,RHS_cpy,incy);
    temp = cblas_ddot(la, RHS_cpy, 1, RHS_cpy,1);
    temp = sqrt(temp);
    printf("Error for dgbmv: %e\n",temp);
  }    

  
  printf("\n INFO DGBSV = %d\n",info);

  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative residual */
  temp = cblas_ddot(la, RHS, 1, RHS,1);
  temp = sqrt(temp);
  cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1);
  relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
  relres = sqrt(relres);
  relres = relres / temp;
  
  printf("\nThe relative residual error for dgbsv is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  free(ipiv);

  printf("\n\n--------- End -----------\n");
}
