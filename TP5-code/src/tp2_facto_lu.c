#include "lib_lutri.h"
#include "lib_poisson1D.h"

int main()
{

	int ku, kl, lab, kv;
	int la;

	kv=0;
	ku=1;
	kl=1;
	lab=kl+ku+1;
	la= 10;
	

	//Préparation de la solution exacte
	double T0=-5.0;
  	double T1=5.0;

  	double *EX_SOL=(double *) malloc(sizeof(double)*la);
  	double *X=(double *) malloc(sizeof(double)*la);

  	set_grid_points_1D(X, &la);
  	set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);


  	//Test résolution de Ax = b par factorisation LU
  	double *RHS=(double *) malloc(sizeof(double)*la);
  	double *x = malloc(sizeof(double)*la);

  	double *poisson1D = malloc(la*lab*sizeof(double));
  	set_GB_operator_rowMajor_poisson1D(poisson1D,&lab,&la,&kv);

  	set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);

  	double * A = malloc(sizeof(double)*lab*la);
	set_GB_operator_rowMajor_poisson1D(A, &lab, &la,&kv);
	
	my_GB_tri_LU(la,la,A,x,RHS);
	A = my_GB_tri_unfacto_lu(A,la,lab);

	//Calcul de l'erreur
	double temp, relres;
	temp = cblas_ddot(la, A, 1, A,1);
  	temp = sqrt(temp);
  	cblas_daxpy(la, -1.0, A, 1, poisson1D, 1);
  	relres = cblas_ddot(la, poisson1D, 1, poisson1D,1);
  	relres = sqrt(relres);
  	relres = relres / temp;
  	printf("The relative residual error is relres = %e\n",relres);


  	free(EX_SOL);
  	free(X);
  	free(RHS);
  	free(x);
  	free(A);
  	free(poisson1D);
	return 0;
}