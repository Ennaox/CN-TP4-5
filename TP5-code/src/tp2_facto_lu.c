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

  	set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);

  	double * A = malloc(sizeof(double)*lab*la);
	set_GB_operator_rowMajor_poisson1D(A, &lab, &la,&kv);

	double **AF = malloc(la*sizeof(double *));

	my_GB_tri_LU(la,la,A,x,RHS);

	write_vec(RHS, &la, "LU.dat");

	//Calcul de l'erreur
	print_GB(la,EX_SOL);
	printf("\n");
	print_GB(la,RHS);

	double temp, relres;
	temp = cblas_ddot(la, RHS, 1, RHS,1);
  	temp = sqrt(temp);
  	cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1);
  	relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
  	relres = sqrt(relres);
  	relres = relres / temp;
  	printf("\nThe relative residual error is relres = %e\n",relres);


  	free(EX_SOL);
  	free(X);
  	free(RHS);
  	free(x);
  	free(A);
	return 0;
}