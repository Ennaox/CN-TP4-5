#include "lib_lutri.h"

void GB_to_full(double **R, double *A, int la)
{
	for(int i=0;i<la-1;i++)
	{
		R[i+1][i] = A[i+2*la];
		R[i][i] = A[i+la];
		R[i][i+1] = A[i+1];
	}
	R[la-1][la-1] = A[(2*la)-1];
}

void print_full(double **R, int la)
{
	for(int i = 0;i<la;i++)
	{
		for(int j= 0;j<la;j++)
		{
			printf("%lf ",R[i][j]);
		}
		printf("\n");
	}
}

void my_GB_tri_facto_lu(double * A, int la)
{
	for(int i = 0; i<la-1;i++)
	{
		A[2*la + i]=A[2*la + i]/A[la + i];
		A[la + i + 1]=A[la+i+1]-A[i + 1] * A[la*2+i];
	}
}

void my_GB_tri_up(int lab, int la, double *A, double *x, double *b)
{
	//Initialisation de la dernière valeur de x
	x[la-1] = b[la-1]/A[2*la-1];	//La formule est x[la-1] = b[la-1]/U[la-1][la-1] mais 
									// U[la-1][la-1] est le dernière élément sur la diagonale de la matrice

	//Ici i>-1 car on veux traité l'élément x[0]
	for(int i = la-2; i>-1; i--)
	{
		x[i] = (b[i] - A[i] * x[i+1])/A[la+i];
	}
}

void my_GB_tri_down(int lab, int la, double *A, double *x, double *b)
{
	//Initialisation de la première valeur de x
	x[0] = b[0]; //La formule est x[0] = b[0]/L[0] mais ici L[0] = 1

	for(int i = 1; i<la;i++)
	{
		x[i] = b[i] - A[(lab-1)*(i-1)] * x[i-1];
	}
}

void my_GB_tri_LU(int lab, int la, double *A, double *x, double *b)
{
	my_GB_tri_facto_lu(A,la);
	my_GB_tri_down(lab,la,A,x,b);
	my_GB_tri_up(lab,la,A,b,x);
}