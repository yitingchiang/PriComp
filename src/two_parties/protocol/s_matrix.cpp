#include <gmp.h>
#include "../../util.h"
#include "../../smc.h"
#include "s_matrix.h"
#include "s_float.h"

void f_matrix_inverse(SFloat **a, int row, int col, int client_type, SFloat **ia)
{
	int i, j, k, m;

	SFloat norm, term;

	mpz_t temp;

	mpz_init(temp);


	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			if(i != j)
				ia[i][j] = fpShare(0.0, client_type);
			else
			{
				if(client_type == ALICE)
					ia[i][j] = fpShare(1.0, client_type);
				if(client_type == BOB)
					ia[i][j] = fpShare(0.0, client_type);
			}
		}
	}

	printf("--- Gauss-Jordan Elimination ---\n");
	// lower triangle
	for(i = 0; i < row; i++)
	{
		norm = a[i][i];
		for(j = 0; j < col; j++)
		{
			ia[i][j] = f_division(ia[i][j], norm, client_type);
			a[i][j]  = f_division(a[i][j] , norm, client_type);
		}
		for(k = i+1; k < row; k++)
		{
			if(client_type == ALICE)
			{
				mpz_set(temp, a[k][i].s);
				mpz_add_ui(temp, temp, 1);
				mpz_mod_ui(temp, temp, 2);
				term.set(temp, a[k][i].e, a[k][i].m);	
			}
			if(client_type == BOB)
			{
				term = a[k][i];
			}

			for(m = 0; m < col; m++)
			{
				a[k][m]  = f_plus(a[k][m],  f_product(a[i][m],  term, client_type), client_type);
				ia[k][m] = f_plus(ia[k][m], f_product(ia[i][m], term, client_type), client_type);
			}
		}
	}
	
  // upper triangle
	for(i = row-1; i >= 0; i--)
	{
		for(j = i-1; j >= 0; j--)
		{
			if(client_type == ALICE)
			{
				mpz_set(temp, a[j][i].s);
				mpz_add_ui(temp, temp, 1);
				mpz_mod_ui(temp, temp, 2);
				term.set(temp, a[j][i].e, a[j][i].m);
			}
			if(client_type == BOB)
			{
				term = a[j][i];
			}
			for(k = col-1; k >= 0; k--)
			{
				a[j][k]  = f_plus(a[j][k],  f_product(a[i][k],  term, client_type), client_type);
				ia[j][k] = f_plus(ia[j][k], f_product(ia[i][k], term, client_type), client_type);
			}
		}
	}
	mpz_clear(temp);
}

void f_transpose(SFloat **a, int rows, int cols, int client_type,  SFloat **trans)
{
	int i, j;

	for(i = 0; i < rows; i++)
	{
		for(j = 0; j < cols; j++)
			trans[i][j] = a[j][i];
	}
}

void f_matrix_mul(SFloat **a, int a_row, int a_col, SFloat **b, int b_row, int b_col, int client_type, SFloat **product)
{
	int i, j, k;

	for(i = 0; i < a_row; i++)
	{
		for(j = 0; j < b_col; j++)
		{
			product[i][j] = fpShare(0.0, client_type);
			for(k = 0; k < a_col; k++)
				product[i][j] = f_plus(product[i][j], f_product(a[i][k], b[k][j], client_type), client_type);
		}
	}
}

void f_dim1_mul(SFloat **a, int a_row, int a_col, SFloat *b, int b_dim, int client_type, SFloat *product)
{
	int i, k;

	for(i = 0; i < a_row; i++)
	{	
		product[i] = fpShare(0.0, client_type);
		for(k = 0; k < a_col; k++)
      // prudtcu[i]=product[i]+a[i][k]*b[k]
			product[i] = f_plus(product[i], f_product(a[i][k], b[k], client_type), client_type);
	}
}

