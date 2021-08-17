/**
 This header file defines secure matrix operations. 
 */
#include "../cSFloat.h"
#ifndef S_MATRIX_H
#define S_MATRIX_H
/**
*  Use Gauss-Jordan Elimination to compute matrix inverse
*  @param a the matrix to compute the inverse
*  @param row number of rows
*  @param col number of columns
*  @param client_type The type of this party:
*     - ALICE (1)
*     - BOB (2)
*  @param ia The inverse of matrix a
*/
void f_matrix_inverse(SFloat **a, int row, int col, int client_type, SFloat **ia);

/**
  Compute the transpose of matrix
  @param a A the matrix to transpose
  @param rows number of rows
  @param cols number of columns
  @param client_type The type of this party:
     - ALICE (1)
     - BOB (2)
   @param trans The result
*/
void f_transpose(SFloat **a, int rows, int cols, int client_type,  SFloat **trans);
/**
* Compute matrix multiplication
* @param a a a_row x a_col matrix
* @param a_row row size of matrix a
* @param a_col col size of matrix a
* @param b a b_row x b_col matrix
* @param b_row row size of matrix b
* @param b_col col size of matrix b
* @param client_type The type of this party:
*     - ALICE (1)
*     - BOB (2)
* @param product the result of a*b
*/
void f_matrix_mul(SFloat **a, int a_row, int a_col, SFloat **b, int b_row, int b_col, int client_type, SFloat **product);

/**
* Compute the multiplication of a*b where a is an a_row x a_col matrix and b is an b_dim x 1 matrix
* @param a a a_row x a_col matrix
* @param a_row row size of matrix a
* @param a_col col size of matrix a
* @param b a b_dim x 1 matrix (a b_dim-dimensional vector)
* @param b_dim dimension of b
* @param client_type The type of this party:
*     - ALICE (1)
*     - BOB (2)
* @param product the result of a*b
*/
void f_dim1_mul(SFloat **a, int a_row, int a_col, SFloat *b, int b_dim, int client_type, SFloat *product);
#endif
