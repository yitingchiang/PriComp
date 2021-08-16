#include "smc.h"
#include "two_parties/protocol/s_float.h"
/*
  to init U.
*/
void init_predict(SFloat* YK, SFloat* U, int n1,int n2, SFloat* alpha, SFloat b, int client_type);

/*
	to update U in SVM training process
*/
void predict_in_training(SFloat* YK_array, SFloat* U, int n1, int n2, int i, int j, SFloat delta_i, SFloat delta_j, SFloat delta_b, mpz_t is_equal,int client_type);

/*
   to predict
*/
void predict_test(SFloat* Y, SFloat* U, double* Xi, SFloat* alpha, int n1,int n2, double* X, int N2, int nCol, double param, double(*kernel)(double*,double*,double,int), SFloat(*Secure_kernel)(double*,double*,double,int,int), SFloat b,int data_holder,int client_type);

