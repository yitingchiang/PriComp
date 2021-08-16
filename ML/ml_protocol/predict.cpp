#include "predict.h"
#include "debug_util.h"

#include <cstring>

/*
compute Y[] * delta * K[]
*/
SFloat predict_product(int i, int l, SFloat YK, SFloat delta, int client_type)
{
	return f_product(YK, delta, client_type);
}

void init_predict(SFloat* YK_array, SFloat* U, int n1,int n2, SFloat* alpha, SFloat b, int client_type)
{
	int N = n1 + n2;

	SFloat sf_v;
	SFloat sf_tmp1, sf_tmp2;

	for(int l=0;l<N;l++)
	{
		sf_v.set_ui(0, 0, 0);
		for(int m = 0; m < N; m++)
		{	
			
			sf_tmp1 = f_product(YK_array[m*N+l],alpha[m],client_type);
			
			sf_tmp2 = sf_v;
			sf_v = f_plus(sf_tmp1, sf_tmp2, client_type);
		}
		U[l] = f_minus(sf_v, b, client_type);
	}
}

void predict_in_training(SFloat* YK_array, SFloat* U, int n1, int n2, int i, int j, SFloat delta_i, SFloat delta_j, SFloat delta_b,mpz_t is_equal, int client_type)
{
	int N = n1 + n2;

	for(int l=0;l<N;l++)
	{
		SFloat sf_tmp, sf_tmp1, sf_tmp2;
		//U[l]+=(Y[i]*delta_i*K[i*N+l])+(Y[j]*delta_j*K[j*N+l])-delta_b;
		sf_tmp1 = f_product(delta_i,YK_array[i*N+l],client_type);
		
		sf_tmp2 = f_product(delta_j,YK_array[j*N+l],client_type);
		
		sf_tmp  = f_plus(sf_tmp1, sf_tmp2, client_type);
		
		sf_tmp1 = f_minus(sf_tmp, delta_b, client_type);

		sf_tmp2 = U[l];
		
		U[l] = f_if_then_else(is_equal,U[l],f_plus(sf_tmp1, sf_tmp2, client_type),client_type);
	}
}

/*
Scenario: after the model is ready, one of the party uses the shared model to do prediction.
Y: the true label vector of the TRAINING data.
U: the prediction.
Xi: feature vector of the training data.
X: feature vector of the testing data. Hold by the data holder, this is assigned by the last parameter.
alpha and b: SVM parameters.
kernel: the normal kernel function.
client_type: Specify being Alice or Bob.
data_holder: Specify who (Alice or Bob) holds X
if data_holder!=client_type, X=NULL.
That is, only the data holder has the testing data set.
*/
void predict_test(SFloat* Y, SFloat* U, double* Xi, SFloat* alpha, int n1,int n2, double* X, int N2, int nCol, double param, double(*kernel)(double*,double*,double,int), SFloat(*Secure_kernel)(double*,double*,double,int,int), SFloat b,int data_holder,int client_type)
{
	int N=n1+n2;
	SFloat sf_v;
	SFloat sf_tmp1, sf_tmp2;

	for(int l=0;l<N2;l++)
	{
		sf_v.set_ui(0,0,0);
		for(int m = 0; m < N; m++)
		{
			SFloat sf_y;
			SFloat k;

			if(m<n1)//Alice has Xi and Y[m]
			{
				if(client_type==ALICE)
				{
					sf_y=Y[m];

					//alice has X and Xi, and can compute k by herself
					if(client_type==data_holder)
						k=float2SFloat(kernel(Xi+m*nCol,X+l*nCol,param,nCol));
					else//if alice has Xi but she does not hold X
						k=Secure_kernel(Xi+m*nCol,X+(l*nCol),param,nCol,client_type);
				}
				else//Bob does not hold training data Xi in this case.
				{
					//alice does not hold X and cannot compute k by herself, so bob joins.
					if(client_type==data_holder)
						k=Secure_kernel(Xi+m*nCol,X+(l*nCol),param,nCol,client_type);
					//else, Alice can compute k. Bob does nothing.
				}
			}
			else //Bob has Y[m] and Xi
			{
				if(client_type==BOB)
				{
					sf_y=Y[m];
					
					//if Bob has X
					if(client_type==data_holder)
						k=float2SFloat(kernel(Xi+m*nCol,X+l*nCol,param,nCol));
					else//if bob does not hold X
						k=Secure_kernel(Xi+m*nCol,X+(l*nCol),param,nCol,client_type);
				}
				else//Alice's case
				{
					//bob can not compute k by herself
					if(client_type==data_holder)
						k=Secure_kernel(Xi+m*nCol,X+(l*nCol),param,nCol,client_type);
				}	
			}
			
			sf_tmp1 = f_product(sf_y, alpha[m], client_type); 
			
			sf_tmp2 = f_product(sf_tmp1, k, client_type);
			
			sf_tmp1 = sf_v;
			
			sf_v = f_plus(sf_tmp1, sf_tmp2, client_type);		

		}

		U[l] = f_minus(sf_v, b, client_type);
		
	}
}
