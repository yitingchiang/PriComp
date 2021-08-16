#include<math.h>
#include<stdlib.h>
#include "kernel.h"
#include "debug_util.h"
#include "time_eval.h"


//param is not used in this kernel
double Linear_Kernel(double* X1, double* X2, double param, int len)
{
  double v=0;
  for(int i=0;i<len;i++) v+=X1[i]*X2[i];
  return v;
}

//K(X1,X2)=e^{(-||X1-X2||^2)/(2*delta)}
double RBF(double* X1, double* X2, double delta, int len)
{

  double v=0;
  for(int i=0;i<len;i++)
  {
    double u=(X1[i]-X2[i]);
    v-=(u*u);  
  }
  return exp(delta*v);
}

SFloat Secure_Linear_Kernel(double* X1, double* X2,double param, int len,int client_type)
{
  SFloat result;

  result=f_product(float2SFloat(X1[0]),float2SFloat(X2[0]),client_type);

  for(int i=1;i<len;i++)
  {
    result=f_plus(f_product(float2SFloat(X1[i]),float2SFloat(X2[i]),client_type),result,client_type);
  }

  return result;
}

/*
   K(X1,X2)=e^{(-||X1-X2||^2)*delta}
   Note that X1 is held by alice, and X2 is held by bob.
   The parameter *X is a length-len double array that Alice and Bob hold their own arrays. Initially, these array are not pair-wised shared.
   */
SFloat Secure_RBF(double* X1, double* X2, double delta, int len, int client_type)
{
  SFloat result;//=float2SFloat(0);
  for(int i=0;i<len;i++)
  {
    SFloat f1,f2,u;

    //the reason of sqrt(delta) is to prevent constant*SFloat
    if(client_type==ALICE) f1=float2SFloat(X1[i]*sqrt(delta));
    else f2=float2SFloat(X2[i]*sqrt(delta));

    u=f_minus(f1,f2,client_type);

    //u=(f1-f2)^2
    u=f_product(u,u,client_type);
    //v=v-u
    result=f_minus(result,u,client_type);
  }

  return f_exp(result,client_type);
}

/*
 * Compute matrix K, which K[i*(len1+len2)+j]=Kernel(X[i],X[j]). 
 * X is the vectors which have (len1+len2)*nCol elements.
 * The i-th element of X1 is {X[i*nCol],X[i*nCol+1],...,X[i*nCol+(len-1)]}, for i=0,...,len1-1.
 * The i-th element of X2 is {X[(i+len1)*nCol],X[(i+len1)*nCol+1],...,X[(i+len1)*nCol+(len-1)]}, for i=0,...,len2-1.
 * Secure_kernel is a secure kernel function pointer.
 * kernel is a kernel function that can be computed locally
 * param is the parameter of the kernel function.
 */
SFloat* compute_K(double* X1,int len1,double* X2,int len2, int nCol, SFloat(*Secure_kernel)(double*,double*,double,int,int),double(*kernel)(double*,double*,double,int),double param,int client_type)
{
  int N=len1+len2;
  SFloat* K=new SFloat[N*N];

  if(client_type==ALICE)
  {
    for(int i=0;i<len1;i++) 
    {
      for(int j=i;j<len1;j++)
      {
        K[j*N+i]=K[i*N+j]=float2SFloat(kernel(X1+(i*nCol),X1+(j*nCol),param,nCol));//local scalar product;
      }
    }
  }

  if(client_type==BOB)
  {
    for(int i=0;i<len2;i++) 
    {
      for(int j=i;j<len2;j++)
      {
        K[(j+len1)*N+(i+len1)]=K[(i+len1)*N+(j+len1)]=float2SFloat(kernel(X2+(i*nCol),X2+(j*nCol),param,nCol));//local scalar product;
      }
    }
  }

  for(int i=0;i<len1;i++)
  {
    for(int j=0;j<len2;j++)
    {
      K[i*N+(j+len1)]=K[(j+len1)*N+i]=Secure_kernel(X1+(i*nCol),X2+(j*nCol),param,nCol,client_type);
    }
  }
  return K;
}
