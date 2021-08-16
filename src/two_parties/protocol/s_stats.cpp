#include <gmp.h>

#include "../../smc.h"
#include "s_float.h"
#include "s_stats.h"

SFloat f_summation(SFloat *operands, int dim, int client_type)
{
  int i;
  SFloat sum;

  for(i = 0; i < dim; i++)
    sum = f_plus(sum, operands[i], client_type);
  return sum;
}

SFloat f_average_2(int party_dim, double* local_data, int local_dim, int client_type)
{
  double l_sum=0;
  SFloat sum;
  SFloat Slsum;
  SFloat zero=float2SFloat(0.0);
  SFloat N=float2SFloat((double)party_dim+local_dim);
  SFloat avg;

  for(int i=0;i<local_dim;i++) 
    l_sum+=local_data[i];

  Slsum=float2SFloat(l_sum);

  if(client_type==ALICE)
  {
    sum=f_plus(Slsum,zero,client_type);
    avg=f_division(sum,N,client_type);
  }
  else
  {
    sum=f_plus(zero,Slsum,client_type);
    avg=f_division(sum,zero,client_type);
  }

  return avg;

}

SFloat f_average(SFloat *operands, int dim, int client_type)
{
  SFloat sum, sets, avg;
  sum  = f_summation(operands, dim, client_type);

  if(dim % 2 == 0)
  {
    if(client_type == ALICE)
      sets = fpShare(dim/2, client_type);
    if(client_type == BOB)
      sets = fpShare(dim/2, client_type);
  }
  else
  {
    if(client_type == ALICE)
      sets = fpShare(dim/3 + 1, client_type);
    if(client_type == BOB)
      sets = fpShare(dim/2, client_type);
  }

  avg = f_division(sum, sets, client_type);

  return avg;
}

SFloat f_variance_2(int party_dim, double* local_data, int local_dim, int client_type)
{
  SFloat S_E2X=f_average_2(party_dim,local_data,local_dim,client_type);//secure E^2(X)
  double sqrt_mean=0;
  SFloat S_sqrt_mean;
  SFloat zero=float2SFloat(0.0);
  double N=(double)party_dim+local_dim;
  SFloat S_EX2;//secure E(X^2)
  SFloat var;

  S_E2X=f_product(S_E2X,S_E2X,client_type);//get E^2(X)
  for(int i=0;i<local_dim;i++)
  {
    double x=local_data[i];
    sqrt_mean+=(x*x);
  }
  sqrt_mean/=N;
  S_sqrt_mean=float2SFloat(sqrt_mean);//local sqare sum

  if(client_type==ALICE)
  {
    S_EX2=f_plus(S_sqrt_mean,zero,client_type);//total square mean
  }
  else 
  {
    S_EX2=f_plus(zero,S_sqrt_mean,client_type);
  }

  var=f_minus(S_EX2,S_E2X,client_type);//E(X^2)-E^2(X)
  return var;
}

SFloat f_variance(SFloat *operands, int dim, int client_type)
{
  int i;
  SFloat avg, var, diff, tmp, sets;

  var=float2SFloat(0.0);

  avg = f_average(operands, dim, client_type);
  for(i = 0; i < dim; i++)
  {
    diff = f_minus(operands[i], avg, client_type);
    tmp  = f_product(diff, diff, client_type);
    var  = f_plus(var, tmp, client_type);
  }

  if(client_type == ALICE)
    sets = float2SFloat(dim);
  else if(client_type == BOB)	
    sets = float2SFloat(0.0);

  var = f_division(var, sets, client_type);
  return var;
}

SFloat f_median(SFloat *data1, SFloat* data2, int dim1, int dim2,int client_type)
{
  SFloat tmp;
  SFloat result;
  mpz_t is_less;
  mpz_init(is_less);

  //index of mid. index starts at 0;
  int is_even=((dim1+dim2) % 2 == 0);
  int mid = is_even ? (dim1+dim2)/2 : (dim1+dim2-1)/2;

  printf("mid=%d\n",mid);

  SFloat TWO;
  if(client_type==ALICE) TWO=float2SFloat(2.0);
  SFloat max;
  max.setMax();

  for(int i = 0; i <= mid ; i++)
  {
    f_less_than(data1[0],data2[0],client_type,is_less);//is data1[0]< data2[0]?
    if(is_even)
    {
      //get the last extracted item
      if(i+1==mid) tmp=f_if_then_else(is_less,data1[0],data2[0],client_type);
      else if(i==mid)
      {
        result=f_plus(f_if_then_else(is_less,data1[0],data2[0],client_type),tmp,client_type);
        result=f_division(result,TWO,client_type);
      }
    }
    else if(i==mid) result=f_if_then_else(is_less,data1[0],data2[0],client_type);

    //if data1[0]<data2[0], then move data1[0]<-data1[1], data1[1]<-data1[2],...
    for(int j=0; j<dim1-1;j++)
    {
      data1[j]=f_if_then_else(is_less,data1[j+1],data1[j],client_type);
    }
    data1[dim1-1]=f_if_then_else(is_less,max,data1[dim1-1],client_type);
    //if data1[0]>=data2[0], then move data2[0]<-data2[1], data2[1]<-data2[2],...
    for(int j=0; j<dim2-1;j++)
      data2[j]=f_if_then_else(is_less,data2[j],data2[j+1],client_type);
    data2[dim2-1]=f_if_then_else(is_less,data2[dim2-1],max,client_type);
  }
  mpz_clear(is_less);  

  return result; 
}
