#include <cstdlib>
#include <cinttypes>
#include "gmp.h"
#include <cmath>
#include "../../util.h"
#include "../sp_init.h"

/**
 A vector that all the elements in it are 1
 */
int* dim_1;
/**
 A vector that all the elements in it are 2
 */
int* dim_2;
/*
 Pre-allocated Memory used to perform computations
*/
mpz_t* tmp_mpz_1;//length: 2*bits
/*
 Pre-allocated Memory used to perform computations
*/
mpz_t* tmp_mpz_2;//length: 2*bits
/*
 Pre-allocated Memory used to perform computations
*/
mpz_t* tmp_mpz_3;//length: 2*bits
/*
 Pre-allocated Memory used to perform computations
*/
mpz_t* tmp_mpz_4;//length: 2*bits
extern int adderIDX;

void Kogge_Stone_zn_to_z2(mpz_t,int,mpz_t,int,mpz_t[]);
void Brent_Kung_zn_to_z2(mpz_t,int,mpz_t,int,mpz_t[]);
void ordinary_zn_to_z2(mpz_t,int,mpz_t,int,mpz_t[]);
void multi_zn_to_z2(mpz_t,int,mpz_t,int,mpz_t[]);
/**
 Cyclic rotate of the array a one time
 @param a the array to rotate
 @param dim the length of a
 */
void rotate(mpz_t a[], int dim);
/**
 Cyclic rotate of the array a offset times
 @param offset the number of times to rotate
 @param a the array to rotate
 @param dim the length of a
 */
void lrotate(int offset, mpz_t a[], int dim);
/**
 Set b[indicate % dim] to 1 and others 0
 @param indicate the index of array b to set to 1
 @param dim the length of array b
 @param b the array to be set
 */
void indicate_string(int indicate, int dim, mpz_t b[]);

void (*adder[4])(mpz_t,int,mpz_t,int,mpz_t[])={ordinary_zn_to_z2,multi_zn_to_z2,Brent_Kung_zn_to_z2,Kogge_Stone_zn_to_z2};

/**
 A function to initialize some global memory such that the functions do not need to allocate and free memory each time.
 */
void init_global(int bits)
{
  tmp_mpz_1=(mpz_t*)malloc(sizeof(mpz_t)*bits*16);
  init_vector(tmp_mpz_1,bits*16);
  tmp_mpz_2=tmp_mpz_1+4*bits;
  tmp_mpz_3=tmp_mpz_1+8*bits;
  tmp_mpz_4=tmp_mpz_1+12*bits;

  dim_1=(int*)malloc(sizeof(int)*bits*2);
  dim_2=dim_1+bits;

  for(int i=0;i<bits;i++) 
  {
    dim_1[i]=1;
    dim_2[i]=2;
  }
}

/**
 A function to free global variables.
 For Kogge Stone and Brent Kung adder.
 */
void free_global(int bits)
{
  free(dim_1);
  clear_vector(tmp_mpz_1,bits*16);
  free(tmp_mpz_1);
}

/**
  For Kogge Stone and Brent Kung adder.
  Perform addition between two vectors which have additively shared elements
dd
  @param x The additively shared variables to negate.
  @param y The additively shared variables to negate.
 @param domain The domain of this operation performs on.
 @param len The number of elements in x and y
 @param result The result that result[i] = x[i]+y[i].
 */
void s_multi_product(mpz_t* x, mpz_t* y, mpz_t domain, int nSC, int client_type, mpz_t* result)
{
  mpz_t* tmp=tmp_mpz_1;
  mpz_t* vec=tmp_mpz_2;
  mpz_t* tmp_result=tmp_mpz_3;

  // alice : scalar_product([x, y], domain)
  // bob   : scalar_product([y, x], domain)
  if(client_type == ALICE)
  {
    int idx=0;
    for(int i=0;i<nSC;i++)
    {
      mpz_set(vec[idx++], x[i]);
      mpz_set(vec[idx++], y[i]);
    }
  }
  if(client_type == BOB)
  {
    int idx=0;
    for(int i=0;i<nSC;i++)
    {
      mpz_set(vec[idx++], y[i]);
      mpz_set(vec[idx++], x[i]);
    }
  }

  multi_scalar_product(vec,client_type,nSC,dim_2,domain,tmp);

  // (x * y + tmp) % domain
  int local_idx=0;
  for(int i=0;i<nSC;i++)
  {
    mpz_set_ui(tmp_result[i],0);
    {
      //tmp_result[i]=tmp_result[i]+x[local_idx]+y[local_idx]
      mpz_addmul(tmp_result[i],x[local_idx], y[local_idx]);
      local_idx++;
    }

    mpz_add(tmp_result[i], tmp_result[i], tmp[i]);

    mpz_mod(tmp_result[i], tmp_result[i], domain);
  }

  for(int i=0;i<nSC;i++) mpz_set(result[i],tmp_result[i]);
}

/**
  For Kogge Stone and Brent Kung adder.
  Perform OR between two vectors which have additively shared elements
dd
  @param x The additively shared variables to negate.
  @param y The additively shared variables to negate.
 @param domain The domain of this operation performs on.
 @param len The number of elements in x and y
 @param result The result that result[i] = x[i]+y[i].
 */
void s_multi_or(mpz_t* x, mpz_t* y, mpz_t domain, int len, int client_type, mpz_t* result)
{
  mpz_t* tmp_result=tmp_mpz_4;

  s_multi_product(x,y,domain,len,client_type,tmp_result);
  for(int i=0;i<len;i++)
  {
    mpz_add(tmp_result[i],tmp_result[i],x[i]);
    mpz_add(result[i],tmp_result[i],y[i]);
  }
}

// return x AND y in domain (=2)
void s_multi_and(mpz_t* x, mpz_t* y, mpz_t domain, int len, int client_type, mpz_t* result)
{
  s_multi_product(x,y,domain,len,client_type,result);
}

// return NOT x in domain (=2)
void s_multi_not(mpz_t* x, mpz_t domain, int len, int client_type, mpz_t* result)
{
  // alice : (x + 1) % domain
  // bob   : x
  if(client_type == ALICE)
  {
    for(int i=0;i<len;i++)
    {
      mpz_add_ui(result[i], x[i], 1);
      mpz_mod(result[i], result[i], domain);
    }
  }
  if(client_type == BOB)
  {
    for(int i=0;i<len;i++) mpz_set(result[i], x[i]);
  }
}

//implement the Kogge-Stone adder as a different zn_to_z2 implementation
void Kogge_Stone_zn_to_z2(mpz_t x, int bits, mpz_t domain, int client_type, mpz_t v[])
{
  mpz_t domain2;
  mpz_t* x_bit=(mpz_t*)malloc(sizeof(mpz_t)*(8*bits));
  mpz_t* p_in=x_bit+2*bits;
  mpz_t* tmp_bit=x_bit+4*bits;
  mpz_t* carry=x_bit+6*bits;
  mpz_t* g_in=carry+1;
  mpz_t* p_out=NULL;
  mpz_t* g_out=NULL;
  mpz_t* p2_in=NULL;
  mpz_t* g2_in=NULL;

  init_global(bits);

  init_vector(x_bit,8*bits);
  mpz_init(domain2);
  mpz_set_ui(domain2, 2);

  //Let y_bit be the bits of another parties' share bits. 
  //set p[i]=x_bit[i] XOR y_bit[i] (don't need SMC).
  for(int i=0;i<bits;i++)
  {
    mpz_set_ui(x_bit[i], mpz_tstbit(x, i));
    mpz_set(p_in[i],x_bit[i]);
  }

  mpz_set_ui(carry[0],0);

  multi_scalar_product(x_bit,client_type,bits,dim_1,domain2,g_in);

  p_out=p_in+1;
  g_out=g_in+1;
  p2_in=p_out;
  g2_in=g_out;

  //run Kogge-Stone algorithm to compute the carry bits
  //each gate, compute:
  // (g_in,p_in) \circ (g2_in, p2_in)= (g_in OR (p_in AND g2_in),p_in AND p2_in)
  int step=1;
  int cnt=0;

  while(step<bits)
  {
    cnt++;
    int len=bits-step;

    //compute G_i
    s_multi_and(p2_in,g_in,domain2,len,client_type,tmp_bit);

    s_multi_or(g2_in,tmp_bit,domain2,len,client_type,g_out);

    //compute P_i
    s_multi_and(p_in,p2_in,domain2,len,client_type,p_out);

    step*=2;

    p2_in=p_in+step;
    g2_in=g_in+step;
    p_out=p2_in;
    g_out=g2_in;

  }

  for(int i = 0; i < bits; i++)
  {
    mpz_add(v[i], x_bit[i], carry[i]);
    mpz_mod_ui(v[i], v[i], 2);
  }

  free_global(bits);

  clear_vector(x_bit,8*bits);
  mpz_clear(domain2);

  free(x_bit);
}

//A modified Brent-Kung adder as a different zn_to_z2 implementation.
void Brent_Kung_zn_to_z2(mpz_t x, int bits, mpz_t domain, int client_type, mpz_t v[])
{
  mpz_t* x_bit=(mpz_t*)malloc(sizeof(mpz_t)*bits*10);
  mpz_t* tmp_bit=x_bit+bits;
  mpz_t* carry=x_bit+2*bits;
  mpz_t* p=x_bit+3*bits;
  mpz_t* p_out=x_bit+4*bits;
  mpz_t* g_out=x_bit+5*bits;
  mpz_t* p_in=x_bit+6*bits;
  mpz_t* g_in=x_bit+7*bits;
  mpz_t* p2_in=x_bit+8*bits;
  mpz_t* g2_in=x_bit+9*bits;
  mpz_t domain2;

  init_vector(x_bit,bits*10);
  mpz_init(domain2);
  mpz_set_ui(domain2,2);

  init_global(bits);
  int* IDXTBL=(int*)malloc(sizeof(int)*bits);

  for(int i=0;i<bits;i++)
  {
    mpz_set_ui(x_bit[i], mpz_tstbit(x, i));
    //set p_i
    mpz_set(p[i], x_bit[i]);
  }
  //compute g_i
  multi_scalar_product(x_bit,client_type,bits,dim_1,domain2,carry);
  //the i-th step
  for(int i=0;;i++)
  {
    //j-th bits
    int idx=0;
    int v=1<<i;
    if(v>bits) break;

    for(int j=1;j<bits;j++)
    {
      if((j>>i)%2!=0)//do operation
      {
        int idx2=(j/v)*v-1;

        mpz_set(p_in[idx],p[j]);
        mpz_set(g_in[idx],carry[j]);
        mpz_set(p2_in[idx],p[idx2]);
        mpz_set(g2_in[idx],carry[idx2]);
        IDXTBL[idx]=j;

        idx++;
      }
    }

    if(idx==0) break;

    s_multi_and(p_in,g2_in,domain2,idx,client_type,tmp_bit);

    s_multi_or(g_in,tmp_bit,domain2,idx,client_type,g_out);
    for(int j=0;j<idx;j++) mpz_set(carry[IDXTBL[j]],g_out[j]);

    s_multi_and(p_in,p2_in,domain2,idx,client_type,p_out);
    for(int j=0;j<idx;j++) mpz_set(p[IDXTBL[j]],p_out[j]);
  }

  mpz_set(v[0], x_bit[0]);
  for(int i = 1; i < bits; i++)
  {
    mpz_add(v[i], x_bit[i], carry[i-1]);
    mpz_mod_ui(v[i], v[i], 2);
  }

  free_global(bits);

  clear_vector(x_bit,bits*10);
  mpz_clear(domain2);

  free(IDXTBL);
  free(x_bit);
}

// input x is in domain.
// return v, and v is vector.
void multi_zn_to_z2(mpz_t x, int bits, mpz_t domain, int client_type, mpz_t v[])
{
  int i;
  int dims[2];
  mpz_t tmp[3];
  mpz_t x_bit[2];
  mpz_t vec[12];
  mpz_t domain2;
  mpz_t carry;  //c^{i-1}

  mpz_init(domain2);
  mpz_init(carry);
  init_vector(x_bit,2);
  init_vector(vec, 12);
  init_vector(tmp,3);
  mpz_set_ui(domain2, 2);

  dims[0]=3;
  dims[1]=9;

// the first carry bit = 0. shared y_i=x_i[0]
  mpz_set_ui(x_bit[0], mpz_tstbit(x, 0));
  mpz_set_ui(x_bit[1], mpz_tstbit(x, 1));
  mpz_set(v[0], x_bit[0]);
  mpz_set_ui(carry,0);

  for(i = 1; i+1 < bits; i+=2)
  {
    //x^{i-1}*x^{i}
    mpz_mul(tmp[0],x_bit[0],x_bit[1]);
    //c^{i-1}*x^i
    mpz_mul(tmp[1],carry,x_bit[1]);
    //c^{i-1}*x^{i-1}
    mpz_mul(tmp[2],carry,x_bit[0]);

    if(client_type == ALICE)
    {
      mpz_set(vec[0], x_bit[0]);
      mpz_set(vec[1], carry);
      mpz_set(vec[2], x_bit[0]);

      mpz_set(vec[3], tmp[0]);
      mpz_set(vec[4], tmp[1]);
      mpz_set(vec[5], x_bit[1]);
      mpz_set(vec[6], tmp[2]);
      mpz_set(vec[7], x_bit[0]);
      mpz_set(vec[8], carry);
      mpz_set(vec[9], tmp[0]);
      mpz_set(vec[10], x_bit[0]);
      mpz_set(vec[11], x_bit[1]);
    }
    else
    {
      mpz_set(vec[0], carry);
      mpz_set(vec[1], x_bit[0]);
      mpz_set(vec[2], x_bit[0]);

      mpz_set(vec[3], carry);
      mpz_set(vec[4], x_bit[0]);
      mpz_set(vec[5], tmp[2]);
      mpz_set(vec[6], x_bit[1]);
      mpz_set(vec[7], tmp[1]);
      mpz_set(vec[8], tmp[0]);
      mpz_set(vec[9], x_bit[0]);
      mpz_set(vec[10], tmp[0]);
      mpz_set(vec[11], x_bit[1]);
    }

    for(int j=0;j<2;j++)
      mpz_set_ui(tmp[j],0);

    // alice : carry * x[i] + scalar_product([carry,x[i],x[i]], 2)
    // bob   : carry * x[i] + scalar_product([x[i],carry,x[i]], 2)
    multi_scalar_product(vec, client_type, 2, dims, domain2, tmp);

    //set the i-th bit
    //1. compute the i-th carry bit c^i=the shared result+c^{i-1}x^{i-1}
    mpz_add(carry,tmp[0],tmp[2]);

    //2. compute the shared y^i=c^i+x^i
    mpz_add(v[i],carry,x_bit[1]);
    mpz_mod_ui(v[i], v[i], 2);

    //set the (i+1)-th bit. compute c^{i-1}*x^{i-1}*x^i
    mpz_mul(tmp[2],tmp[2],x_bit[1]);
    //compute c^{i+1}=c^{i-1}*x^{i-1}*x^i+the shared result
    mpz_add(carry,tmp[1],tmp[2]);
    //save carry as the c^{i-1} of the next loop
    //mpz_set(carry,v[i+1]);
    mpz_mod_ui(carry, carry, 2);

    //set the next two bits as x^{i-1} and x^i of the next loop
    mpz_set_ui(x_bit[0], mpz_tstbit(x, i+1));
    mpz_set_ui(x_bit[1], mpz_tstbit(x, i+2));

    //compute y^{i+1}
    mpz_add(v[i+1],carry,x_bit[0]);
    //last: mod
    mpz_mod_ui(v[i+1],v[i+1],2);
  }

  if(i<bits)//so i+1==bits. need to set the last carry bit
  {
    if(client_type == ALICE)
    {
      mpz_set(vec[0], x_bit[0]);
      mpz_set(vec[1], carry);
      mpz_set(vec[2], x_bit[0]);
    }
    else
    {
      mpz_set(vec[0], carry);
      mpz_set(vec[1], x_bit[0]);
      mpz_set(vec[2], x_bit[0]);
    }

    scalar_product(vec, client_type, 3, domain2, v[i]);

    //compute c^{i-1}*x^{i-1}
    mpz_mul(tmp[2],carry,x_bit[0]);
    //compute carry bit
    mpz_add(v[i],v[i],tmp[2]);

    //compute current y bit
    mpz_add(v[i],v[i],x_bit[1]);


    mpz_mod_ui(v[i], v[i], 2);
  }

  mpz_clear(carry);
  mpz_clear(domain2);
  clear_vector(x_bit,2);
  clear_vector(tmp,3);
  clear_vector(vec, 12);
}

//run different versions of zn_to_z2 according to the setting.
void zn_to_z2(mpz_t x, int bits, mpz_t domain, int client_type, mpz_t v[])
{
  adder[adderIDX](x,bits,domain,client_type,v);
}

// input x is in domain.
// return v, and v is vector.
void ordinary_zn_to_z2(mpz_t x, int bits, mpz_t domain, int client_type, mpz_t v[])
{
  int i;
  mpz_t carry;
  mpz_t tmp;
  mpz_t x_bit;
  mpz_t vec[3];
  mpz_t domain2;

  mpz_init(carry); 
  mpz_init(tmp);
  mpz_init(x_bit);
  mpz_init(domain2);
  init_vector(vec, 3);

  mpz_set_ui(carry, 0);
  mpz_set_ui(domain2, 2);

  for(i = 0; i < bits-1; i++)
  {
    // v[i] = (x[i] + carry) % 2
    mpz_set_ui(x_bit, mpz_tstbit(x, i));
    mpz_add(v[i], x_bit, carry);
    mpz_mod_ui(v[i], v[i], 2);

    // alice : carry * x[i] + scalar_product([carry,x[i],x[i]], 2)
    // bob   : carry * x[i] + scalar_product([x[i],carry,x[i]], 2)
    if(client_type == ALICE)
    {
      mpz_set(vec[0], carry);
      mpz_set(vec[1], x_bit);
      mpz_set(vec[2], x_bit);
    }
    if(client_type == BOB)
    {
      mpz_set(vec[0], x_bit);
      mpz_set(vec[1], carry);
      mpz_set(vec[2], x_bit);
    }

    scalar_product(vec, client_type, 3, domain2, tmp);

    mpz_mul(carry, carry, x_bit);
    mpz_add(carry, carry, tmp);

  }

  // v[bits-1] = (x[bits-1] + carry) % 2
  mpz_set_ui(x_bit, mpz_tstbit(x, bits-1));
  mpz_add(v[bits-1], x_bit, carry);
  mpz_mod_ui(v[bits-1], v[bits-1], 2);

  mpz_clear(domain2);
  mpz_clear(carry);
  mpz_clear(tmp);
  mpz_clear(x_bit);
  clear_vector(vec, 3);

}

// input v is vector.
// return a number, it is in domain.
void z2_to_zn(mpz_t v[], int bits, mpz_t domain, int client_type, mpz_t result)
{
  // bits = v.length

  int i;
  mpz_t t;
  mpz_t y;
  mpz_t n;
  mpz_t vec[bits];

  mpz_init(t);
  mpz_init(y);
  mpz_init(n);
  init_vector(vec, bits);

  // alice : vec = v
  // bob   : vec[i] = 2**(i+1) * v[i]
  if(client_type == ALICE)
  {
    for(i = 0; i < bits; i++)
      mpz_set(vec[i], v[i]);
  }
  if(client_type == BOB)
  {
    for(i = 0; i < bits; i++)
      mpz_mul_2exp(vec[i], v[i], i+1);
  }

  scalar_product(vec, client_type, bits, domain, t);

  mpz_set_ui(y, 0);

  // y += v[i] * (2**i) 
  for(i = 0; i < bits; i++)
  {
    mpz_mul_2exp(n, v[i], i);
    mpz_add(y, y, n);
  }

  // (y - t) % domain
  mpz_sub(result, y, t);
  mpz_mod(result, result, domain);

  mpz_clear(t);
  mpz_clear(y);
  mpz_clear(n);
  clear_vector(vec, bits);
}

// return (x + y) in domain    
void s_plus(mpz_t x, mpz_t y, mpz_t domain, mpz_t result)
{
  // (x + y) % domain
  mpz_add(result, x, y);
  mpz_mod(result, result, domain);
}

// return (x - y) in domain
void s_minus(mpz_t x, mpz_t y, mpz_t domain, mpz_t result)
{
  // (x - y) % domain
  mpz_sub(result, x, y);
  mpz_mod(result, result, domain);
}

// return (x * y) in domain
void s_product(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result)
{

  mpz_t tmp;
  mpz_t vec[2];

  mpz_init(tmp);
  init_vector(vec, 2);

  // alice : scalar_product([x, y], domain)
  // bob   : scalar_product([y, x], domain)
  if(client_type == ALICE)
  {
    mpz_set(vec[0], x);
    mpz_set(vec[1], y);
  }
  if(client_type == BOB)
  {
    mpz_set(vec[0], y);
    mpz_set(vec[1], x);
  }

  scalar_product(vec, client_type, 2, domain, tmp);

  // (x * y + tmp) % domain
  mpz_mul(result, x, y);
  mpz_add(result, result, tmp);
  mpz_mod(result, result, domain);

  mpz_clear(tmp);
  clear_vector(vec, 2);

}

// return x**2 in domain
void s_square(mpz_t x, mpz_t domain, int client_type, mpz_t result)
{
  mpz_t tmp;
  mpz_t vec[1];

  mpz_init(tmp);
  mpz_init(vec[0]);

  // scalar_product([x],domain)
  mpz_set(vec[0], x);
  scalar_product(vec, client_type, 1, domain, tmp);

  //(x ** 2 + 2*tmp) % domain
  mpz_mul(result, x, x);
  mpz_mul_ui(tmp, tmp, 2);
  mpz_add(result, result, tmp);
  mpz_mod(result, result, domain);

  mpz_clear(tmp);
  clear_vector(vec, 1);
}

// b: 1 or 0
// set result to x (if b==1) or y (if b==0)
// it is ok if one of x and y is actually result
void s_if_then_else(mpz_t b, mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result)
{
  mpz_t tmp;
  mpz_t n1;
  mpz_t n2;
  mpz_t vec[1];

  mpz_init(tmp);
  mpz_init(n1);
  mpz_init(n2);
  mpz_init(vec[0]);

  // z2_to_zn([b], domain)
  mpz_set(vec[0], b);

  z2_to_zn(vec, 1, domain, client_type, n1);

  // s_product(b, (x - y) % domain, domain)
  mpz_sub(n2, x, y);
  mpz_mod(n2, n2, domain);

  s_product(n1, n2, domain, client_type, tmp);

  // (tmp + y) % domain
  mpz_add(tmp, tmp, y);
  mpz_mod(result, tmp, domain);

  mpz_clear(tmp);
  mpz_clear(n1);
  mpz_clear(n2);
  clear_vector(vec, 1);
}

// b: 1 or 0
// set multiple results to x (if b==1) or y (if b==0)
void s_multiple_if_then_else(mpz_t b, mpz_t* x, mpz_t* y, mpz_t domain, int client_type, mpz_t* result, int nResults)
{
  mpz_t tmp;
  mpz_t n1;
  mpz_t n2;
  mpz_t vec[1];

  mpz_init(tmp);
  mpz_init(n1);
  mpz_init(n2);
  mpz_init(vec[0]);

  // z2_to_zn([b], domain)
  mpz_set(vec[0], b);

  z2_to_zn(vec, 1, domain, client_type, n1);

  for(int i=0;i<nResults;i++)
  {
    // s_product(b, (x - y) % domain, domain)
    mpz_sub(n2, x[i], y[i]);
    mpz_mod(n2, n2, domain);

    s_product(n1, n2, domain, client_type, tmp);

    // (tmp + y) % domain
    mpz_add(tmp, tmp, y[i]);
    mpz_mod(result[i], tmp, domain);
  }
  mpz_clear(tmp);
  mpz_clear(n1);
  mpz_clear(n2);
  clear_vector(vec, 1);
}

/* 
   return (negative?) ? 1 : 0
 */
void s_negative(mpz_t x, mpz_t domain, int client_type, mpz_t result)
{
  int bits = binary_size(domain); 

  mpz_t vec[bits];
  init_vector(vec, bits);

  // zn_to_z2(x, domain).last
  zn_to_z2(x, bits, domain, client_type, vec);

  mpz_set(result, vec[bits - 1]);

  clear_vector(vec, bits);

}


//binary equal output:
//step 1: 
//  each party locally do XOR on their local data
//step 2:
//  one of the party negate its output.
void binary_equal(mpz_t x,mpz_t y,int client_type,mpz_t result)
{
  mpz_t temp;//;,d;
  mpz_init(temp);

  //local XOR
  mpz_add(temp, x, y);
  mpz_mod_ui(temp,temp,2);

  //negative by adding 1.
  if(client_type==ALICE)
  {
    mpz_add_ui(temp,temp,1);
    mpz_mod_ui(result,temp,2);
  }
  else mpz_set(result,temp);

  //set temp result as a flag. this result will be applied if b==0
  mpz_clear(temp);
}

void binary_less_than(mpz_t x,mpz_t y,mpz_t domain,int client_type,mpz_t result)
{
  mpz_t temp,b;//;,d;
  mpz_init(temp);
  mpz_init(b);

  mpz_add(temp, x, y);

  int x0=mpz_get_si(x);
  int y0=mpz_get_si(y);

  //set temp result as a flag. this result will be applied if b==0
  mpz_mod_ui(b,temp,2);

  //set possible result if b==1
  if(x0<y0||(x0==1&&y0==1)) mpz_set_ui(result,1);
  else mpz_set_ui(result,0);
  //if b==0, apply temp. otherwise, apply new temp
  s_if_then_else(b,result,temp,domain,client_type,result);
  mpz_clear(temp);
  mpz_clear(b);
}

// return (x < y) ? 1 : 0
void s_less_than(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result)
{
  //double n;
  mpz_t temp;//;,d;
  mpz_t d;

  if(mpz_cmp_ui(domain,2) == 0)//for binary case only
  {
    binary_less_than(x,y,domain,client_type,result);
    return;
  }

  mpz_init(temp);
  mpz_init(d);
  mpz_set(d,domain);

  mpz_sub(temp, x, y);//temp=x-y

  s_negative(temp, d, client_type, result);

  mpz_clear(temp);
  mpz_clear(d);
}

// return (x > y) ? 1 : 0
void s_great_than(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result)
{
  double n;
  mpz_t temp;
  mpz_t d;

  if(mpz_cmp_ui(domain, 2) == 0)//for binary case
  {
    binary_less_than(y,x,domain,client_type,result);
    return; 
  }

  mpz_init(temp);
  mpz_init(d);
  mpz_set(d,domain);
  // s_negative?(y - x, domain)
  mpz_sub(temp, y, x);//
  s_negative(temp, d, client_type, result);

  mpz_clear(d);
  mpz_clear(temp);
}

// return -x in domain
void s_neg(mpz_t x, mpz_t domain, mpz_t result)
{
  // -x
  mpz_neg(result, x);
}

// return x OR y in domain (=2)
void s_or(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result)
{
  mpz_t temp1;
  mpz_t temp2;
  mpz_init(temp1);
  mpz_init(temp2);

  // tmp1 = s_plus(x, y, domain)
  s_plus(x, y, domain, temp1);
  // tmp2 = s_product(x, y, domain)
  s_product(x, y, domain, client_type, temp2);

  // s_plus(tmp1, tmp2, domain) 
  s_plus(temp1, temp2, domain, result);

  mpz_clear(temp1);
  mpz_clear(temp2);
}

// return x AND y in domain (=2)
void s_and(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result)
{
  // s_product(x, y, domain)
  s_product(x, y, domain, client_type, result);
}

// return NOT x in domain (=2)
void s_not(mpz_t x, mpz_t domain, int client_type, mpz_t result)
{
  // alice adds 1 to x to get (x + 1) % domain
  // bob sets the result to the value he holds
  if(client_type == ALICE)
  {
    mpz_add_ui(result, x, 1);
    mpz_mod(result, result, domain);
  }
  if(client_type == BOB)
  {
    mpz_set(result, x);
  }
}

// return (x == y) ? 1 : 0
void s_equal(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result)
{
  mpz_t great;
  mpz_t less;
  mpz_t domain2;

  if(mpz_cmp_ui(domain, 2) == 0)//for binary case
  {
    binary_equal(x,y,client_type,result);
    return;
  }

  mpz_init(great);
  mpz_init(less);
  mpz_init(domain2);

  mpz_set_ui(domain2, 2);

  s_great_than(x, y, domain, client_type, great);

  s_less_than(x, y, domain, client_type, less);

  //result = s_and(s_not(great, 2), s_not(less, 2), 2)
  s_not(great, domain2, client_type, great);
  s_not(less,  domain2, client_type, less);
  s_and(great, less, domain2, client_type, result);

  mpz_clear(great);
  mpz_clear(less);
  mpz_clear(domain2);
}

// Dividend: x      
// Divisor: y       
// return quotient  
void s_division(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result)
{
  int i;
  int bits = 32;
  int dim = binary_size(domain);

  mpz_t d, u, v, t, s;
  mpz_t vec_x[dim];
  mpz_t vec_y[dim];
  mpz_t q[bits];

  mpz_init(d);
  mpz_init(u);
  mpz_init(v);
  mpz_init(t);
  mpz_init(s);
  init_vector(vec_x, dim);
  init_vector(vec_y, dim);
  init_vector(q, bits);

  // d = 2**(bits*2)
  mpz_set_d(d, pow(2, bits*2));
  // u = z2_to_zn(zn_to_z2(x, domain), d)
  zn_to_z2(x, dim, domain, client_type, vec_x);
  z2_to_zn(vec_x, dim, d, client_type, u);
  // v = z2_to_zn(zn_to_z2(y, domain), d)
  zn_to_z2(y, dim, domain, client_type, vec_y);
  z2_to_zn(vec_y, dim, d, client_type, v);

  for(i = bits - 1; i >= 0; i--)
  {
    // t = (u - v * 2**i) % d
    mpz_mul_2exp(t, v, i);
    mpz_sub(t, u, t);
    mpz_mod (t, t, d);

    // s = s_negative?(t, d)
    s_negative(t, d, client_type, s);

    // alice : q[i] = (1+s) % 2
    // bob   : q[i] = s
    if(client_type == ALICE)
    {
      mpz_add_ui(q[i], s, 1);
      mpz_mod_ui(q[i], q[i], 2);
    }
    if(client_type == BOB)
      mpz_set(q[i], s);

    // u = s_if_then_else(s, u, t, d)
    s_if_then_else(s, u, t, d, client_type, u);
  }

  // z2_to_zn(q, domain)
  z2_to_zn(q, bits, domain, client_type, result);

  mpz_clear(d);
  mpz_clear(u);
  mpz_clear(v);
  mpz_clear(t);
  mpz_clear(s);
  clear_vector(vec_x, dim);
  clear_vector(vec_y, dim);
  clear_vector(q, bits);
}

// Dividend: x
// Divisor: y
// return remainder      
void s_remainder(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result)
{
  int i;
  int bits = 32;
  int dim = binary_size(domain);

  mpz_t d, u, v, t, s;
  mpz_t vec_x[dim];
  mpz_t vec_y[dim];
  mpz_t q[bits];

  mpz_init(d);
  mpz_init(u);
  mpz_init(v);
  mpz_init(t);
  mpz_init(s);
  init_vector(vec_x, dim);
  init_vector(vec_y, dim);
  init_vector(q, bits);

  // d = 2**(bits*2)
  mpz_set_d(d, pow(2, bits*2));
  // u = z2_to_zn(zn_to_z2(x, domain), d)
  zn_to_z2(x, dim, domain, client_type, vec_x);
  z2_to_zn(vec_x, dim, d, client_type, u);
  // v = z2_to_zn(zn_to_z2(y, domain), d)
  zn_to_z2(y, dim, domain, client_type, vec_y);
  z2_to_zn(vec_y, dim, d, client_type, v);

  for(i = bits - 1; i >= 0; i--)
  {
    // t = (u - v * 2**i) % d
    mpz_mul_2exp(t, v, i);
    mpz_sub(t, u, t);
    mpz_mod (t, t, d);

    // s = s_negative?(t, d)
    s_negative(t, d, client_type, s);

    // alice : q[i] = (1+s) % 2
    // bob   : q[i] = s
    if(client_type == ALICE)
    {
      mpz_add_ui(q[i], s, 1);
      mpz_mod_ui(q[i], q[i], 2);
    }
    if(client_type == BOB)
      mpz_set(q[i], s);

    // u = s_if_then_else(s, u, t, d)
    s_if_then_else(s, u, t, d, client_type, u);
  }

  mpz_set(result, u);

  mpz_clear(d);
  mpz_clear(u);
  mpz_clear(v);
  mpz_clear(t);
  mpz_clear(s);
  clear_vector(vec_x, dim);
  clear_vector(vec_y, dim);
  clear_vector(q, bits);
}

// return (base**exp) in domain
void s_power(mpz_t base, mpz_t exp, mpz_t domain, int client_type, mpz_t result)
{
  int i;
  int bits = binary_size(domain);

  mpz_t v, t;
  mpz_t vector[bits];

  mpz_init(v);
  mpz_init(t);
  init_vector(vector, bits);

  // vector = zn_to_z2(exp, domain);
  zn_to_z2(exp, bits, domain, client_type, vector);

  // alice : v = 1
  // bob   : v = 0
  if(client_type == ALICE)
    mpz_set_ui(v, 1);
  if(client_type == BOB)
    mpz_set_ui(v, 0);

  for(i = bits-1; i >= 0; i--)
  {
    // v = s_product(v, v, domain)
    s_product(v, v, domain, client_type, v);
    // t = s_product(v, base, domain)
    s_product(v, base, domain, client_type, t);
    // v = s_if_then_else(vector[i], t, v, domain)
    s_if_then_else(vector[i], t, v, domain, client_type, v);
  }

  // return v
  mpz_set(result, v);

  mpz_clear(v);
  mpz_clear(t);
  clear_vector(vector, bits);

}

// secret: x, n
void s_shift_left(mpz_t x, mpz_t n, mpz_t domain, int client_type, mpz_t result)
{
  mpz_t t;
  mpz_t u;

  mpz_init(t);
  mpz_init(u);

  // alice : t = s_power(2, n, domain)
  // bob   : t = s_power(0, n, domain)
  if(client_type == ALICE)
    mpz_set_ui(u, 2);
  if(client_type == BOB)
    mpz_set_ui(u, 0);
  //t = 2^u
  s_power(u, n, domain, client_type, t);
  
  // s_product(x, t, domain)
  s_product(x, t, domain, client_type, result);

  mpz_clear(t);
  mpz_clear(u);
}

// secret: x, n
void s_shift_right(mpz_t x, mpz_t n, mpz_t domain, int client_type, mpz_t result)
{
  int i;
  int bits = binary_size(domain);

  mpz_t r, t;
  mpz_t vec[bits];
  mpz_t array_x[bits];

  mpz_init(r);
  mpz_init(t);
  init_vector(vec, bits);
  init_vector(array_x, bits);

  // array_x = zn_to_z2(x, domain).reverse
  zn_to_z2(x, bits, domain, client_type, vec);
  for(i = 0; i < bits; i++)
    mpz_set(array_x[bits-1-i], vec[i]);

  // r = z2_to_zn(array_x, domain)
  z2_to_zn(array_x, bits, domain, client_type, r);

  // t = s_shift_left(r, n, domain)
  s_shift_left(r, n, domain, client_type, t);

  // array_x = zn_to_z2(t, domain).reverse
  zn_to_z2(t, bits, domain, client_type, vec);
  for(i = 0; i < bits; i++)
    mpz_set(array_x[bits-1-i], vec[i]);

  // z2_to_zn(array_x, domain)
  z2_to_zn(array_x, bits, domain, client_type, result);

  mpz_clear(r);
  mpz_clear(t);
  clear_vector(vec, bits);
  clear_vector(array_x, bits);
}

void rotate(mpz_t a[], int dim)
{
  int i;
  mpz_t vector[dim];
  init_vector(vector, dim);

  // a.push a.shift
  for(i = 1; i < dim; i++)
    mpz_set(vector[i-1], a[i]);
  mpz_set(vector[dim-1], a[0]);

  for(i = 0; i < dim; i++)
    mpz_set(a[i], vector[i]);

  clear_vector(vector, dim);
}

void lrotate(int offset, mpz_t a[], int dim)
{
  int i;

  for(i = 0; i < offset; i++)
    rotate(a, dim);
}

void indicate_string(int indicate, int dim, mpz_t b[])
{
  // dim = a.length
  int i;

  for(i = 0; i < dim; i++)
    mpz_set_ui(b[i], 0);

  mpz_set_ui(b[indicate % dim], 1);
}

void s_access_array(mpz_t i, int dim, mpz_t array[], mpz_t domain, int client_type, mpz_t item)
{
  int t, n;
  mpz_t tmp;
  mpz_t ary[dim];

  mpz_init(tmp);
  init_vector(ary, dim);

  mpz_set(tmp, i);
  for(t = 0; t < dim; t++)
    mpz_set(ary[t], array[t]);


  // dim = array.length
  //i = i % dim;
  mpz_mod_ui(tmp, i, dim);
  n = mpz_get_ui(tmp);

  // alice : r_array = lrotate(i , array)
  //		   item = scalar_product(r_array, domain)
  // bob   : bit_array = indicate_string(i , array) 
  //         item = scalar_product(bit_array, domain)
  if(client_type == ALICE)
    lrotate(n, ary, dim);
  if(client_type == BOB)
    indicate_string(n, dim, ary);

  scalar_product(ary, client_type, dim, domain, item);

  clear_vector(ary, dim);
  mpz_clear(tmp);

}


