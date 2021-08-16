#include <cstdlib>
#include <cinttypes>
#include <cmath>

#include "s_float.h"
extern mpz_t *logtable;

/**
 Truncate the fractional part of a SFloat value.
 Notice: This function can only get correct results under some specific domains. So don't use this function outside this file.
 @param f The floatint-point value to perfoem round down.
 @param The domain the the result integer value
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result The result additive shared by the two parties.
 */
void f_float2int(SFloat f, mpz_t domain, int client_type, mpz_t result);

/**
 Convert an additive shared integer to SFloat.
 Notice: This function can only get correct results under domains=2^k. So don't use this function outside this file.
 @param mun The additive shared value.
 @param domain The domain of num.
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @return A secure floating-point value.
 */
SFloat f_int2float(mpz_t num, mpz_t domain, int client_type);

/**
 Compute log_2(v) (v is an additive-sharing int) and only get the integer part (additively shared by the two parties)
 Notice: This function can only get correct results under domains=2^k. So don't use this function outside this file.
 @param v The additive shared input value.
 @param logtable The table initialized in SP_init.
 @param domain The domain of v.
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result the result
 */
void s_log2(mpz_t v, mpz_t logtable[], mpz_t domain, int client_type, mpz_t result);
// return 1 if v>0. Otherwise 0
SFloat f_is_positive(SFloat v, int client_type)
{
  mpz_t isPositive;
  mpz_init(isPositive);
  SFloat ZERO=float2SFloat(0.0);
  //Set this value as a constant instead of calling fpShare for speed up.
  SFloat ONE=(client_type==ALICE) ? float2SFloat(0.0): float2SFloat(1.0);
  //SFloat ONE=fpShare(0.5,client_type);//two parties share 0.5+0.5

  f_is_positive(v,client_type,isPositive);
  SFloat result=f_if_then_else(isPositive,ONE,ZERO,client_type);
  mpz_clear(isPositive);
  return result;
}

SFloat float2SFloat(double f)
{
  int i, size, field_m;
  int64_t int_f;
  int e_count;
  int v[120]={0};
  double x;

  mpz_t m;
  mpz_t sign;
  mpz_t mpz_f;
  mpz_t mantissa;
  mpz_t exponent;
  mpz_t offset;
  SFloat s;

  if(f>=0&&f<=0)
  {
    s.set_ui(0,0,0);
    return s;
  }

  mpz_init(m);
  mpz_init(sign);
  mpz_init(mpz_f);
  mpz_init(mantissa);
  mpz_init(offset);
  mpz_init(exponent);


  if(f < 0)
    mpz_set_ui(sign, 1);
  else
    mpz_set_ui(sign, 0);


  f  = (f < 0)? -f:f;
  int_f = floor(f);//the integer part
  mpz_set_ui(mpz_f, int_f);
  x = f - int_f; //the fracional part

  for(i = 0; i < 120 && x > 0; i++)
  {
    x *= 2;
    v[i] = floor(x);
    x -= v[i];
  }

  if(int_f != 0)
  {
    // mantissa
    field_m = MANTISSA_FIELD - 2;
    size = mpz_sizeinbase(mpz_f, 2);

    for(i = size-1; i >= 0; i--)
    {
      field_m--;

      if(mpz_tstbit(mpz_f, i))
        mpz_ui_pow_ui(m, 2, field_m);
      else
        mpz_set_ui(m, 0);
      mpz_add(mantissa, mantissa, m);

      if(!field_m)
        break;
    }

    i = 0;
    while(field_m)
    {
      field_m--;
      if(v[i])
        mpz_ui_pow_ui(m, 2, field_m);
      else
        mpz_set_ui(m, 0);
      mpz_add(mantissa, mantissa, m);
      i++;
    }

    //exponent
    mpz_ui_pow_ui(exponent, 2, EXPONENT_FIELD-1);//exponent=2^{EXPONENT_FIELD-1}
    mpz_sub_ui(exponent, exponent, 1);//exponent=2^{EXPONENT_FIELD-1}-1
    mpz_add_ui(exponent, exponent, size-1);
  }
  else
  {
    e_count = 0;
    for(i = 0; i < 120; i++)
    {
      e_count--;
      if(v[i])
        break;
    }

    field_m = MANTISSA_FIELD - 2;
    for(i = (-e_count)-1; i < MANTISSA_FIELD; i++)
    {
      field_m--;
      if(v[i])
        mpz_ui_pow_ui(m, 2, field_m);
      else
        mpz_set_ui(m, 0);
      mpz_add(mantissa, mantissa, m);

      if(!field_m)
        break;
    }

    mpz_ui_pow_ui(exponent, 2, EXPONENT_FIELD-1);
    mpz_sub_ui(exponent, exponent, 1);
    mpz_sub_ui(exponent, exponent, -e_count);			
  }
  //SFloat s;
  s.set(sign,exponent,mantissa);

  mpz_clear(m);
  mpz_clear(sign);
  mpz_clear(mpz_f);
  mpz_clear(mantissa);
  mpz_clear(offset);
  mpz_clear(exponent);

  return s;

}

double sFloat2float(SFloat f)
{
  int size, size_v;
  int64_t i, iexp;
  double dans;

  mpz_t m;
  mpz_t exp;
  mpz_t v;
  mpz_t ans;

  mpz_init(m);
  mpz_init(exp);
  mpz_init(v);
  mpz_init(ans);

  //if mantissa==0, set result to be 0.
  if(mpz_cmp_ui(f.m,0)==0)
  {
    if(mpz_cmp_ui(f.s,1)==0) return -0.0;
    else return 0.0;
  }

  mpz_ui_pow_ui(exp, 2, EXPONENT_FIELD-1);
  mpz_sub_ui(exp, exp, 1);
  mpz_sub(exp, f.e, exp);

  mpz_set(v, f.m);
  mpz_set_ui(ans, 0);
  iexp = mpz_get_si(exp);
  size_v = mpz_sizeinbase(v, 2);
  dans = 0;

  if(iexp < 0) // exp < 0
  {
    iexp = -iexp;

    size = size_v;
    dans = 0;
    for(i = iexp; i < size_v + iexp; i++)
    {
      size--;
      if(mpz_tstbit(v, size))
        dans += pow(2, -i);
    }
  }
  else if(iexp > MANTISSA_FIELD-3) // exp > Mantissa_field-3
  {
    size = size_v;
    dans = 0;
    for(i = 0; i < size_v; i++)
    {
      size--;
      if(mpz_tstbit(v, size))
        dans += pow(2, -i);
    }
    dans = dans * pow(2, iexp);
  }
  else
  {
    size = size_v;

    for(i = iexp; i >= 0; i--)
    {
      size--;		
      if(mpz_tstbit(v, size))
        mpz_ui_pow_ui(m, 2, i);
      else
        mpz_set_ui(m, 0);
      mpz_add(ans, ans, m);	
    }

    size = size_v;
    dans = mpz_get_d(ans);

    i = 0;
    for(size = size_v-iexp-1-1; size >= 0; size--)
    {
      if(mpz_tstbit(v, size))
        dans += pow(2, -i-1);
      i++;
    }
  }


  if(mpz_cmp_ui(f.s, 1) == 0)
    dans = -dans;

  mpz_clear(m);
  mpz_clear(exp);
  mpz_clear(v);
  mpz_clear(ans);

  return dans;
}


SFloat f_ln(SFloat f, int client_type)
{
  SFloat* v=new SFloat[6];
  SFloat ln2=float2SFloat(0.0);

  mpz_t tmp;
  mpz_init(tmp);
  int offset=(DOMAIN_E/2-1);
  //secure float is (2^{e1+e2})*(m1+m2). 
  SFloat v1;
  SFloat s_offset=float2SFloat(offset);
  SFloat s_exp=float2SFloat(mpz_get_ui(f.e));

  SFloat s_max_exp=float2SFloat(DOMAIN_E-1);
  SFloat s_domain_e=float2SFloat(DOMAIN_E);
  SFloat s_zero=float2SFloat(0.0);
  v1.set(f.s,f.e,f.m);

  //make f.e1+f.e2 < DOMAIN_E 

  if(client_type==ALICE)
  {
    v[5]=float2SFloat(-1.9410661);
    v[4]=float2SFloat(3.5293114);
    v[3]=float2SFloat(-2.4612314);
    v[2]=float2SFloat(1.1306327);
    v[1]=float2SFloat(-.28874220);
    v[0]=float2SFloat(0.31104322e-1);
    //v1.e=offset=2^{exp_field-1}-1
    mpz_set_ui(v1.e,offset/2-1); //v1.e=offset/2-1
    ln2=float2SFloat(0.69314718);//ln(2)

    //if f.e1+f.e2 (sum of exp shared by alice and bob) > DOMAIN_E, then subtract DOMAIN_E form f.e1+f.e2
    //Note: s_great_theni and f_great_than cannot return correct result. NEED TO CHECK!
    s_exp=f_plus(s_exp,s_zero,client_type);
    SFloat exp_sum;
    exp_sum.set(s_exp.s,s_exp.e,s_exp.m);
    exp_sum=f_minus(exp_sum,s_domain_e,client_type);
    f_is_negative(exp_sum,client_type,tmp);
    s_exp=f_if_then_else(tmp,s_exp,exp_sum,client_type);
    s_exp=f_minus(s_exp,s_offset,client_type);
  }
  else
  {
    //s_e=f_plus(s_zero,s_e,client_type);
    for(int i=0;i<=5;i++) v[i]=float2SFloat(0.0);
    mpz_set_ui(v1.e,offset/2); //v1.e=offset/2

    //if f.e1+f.e2 (sum of exp shared by alice and bob) > DOMAIN_E, then subtract DOMAIN_E form f.e1+f.e2
    //Note: s_great_theni and f_great_than cannot return correct result. NEED TO CHECK!
    s_exp=f_plus(s_zero,s_exp,client_type);
    SFloat exp_sum;
    exp_sum.set(s_exp.s,s_exp.e,s_exp.m);
    exp_sum=f_minus(exp_sum,s_zero,client_type);
    f_is_negative(exp_sum,client_type,tmp);
    s_exp=f_if_then_else(tmp,s_exp,exp_sum,client_type);
    s_exp=f_minus(s_exp,s_zero,client_type);

  }

  //compute ln(x) = -1.9410661+(3.5293114+(-2.4612314+(1.1306327+(-.28874220+0.31104322e-1*x)*x)*x)*x)*x
  for(int i=1;i<=5;i++)
  {
    v[0]=f_product(v[0],v1,client_type);
    v[0]=f_plus(v[i],v[0],client_type);
  }

  v1=f_product(s_exp,ln2,client_type);//(e1+e2)*ln2

  v1=f_plus(v1,v[0],client_type);//(e1+e2)*ln2+ln(x)

  mpz_clear(tmp);

  return v1;
}


SFloat f_plus(SFloat f1, SFloat f2, int client_type)
{
  int i;
  int bits_e, bits_m;
  mpz_t e1, e2;
  mpz_t domain_e, domain_m;
  mpz_t de;
  mpz_t cond;
  mpz_t tmp;
  mpz_t m1, m2;
  mpz_t n1, n2;
  mpz_t m, n;
  mpz_t left, right;
  mpz_t sign, exponent, mantissa;
  mpz_t el, lf, dis, ml, er, mr;

  mpz_t zero, domain_2;
  mpz_init(zero);
  mpz_init(domain_2);
  mpz_set_ui(domain_2, 2);
  mpz_set_ui(zero, 0);


  mpz_init(exponent);
  mpz_init(mantissa);
  mpz_init(sign);

  mpz_init(left);
  mpz_init(right);
  mpz_init(el);
  mpz_init(lf);
  mpz_init(dis);
  mpz_init(ml);
  mpz_init(er);
  mpz_init(mr);

  mpz_init(e1);
  mpz_init(e2);
  mpz_init(domain_e);
  mpz_init(domain_m);
  mpz_init(de);
  mpz_init(cond);
  mpz_init(tmp);
  mpz_init(m1);
  mpz_init(m2);
  mpz_init(n1);
  mpz_init(n2);
  mpz_init(m);
  mpz_init(n);

  mpz_set_ui(domain_e, DOMAIN_E);
  mpz_set_ui(domain_m, DOMAIN_M);

  /*
     fprintf(stderr,"f1:%f\n",fpMerge(f1,client_type));
     fprintf(stderr,"f2:%f\n\n",fpMerge(f2,client_type));
     fprintf(stderr,"f1 detail:\n");
     f1.to_s();
     fprintf(stderr,"f2 detail:\n");
     f2.to_s();
   */

  bits_e = EXPONENT_FIELD; //binary_size(domain_e);
  bits_m = MANTISSA_FIELD; //binary_size(domain_m);

  mpz_t vec[bits_e];
  mpz_t m_array[bits_m];

  init_vector(vec, bits_e);
  init_vector(m_array, bits_m);

  zn_to_z2(f1.e, bits_e, domain_e, client_type, vec);
  z2_to_zn(vec , bits_e, domain_m, client_type, e1);

  init_vector(vec, bits_e);
  zn_to_z2(f2.e, bits_e, domain_e, client_type, vec);
  z2_to_zn(vec , bits_e, domain_m, client_type, e2 );

  // de = (e1 - e2) % SFloat.domain_m
  mpz_sub(de, e1, e2);
  mpz_mod(de, de, domain_m);

  // cond = s_negative(de, SFloat.domain_m);
  s_negative(de, domain_m, client_type, cond);

  mpz_neg(de, de);
  s_shift_right(f1.m, de, domain_m, client_type, tmp);
  s_if_then_else(cond, tmp, f1.m, domain_m, client_type, m1);
  mpz_neg(de, de);
  s_shift_right(f2.m, de, domain_m, client_type, tmp);
  s_if_then_else(cond, f2.m, tmp, domain_m, client_type, m2);

  mpz_neg(n1, m1);
  mpz_neg(n2, m2);

  s_if_then_else(f1.s, n1, m1, domain_m, client_type, m1);
  s_if_then_else(f2.s, n2, m2, domain_m, client_type, m2);

  mpz_add(m, m1, m2);
  mpz_mod(m, m, domain_m);
  mpz_neg(n, m);

  s_negative(m, domain_m, client_type, sign);
  s_if_then_else(sign, n, m, domain_m, client_type, m);
  mpz_neg(n, m);

  zn_to_z2(m, bits_m, domain_m, client_type, m_array);

  if(client_type == ALICE)
  {
    mpz_add_ui(m1, m_array[bits_m-3], 1);
    mpz_mod_ui(m1, m1, 2);

    mpz_add_ui(m2, m_array[bits_m-2], 1);
    mpz_mod_ui(m2, m2, 2);
  }
  if(client_type == BOB)
  {
    mpz_set(m1, m_array[bits_m-3]);

    mpz_set(m2, m_array[bits_m-2]);
  }

  s_product(m1, m2, domain_2, client_type, left);
  mpz_set(right, m_array[bits_m-2]);

  s_negative(m, domain_m, client_type, m1);
  s_negative(n, domain_m, client_type, m2);

  mpz_add(zero, m1, m2);
  mpz_mod_ui(zero, zero, 2);

  s_product(left, zero, domain_2, client_type, left);
  s_if_then_else(cond, f2.e, f1.e, domain_e, client_type, exponent);

  mpz_set(el, exponent);
  if(client_type == ALICE)
  {
    mpz_add_ui(lf, left, 1);
    mpz_mod_ui(lf, lf, 2);
  }
  if(client_type == BOB)
  {
    mpz_set(lf, left);
  }

  mpz_set_ui(dis, 0);
  for(i = bits_m-3; i >= 0; i--)
  {
    s_product(lf, m_array[i], domain_2, client_type, tmp);
    mpz_add(lf, lf, m_array[i]);
    mpz_mod_ui(lf, lf, 2);
    mpz_sub(lf, lf, tmp);
    mpz_mod_ui(lf, lf, 2);

    if(client_type == ALICE)
    {	
      mpz_add_ui(tmp, dis, 1);
      s_if_then_else(lf, dis, tmp, domain_m, client_type, dis);
    }
    if(client_type == BOB)
    {
      s_if_then_else(lf, dis, dis, domain_m, client_type, dis);	
    }
  }

  mpz_sub(tmp, el, dis);
  s_if_then_else(left, tmp, el, domain_e, client_type, el);
  s_shift_left(m, dis, domain_m, client_type, ml);

  mpz_set(er, exponent);
  if(client_type == ALICE)
  {	
    mpz_set_ui(tmp, 1);
    s_shift_right(m, tmp, domain_m, client_type, tmp);
    s_if_then_else(right, tmp, m, domain_m, client_type, mr);

    mpz_add_ui(tmp, er, 1);
    s_if_then_else(right, tmp, er, domain_e, client_type, er);

  }
  if(client_type == BOB)
  {
    mpz_set_ui(tmp, 0);
    s_shift_right(m, tmp, domain_m, client_type, tmp);
    s_if_then_else(right, tmp, m, domain_m, client_type, mr);

    s_if_then_else(right, er, er, domain_e, client_type, er);
  }


  s_if_then_else(right, mr, ml, domain_m, client_type, mantissa);

  s_if_then_else(right, er, el, domain_e, client_type, exponent);

  //mantissa<0?
  s_negative(mantissa,domain_m,client_type,n1);
  mpz_neg(m,mantissa);
  //mantissa>0?
  s_negative(m,domain_m,client_type,n2);
  //mantissa!=0
  s_or(n1,n2,domain_2,client_type,m);
  mpz_set_ui(zero,0);
  //if mantissa=0, exponent=0.
  s_if_then_else(m,exponent,zero,domain_e,client_type,exponent);

  SFloat s;
  s.set(sign,exponent,mantissa);

  mpz_clear(e1);
  mpz_clear(e2);
  mpz_clear(domain_e);
  mpz_clear(domain_m);
  mpz_clear(domain_2);
  mpz_clear(de);
  mpz_clear(cond);
  mpz_clear(tmp);
  mpz_clear(m1);
  mpz_clear(m2);
  mpz_clear(n1);
  mpz_clear(n2);
  mpz_clear(m);
  mpz_clear(n);
  mpz_clear(left);
  mpz_clear(right);
  mpz_clear(zero);
  mpz_clear(sign);
  mpz_clear(exponent);
  mpz_clear(mantissa);
  mpz_clear(el);
  mpz_clear(lf);
  mpz_clear(dis);
  mpz_clear(ml);
  mpz_clear(er);
  mpz_clear(mr);
  clear_vector(vec, bits_e);
  clear_vector(m_array, bits_m);
  return s;
}

// the two parties shares the summation of their input
SFloat fpShare(double raw, int client_type)
{
  SFloat s1 = float2SFloat(raw);
  SFloat s2; //zero

  //Alice and Bob get shared (s1+s2)
  if(client_type == ALICE)
    return f_plus(s1, s2, client_type);
  if(client_type == BOB)
    return f_plus(s2, s1, client_type);
}

//return f1-f2
SFloat f_minus(SFloat f1, SFloat f2, int client_type)
{
  SFloat result;
  if(client_type == ALICE)
  {
    result=f_plus(f1, f2, client_type);	
  }
  if(client_type == BOB)
  {
    SFloat s;//=f2;
    s.set(f2.s,f2.e,f2.m);
    mpz_add_ui(s.s, s.s, 1);
    mpz_mod_ui(s.s, s.s, 2);
    result=f_plus(f1, s, client_type);
  }

  return result;
}

/*
   compute x=f1*f2, where f1 and f2 are both shared sfloat.
return: the shared sfloat x.
 */
SFloat f_product(SFloat f1, SFloat f2, int client_type)
{
  int i, count;
  int bits;
  int bits_m;
  int bits_m2;
  int num;

  mpz_t m, m1, m2;
  mpz_t sign;
  mpz_t exponent;
  mpz_t mantissa1, mantissa2, mantissa;
  mpz_t tmp;
  mpz_t domain;
  mpz_t domain_m, domain_e;
  mpz_t msb;
  mpz_t vec[MANTISSA_FIELD];

  mpz_t zero,domain_2;
  mpz_init(zero);
  mpz_init(domain_2);
  mpz_set_ui(zero,0);
  mpz_set_ui(domain_2,2);


  mpz_init(m);
  mpz_init(m1);
  mpz_init(m2);
  mpz_init(sign);
  mpz_init(exponent);
  mpz_init(tmp);
  mpz_init(domain);
  mpz_init(domain_m);
  mpz_init(domain_e);
  mpz_init(msb);
  mpz_init(mantissa1);
  mpz_init(mantissa2);
  mpz_init(mantissa);
  init_vector(vec, MANTISSA_FIELD);


  mpz_set_ui(domain_e, DOMAIN_E);

  bits = MANTISSA_FIELD-2;

  mpz_add(tmp, f1.s, f2.s);
  mpz_mod_ui(sign, tmp, DOMAIN_S);

  if(client_type == ALICE)
  {
    mpz_ui_pow_ui(tmp, 2, EXPONENT_FIELD-1);
    mpz_sub_ui(tmp, tmp, 1);
    mpz_add(exponent, f1.e, f2.e);
    mpz_sub(exponent, exponent, tmp);
    mpz_mod_ui(exponent, exponent, DOMAIN_E);
  }
  if(client_type == BOB)
  {
    mpz_add(exponent, f1.e, f2.e);
    mpz_mod_ui(exponent, exponent, DOMAIN_E);
  }

  mpz_set_ui(domain, DOMAIN_M);
  mpz_pow_ui(domain, domain, 2);

  mpz_set_ui(domain_m, DOMAIN_M);
  bits_m = MANTISSA_FIELD;
  zn_to_z2(f1.m, bits_m, domain_m, client_type, vec);
  z2_to_zn(vec,  bits_m, domain,   client_type, m1);

  init_vector(vec, MANTISSA_FIELD);
  zn_to_z2(f2.m, bits_m, domain_m, client_type, vec);

  z2_to_zn(vec,  bits_m, domain  , client_type, m2);

  s_product(m1, m2, domain, client_type, m);

  bits_m2 = MANTISSA_FIELD * 2;
  mpz_t m_array[bits_m2];
  init_vector(m_array, bits_m2);

  zn_to_z2(m, bits_m2, domain, client_type, m_array);
  mpz_set(msb, m_array[bits*2-1]);

  mpz_t v[bits];
  init_vector(v, bits);

  count = bits;
  for(i = 0; i < bits; i++)
    mpz_set(v[i], m_array[count++]);

  z2_to_zn(v, bits, domain_m, client_type, mantissa1);

  count = bits-1;
  for(i = 0; i < bits; i++)
    mpz_set(v[i], m_array[count++]);

  z2_to_zn(v, bits, domain_m, client_type, mantissa2);

  s_if_then_else(msb, mantissa1, mantissa2, domain_m, client_type, mantissa);

  if(client_type == ALICE)
  {	
    mpz_add_ui(tmp, exponent, 1);
    s_if_then_else(msb, tmp, exponent, domain_e, client_type, exponent);
  }
  if(client_type == BOB)
  {
    s_if_then_else(msb, exponent, exponent, domain_e, client_type, exponent);	
  }

  //mantissa<0?
  s_negative(mantissa,domain_m,client_type,m1);
  mpz_neg(m,mantissa);
  //mantissa>0?
  s_negative(m,domain_m,client_type,m2);
  //mantissa!=0
  s_or(m1,m2,domain_2,client_type,m);
  //if mantissa=0, exponent=0.
  s_if_then_else(m,exponent,zero,domain_e,client_type,exponent);
  mpz_clear(zero);
  mpz_clear(domain_2);

  SFloat s;
  s.set(sign, exponent, mantissa);

  mpz_clear(m);
  mpz_clear(m1);
  mpz_clear(m2);
  mpz_clear(sign);
  mpz_clear(exponent);
  mpz_clear(tmp);
  mpz_clear(domain);
  mpz_clear(domain_m);
  mpz_clear(domain_e);
  mpz_clear(msb);
  mpz_clear(mantissa1);
  mpz_clear(mantissa2);
  mpz_clear(mantissa);

  clear_vector(vec, MANTISSA_FIELD);
  clear_vector(m_array, bits_m2);
  clear_vector(v, bits);
  return s;
}

/*
   Compute d=f1/f2. return sfloat d.
 */
SFloat f_division(SFloat f1, SFloat f2, int client_type)
{
  int64_t l;
  int i;
  int b;
  int count;

  int bits;
  int bits_m;

  mpz_t zero,domain_2;
  mpz_init(zero);
  mpz_init(domain_2);
  mpz_set_ui(zero,0);
  mpz_set_ui(domain_2,2);


  mpz_t m1, m2;
  mpz_t t, s;
  mpz_t msb;
  mpz_t mantissa;
  mpz_t mantissa1;
  mpz_t mantissa2;
  mpz_t tmp;
  mpz_t domain;
  mpz_t domain_m;
  mpz_t domain_e;
  mpz_t domain_binary;
  mpz_t vec_m[MANTISSA_FIELD];
  mpz_t sign;
  mpz_t exponent;

  mpz_init(m1);
  mpz_init(m2);
  mpz_init(t);
  mpz_init(s);
  mpz_init(msb);
  mpz_init(mantissa);
  mpz_init(mantissa1);
  mpz_init(mantissa2);
  mpz_init(tmp);
  mpz_init(domain);
  mpz_init(domain_m);
  mpz_init(domain_e);
  mpz_init(domain_binary);
  mpz_init(sign);
  mpz_init(exponent);


  init_vector(vec_m, MANTISSA_FIELD);

  bits = MANTISSA_FIELD - 2;
  bits_m = MANTISSA_FIELD;

  mpz_set_ui(domain, DOMAIN_M);
  mpz_pow_ui(domain, domain, 2);

  mpz_set_ui(domain_m, DOMAIN_M);
  mpz_set_ui(domain_binary, 2);

  b = bits + MANTISSA_FIELD;

  mpz_t v[b];
  mpz_t v2[bits];

  init_vector(v, b);
  init_vector(v2, bits);

  zn_to_z2(f1.m, bits_m, domain_m, client_type, vec_m);

  for(i = 0; i < bits; i++)
    mpz_set_ui(v[i], 0);
  for(i = bits; i < b; i++)
    mpz_set(v[i], vec_m[i-bits]);

  z2_to_zn(v, b, domain, client_type, m1);

  init_vector(vec_m, bits_m);

  zn_to_z2(f2.m, bits_m, domain_m, client_type, vec_m);
  z2_to_zn(vec_m, bits_m, domain, client_type, m2);

  mpz_t q_array[bits + 1];
  init_vector(q_array, bits+1);

  for(i = bits; i >= 0; i--)
  {
    mpz_mul_2exp(tmp, m2, i);
    mpz_sub(t, m1, tmp);
    mpz_mod(t, t, domain);
    s_negative(t, domain, client_type, s);

    if(client_type == ALICE)
    {
      mpz_add_ui(q_array[i], s, 1);
      mpz_mod_ui(q_array[i], q_array[i], 2);
    }
    if(client_type == BOB)
    {
      mpz_set(q_array[i], s);
    }

    s_if_then_else(s, m1, t, domain, client_type, m1);	
  }

  mpz_set(msb, q_array[bits]);
  count = 1;

  for(i = 0; i < bits; i++)
    mpz_set(v2[i], q_array[count++]);
  z2_to_zn(v2, bits, domain_m, client_type, mantissa1);

  count = 0;
  for(i = 0; i < bits; i++)
    mpz_set(v2[i], q_array[count++]);

  z2_to_zn(v2, bits, domain_m, client_type, mantissa2);

  s_if_then_else(msb, mantissa1, mantissa2, domain_m, client_type, mantissa);

  mpz_add(sign, f1.s, f2.s);
  mpz_mod_ui(sign, sign, DOMAIN_S);

  mpz_set_ui(domain_e, DOMAIN_E);

  if(client_type == ALICE)
  {
    mpz_ui_pow_ui(tmp, 2, EXPONENT_FIELD-1);
    mpz_sub_ui(tmp, tmp, 1);
    mpz_sub(exponent, f1.e, f2.e);
    mpz_add(exponent, exponent, tmp);
    mpz_mod_ui(exponent, exponent, DOMAIN_E);

    mpz_sub_ui(tmp, exponent, 1);
    s_if_then_else(msb, exponent, tmp, domain_e, client_type, exponent);	
  }
  if(client_type == BOB)
  {
    mpz_sub(exponent, f1.e, f2.e);
    mpz_mod_ui(exponent, exponent, DOMAIN_E);

    s_if_then_else(msb, exponent, exponent, domain_e, client_type, exponent);	
  }

  //mantissa<0?
  s_negative(mantissa,domain_m,client_type,m1);
  mpz_neg(tmp,mantissa);
  //mantissa>0?
  s_negative(tmp,domain_m,client_type,m2);
  //mantissa!=0
  s_or(m1,m2,domain_2,client_type,tmp);
  //if mantissa=0, exponent=0.
  s_if_then_else(tmp,exponent,zero,domain_e,client_type,exponent);
  mpz_clear(zero);
  mpz_clear(domain_2);


  SFloat sfloat;
  sfloat.set(sign, exponent, mantissa);

  mpz_clear(m1);
  mpz_clear(m2);
  mpz_clear(t);
  mpz_clear(s);
  mpz_clear(msb);
  mpz_clear(mantissa);
  mpz_clear(mantissa1);
  mpz_clear(mantissa2);
  mpz_clear(tmp);
  mpz_clear(domain);
  mpz_clear(domain_m);
  mpz_clear(domain_e);
  mpz_clear(domain_binary);
  mpz_clear(sign);
  mpz_clear(exponent);

  clear_vector(vec_m, bits_m);

  clear_vector(v, b);
  clear_vector(v2, bits);
  clear_vector(q_array, bits+1);

  return sfloat;
}


// return 1 if v>0. Otherwise 0
SFloat f_relu(SFloat v, int client_type)
{
  mpz_t isPositive;
  mpz_init(isPositive);
  SFloat ZERO=float2SFloat(0.0);

  f_is_positive(v,client_type,isPositive);
  SFloat result=f_if_then_else(isPositive,v,ZERO,client_type);
  mpz_clear(isPositive);
  return result;
}

void f_equal(SFloat f1,SFloat f2,int client_type,mpz_t result)
{
  mpz_t s_eqflag, exp_eqflag, m_eqflag;
  mpz_t domain_s, domain_e, domain_m, bin_domain;
  mpz_t tmp;

  mpz_init(s_eqflag);
  mpz_init(m_eqflag);
  mpz_init(exp_eqflag);
  mpz_init(domain_s);
  mpz_init(domain_e);
  mpz_init(domain_m);
  mpz_init(bin_domain);
  mpz_init(tmp);

  mpz_set_ui(domain_s,DOMAIN_S);
  mpz_set_ui(domain_e,DOMAIN_E);
  mpz_set_ui(domain_m,DOMAIN_M);
  mpz_set_ui(bin_domain,2);

  s_equal(f1.s,f2.s,domain_s,client_type,s_eqflag);
  s_equal(f1.e,f2.e,domain_e,client_type,exp_eqflag);
  s_equal(f1.m,f2.m,domain_m,client_type,m_eqflag);
  //f1.s==f2.s and f1.e==f2.e?
  s_and(s_eqflag,exp_eqflag,bin_domain,client_type,tmp);
  //f1.s==f2.s and f1.e==f2.e and f1.m==f2.m?
  s_and(m_eqflag,tmp,bin_domain,client_type,result);
  mpz_clear(tmp);
  mpz_clear(domain_s);
  mpz_clear(domain_e);
  mpz_clear(domain_m);
  mpz_clear(bin_domain);
  mpz_clear(s_eqflag);
  mpz_clear(m_eqflag);
  mpz_clear(exp_eqflag);
}

//only check if mantissa is 0. This implementation is different with IEEE floating-point values because we have one more bit to explicitly represent the integer bit 1 on mantissa.  
void f_is_zero(SFloat f,int client_type,mpz_t result)
{
  mpz_t mantissa;
  mpz_t n1,n2,m1;
  mpz_t domain_m, domain_2;
  mpz_t ZERO;
 
  mpz_init(mantissa);
  mpz_init(m1);
  mpz_init(n1);
  mpz_init(n2);
  mpz_init(domain_m);
  mpz_init(domain_2);
  mpz_init(ZERO);
  
  mpz_set_ui(domain_m,DOMAIN_M);
  mpz_set_ui(domain_2,2);
  mpz_set(mantissa,f.m);
  mpz_set_ui(ZERO,0);

  //n1=1 if mantissa<0
  s_negative(mantissa,domain_m,client_type,n1);
  mpz_neg(m1,mantissa); // m1 = -mantissa

  //n2=1 if mantissa>0
  s_negative(m1,domain_m,client_type,n2);

  //m1=1 if mantissa!=0 (n1=1 or n2=1)
  s_or(n1,n2,domain_2,client_type,m1);

  mpz_set(result,m1);
  if(client_type==ALICE)
    mpz_binnegate(result,m1);//result = (m1+1)%2
 
  //If result is true (1), make f.s=0. 
  //To make -0 becomes +0
  s_if_then_else(result,ZERO,f.s,domain_2,client_type,f.s);

  mpz_clear(ZERO);
  mpz_clear(mantissa);
  mpz_clear(m1);
  mpz_clear(n1);
  mpz_clear(n2);
  mpz_clear(domain_m);
  mpz_clear(domain_2);
}

//Note that this function returns 0 even if f = +0
void f_is_positive(SFloat f, int client_type, mpz_t result)
{
  mpz_t tmp;
  mpz_t not_zero,is_zero;
  mpz_t domain_2;
  mpz_init(tmp);
  mpz_init(not_zero);
  mpz_init(is_zero);
  mpz_init(domain_2);

  mpz_set_ui(domain_2,2);

  //first flag: f!=0?
  //Note that if f is zero, f_is_zero will set -0 to +0
  f_is_zero(f,client_type,is_zero);
  mpz_set(not_zero,is_zero);
  if(client_type==ALICE)
    mpz_binnegate(not_zero,is_zero);

  //second flag: f.s>=0?
  //get sign bit. if shared 0+0 or 1+1, after alice negates the result, we get 1.
  mpz_set(tmp,f.s);
  if(client_type==ALICE)
    mpz_binnegate(tmp,f.s);

  //f is positive if f.s > 0 and f is not zero
  s_and(tmp,not_zero,domain_2,client_type,result);

  mpz_clear(domain_2);
  mpz_clear(is_zero);
  mpz_clear(not_zero);
  mpz_clear(tmp);
}

//Note that this function returns 0 even if f = -0
void f_is_negative(SFloat f,int client_type,mpz_t result)
{
  f_is_positive(f,client_type,result);
  if(client_type==ALICE) 
  {
    mpz_add_ui(result,result,1);
    mpz_mod_ui(result,result,2);
  }
  
  //mpz_set(result,f.s);
}

void f_great_than(SFloat f1, SFloat f2, int client_type,  mpz_t result)
{
  f_less_than(f2,f1,client_type,result);
}

void f_less_than(SFloat f1, SFloat f2, int client_type, mpz_t result)
{
  
  mpz_t sign_flag, exponent_flag, mantissa_flag;
  mpz_t s_eqflag, exp_eqflag,m_eqflag;
  mpz_t domain_s, domain_e, domain_m, domain;
  mpz_t tmp,isZERO;

  mpz_init(sign_flag);
  mpz_init(s_eqflag);
  mpz_init(exponent_flag);
  mpz_init(exp_eqflag);
  mpz_init(m_eqflag);
  mpz_init(mantissa_flag);
  mpz_init(tmp);

  mpz_init(domain_s);
  mpz_init(domain_e);
  mpz_init(domain_m);
  mpz_init(isZERO);

  mpz_set_ui(domain_s,DOMAIN_S);//2^1 = 2
  mpz_set_ui(domain_e,DOMAIN_E);//2^11 = 2048
  mpz_set_ui(domain_m,DOMAIN_M);//2^55

  mpz_init(domain);
  mpz_set_ui(domain, 2);

  SFloat backup_f1=f1;
  SFloat backup_f2=f2;

  //check if f1 is zero
  f_is_zero(f1,client_type,isZERO);
  if(client_type==ALICE) mpz_add_ui(tmp,f1.e,DOMAIN_E/2-1);//exponent bias

  s_if_then_else(isZERO,tmp,f1.e,domain_e,client_type,f1.e);

  f_is_zero(f2,client_type,isZERO);
  if(client_type==ALICE) mpz_add_ui(tmp,f2.e,DOMAIN_E/2-1);//exponent bias

  s_if_then_else(isZERO,tmp,f2.e,domain_e,client_type,f2.e);


  //Step 1: check sign bit
  //check if f1.s > f2.s (which means f1 < f2)
  //sign_flag = 1 if f1 is negative and f2 is positive
  s_great_than(f1.s, f2.s, domain_s, client_type, sign_flag);
  
  //if f1.s==f2.s?
  s_less_than(f1.s, f2.s, domain_s, client_type, s_eqflag);//check f1.s < f2.s
  s_not(s_eqflag,domain,client_type,s_eqflag);//f1.s!=f2.s
  s_not(sign_flag,domain,client_type,tmp);//f1.s<=f2.s
  s_and(tmp,s_eqflag,domain,client_type,s_eqflag); //f1.s<=f2,s and f1.s != f2.s
  //post condition:
  //sign_flag = 1 if f1 is negative and f2 is positive (f1 < f2)
  //s_eqflag = 1 if f1.s = f2.s
 
  //Step 2: check exponent
  //f1.e < f2.e ?
  s_less_than(f1.e,f2.e,domain_e,client_type,exponent_flag);
  //f1.e > f2.e ?
  s_great_than(f1.e, f2.e, domain_e, client_type, exp_eqflag);
  //decide if e1=e2
  s_not(exponent_flag,domain,client_type,tmp);//e1<=e2 ?
  s_not(exp_eqflag,domain,client_type,exp_eqflag);//e1>=e2 ?
  s_and(tmp,exp_eqflag,domain,client_type,exp_eqflag); // (e1<=e2 and e1>=e2) ?
  //post condition:
  //if f1.e < f2.e, exponent = 1
  //if f1.e = f2.e, exp_eqflag = 1

  //Step 3: check mantissa
  //f1.m < f2.m ?
  s_less_than(f1.m, f2.m, domain_m, client_type, mantissa_flag); 
  //f1.m > f2.m ?
  s_great_than(f1.m, f2.m, domain_m, client_type, m_eqflag);
  s_not(mantissa_flag,domain,client_type,tmp);//f1.m >= f2.m ?
  s_not(m_eqflag,domain,client_type,m_eqflag);// f1.m<=f2.m?
  //if f1.m==f2.m?
  s_and(tmp,m_eqflag,domain,client_type,m_eqflag);	
  //post condition:
  //if f1.m < f2.m, mantissa_flag = 1
  //if f1.m = f2.m, m_eqflag =1

 
  //To decide if f1<f2, check if any of the three cases is satisfied:
  //1. f1.s > f2.s (f1 is negative and f2 is possitive)
  //2. f1.s = f2.s and f1.e < f2.e
  //3. f1.s = f2.s and f1.e = f2.e and f1.m < f2.m

  //check if s1=s2 and f1.e < f2.e?
  s_and(s_eqflag,exponent_flag,domain,client_type,tmp);

  //store the result of case 1 or case 2 to result
  s_or(sign_flag,tmp,domain,client_type,result);

  //check if s1=s2 and e1=e2 and f1.m < f2.m
  //s1=s2 and e1=e2?
  s_and(s_eqflag,exp_eqflag,domain,client_type,tmp);

  //check if s1==s2 and e1==e2 and f1.m<f2.m
  s_and(tmp,mantissa_flag,domain,client_type,tmp);//step 2

  //check if f1.s>f2.s or (f1.s=f2.s and f1.e<f2.e) or (f1.s=f2.s and f1.e=f2.e and f1.m<f2.m)
  s_or(result, tmp, domain, client_type, result);
  
   
  //the last step: result XOR (f1.s and f2.s). Since if the two values are both negative, the result should be negated. 
  s_and(f1.s,f2.s,domain,client_type,tmp);

  mpz_add(tmp,tmp,result);

  mpz_clear(sign_flag);
  mpz_clear(exponent_flag);
  mpz_clear(s_eqflag);
  mpz_clear(exp_eqflag);
  mpz_clear(m_eqflag);
  mpz_clear(mantissa_flag);

  mpz_clear(tmp);
  mpz_clear(domain_s);
  mpz_clear(domain_e);
  mpz_clear(domain_m);

  mpz_clear(domain);
  f1=backup_f1;
  f2=backup_f2;
}

SFloat f_if_then_else(mpz_t b, SFloat f1, SFloat f2, int client_type)
{
  mpz_t domain_s, domain_e, domain_m;

  mpz_init(domain_s);
  mpz_init(domain_e);
  mpz_init(domain_m);

  mpz_set_ui(domain_s,DOMAIN_S);
  mpz_set_ui(domain_e,DOMAIN_E);
  mpz_set_ui(domain_m,DOMAIN_M);

  SFloat result;

  s_if_then_else(b, f1.s, f2.s, domain_s, client_type, result.s);
  s_if_then_else(b, f1.e, f2.e, domain_e, client_type, result.e);
  s_if_then_else(b, f1.m, f2.m, domain_m, client_type, result.m);

  mpz_clear(domain_s);
  mpz_clear(domain_e);
  mpz_clear(domain_m);

  return result;
}

SFloat f_exp(SFloat f1, int client_type)
{
  // set_no = 0
  int i, j;
  mpz_t positive_d;
  mpz_t is_positive; // 0: f1 <=0. 1: f1>0.

  mpz_t domain; // 0: f1 <=0. 1: f1>0.

  SFloat zero, one;
  SFloat x, y;
  SFloat inversed_y; // compute 1/y.
  SFloat temp, tmp;
  SFloat f_int;

  mpz_t int_part;

  mpz_init(is_positive);
  mpz_init(positive_d);
  mpz_init(int_part);

  mpz_init(domain);
  //domain: 2^64
  //mpz_set_str(domain,"18446744073709551616",10);
  mpz_set_ui(domain, 4294967296); //2^32

  //TING-YU gen logtable
  SFloat negated_f1;
  negated_f1.set(f1.s,f1.e,f1.m);

  if(client_type==ALICE)
    negated_f1.negate();

  mpz_set(is_positive, negated_f1.s);

  //if the value is a nagetive number, make it positive
  //x =f1;
  x = f_if_then_else(is_positive,f1,negated_f1,client_type);

  //f1*log2(e), by TY
  //log2(e) = 1.442695040888960
  if(client_type == ALICE)	tmp = float2SFloat(1.442695040888960);
  if(client_type == BOB)		tmp = float2SFloat(0);
  x = f_product(x, tmp, client_type);

  f_float2int(x,domain,client_type,int_part);

  f_int = f_int2float(int_part, domain, client_type);

  x = f_minus(x, f_int,client_type);

  //0.693147180559945
  if(client_type == ALICE)	tmp = float2SFloat(0.693147180559945);
  if(client_type == BOB)		tmp = float2SFloat(0);
  x = f_product(x, tmp, client_type);

  //Get int and feaction parts of x. Start main procedure.............	
  int order=6;
  int v=6+order*4;
  SFloat sv=(client_type==ALICE)? float2SFloat((float)v):float2SFloat(0);
  //z=x^2
  SFloat z=f_product(x,x,client_type);
  SFloat d;
  SFloat n;
  SFloat tmpd;

  d=z;
  n=sv;

  for(int i=order;i>0;i--)
  {
    v-=4;
    sv=(client_type==ALICE)? float2SFloat((float)v): float2SFloat(0);
    //tmpd=v*d+n;
    tmpd=f_plus(f_product(sv,d,client_type),n,client_type);
    //n=z*d;
    n=f_product(z,d,client_type);
    d=tmpd;
  }
  //compute: 
  //1. d=2*d-d*X+n;
  //2. n=2*d+d*X+n;
  SFloat d2;
  d2=d;
  if(client_type==ALICE)
    mpz_add_ui(d2.e,d2.e,1);//compute d2=2*d

  temp=f_plus(d2,n,client_type);//2*d+n
  tmp=f_product(d,x,client_type);//d*X

  n=f_plus(temp,tmp,client_type);
  d=f_minus(temp,tmp,client_type);

  y=f_division(n,d,client_type);
  //end main procedure

  mpz_add(y.e, y.e, int_part);

  one=float2SFloat(1.0);
  zero.set_ui(0,0,0);

  if(client_type==ALICE) 
    inversed_y=f_division(one,y,client_type);
  else 
    inversed_y=f_division(zero,y,client_type);

  temp=f_if_then_else(is_positive,y,inversed_y,client_type);

  mpz_clear(int_part);
  mpz_clear(positive_d);
  mpz_clear(is_positive);

  return temp;
}


void f_float2int(SFloat f, mpz_t domain, int client_type, mpz_t result)
{
  int i;
  int n, m;
  mpz_t bias;
  mpz_t exp;
  mpz_t domain_e, domain_m;
  mpz_t k[MANTISSA_FIELD];
  mpz_t v[MANTISSA_FIELD-2];
  mpz_t bs[MANTISSA_FIELD-2];
  mpz_t cl[1];
  mpz_t base;
  mpz_t n1;
  mpz_t t;
  mpz_t q;
  mpz_t d;
  mpz_t tmp;
  mpz_t m_int;
  mpz_t n_int;

  mpz_init(bias);
  mpz_init(exp);
  mpz_init(domain_e);
  mpz_init(domain_m);
  mpz_init(base);
  mpz_init(n1);
  mpz_init(t);
  mpz_init(q);
  mpz_init(d);
  mpz_init(tmp);
  mpz_init(m_int);
  mpz_init(n_int);

  init_vector(k,  MANTISSA_FIELD);
  init_vector(v,  MANTISSA_FIELD-2);
  init_vector(bs, MANTISSA_FIELD-2);
  init_vector(cl, 1);

  mpz_set_ui(domain_e, DOMAIN_E);
  mpz_set_ui(domain_m, DOMAIN_M);

  if(client_type == ALICE)
  {
    mpz_ui_pow_ui(bias, 2, EXPONENT_FIELD-1);
    mpz_sub_ui(bias, bias, 1);
  }
  if(client_type == BOB)
  {
    mpz_set_ui(bias, 0);
  }

  s_minus(f.e, bias, domain_e, exp);

  zn_to_z2(f.m, MANTISSA_FIELD, domain_m, client_type, k);
  n = MANTISSA_FIELD;
  m = MANTISSA_FIELD-2; // k[0..n-3] k[n-2]

  for(i = 0; i < m; i++)
    mpz_set(v[i], k[i]);

  mpz_set_ui(m_int, 0);

  if(client_type == ALICE)
  {
    mpz_set_ui(base, 1);
    mpz_set_ui(n1, 1);
  }
  if(client_type == BOB)
  {
    mpz_set_ui(base, 1);
    mpz_set_ui(n1, 0);
  }
  mpz_set_ui(t, 0);
  mpz_set_ui(q, 0);

  for(i = m-1; i >= 0; i--)
  {
    mpz_set(cl[0], v[i]);
    s_negative(exp, domain_e, client_type, d);
    mpz_set(bs[i], d);
    z2_to_zn(cl, 1, domain, client_type, q);
    s_product(base, m_int, domain, client_type, tmp);
    s_plus(tmp, q, domain, t);
    s_if_then_else(d, m_int, t, domain, client_type, m_int);
    s_minus(exp, n1, domain_e, exp);
  }

  mpz_neg(n_int, m_int);
  s_if_then_else(f.s, n_int, m_int, domain, client_type, result);

  mpz_clear(bias);
  mpz_clear(exp);
  mpz_clear(domain_e);
  mpz_clear(domain_m);
  mpz_clear(base);
  mpz_clear(n1);
  mpz_clear(t);
  mpz_clear(q);
  mpz_clear(d);
  mpz_clear(tmp);
  mpz_clear(m_int);
  mpz_clear(n_int);

  clear_vector(k,  MANTISSA_FIELD);
  clear_vector(v,  MANTISSA_FIELD-2);
  clear_vector(bs, MANTISSA_FIELD-2);
  clear_vector(cl, 1);
}

SFloat f_int2float(mpz_t num, mpz_t domain, int client_type)
{
  int i;
  int bits = binary_size(domain); 

  mpz_t sign, exponent, mantissa;
  mpz_init(sign);
  mpz_init(exponent);
  mpz_init(mantissa);

  mpz_t d, on, bias, lg, expe, disp, n, tmp;
  mpz_init(d);
  mpz_init(on);
  mpz_init(bias);
  mpz_init(lg);
  mpz_init(expe);
  mpz_init(disp);
  mpz_init(tmp);
  mpz_init(n);

  mpz_t zero, domain_2, domain_s, domain_e, domain_m;
  mpz_init_set_ui(zero, 0);
  mpz_init_set_ui(domain_2, 2);
  mpz_init_set_ui(domain_s, DOMAIN_S);
  mpz_init_set_ui(domain_e, DOMAIN_E);
  mpz_init_set_ui(domain_m, DOMAIN_M);

  mpz_t v[bits];
  mpz_t v2[MANTISSA_FIELD-2];
  init_vector(v, bits);
  init_vector(v2, MANTISSA_FIELD-2);

  //ting-yu 
  mpz_set(tmp, domain);
  mpz_mul_ui(tmp, tmp, 2);
  mpz_mul_ui(tmp, tmp, 2);

  s_negative(num, tmp, client_type, d);

  if(client_type == ALICE)
    mpz_set_ui(on, 1);
  if(client_type == BOB)
    mpz_set_ui(on, 0);//0717(2->0)
  s_if_then_else(d, on, zero, domain_s, client_type, sign);
  if(client_type == ALICE)
  {
    mpz_set_ui(tmp, 1);
    mpz_mul_2exp(bias, tmp, EXPONENT_FIELD-1);
    mpz_sub_ui(bias, bias, 1);
  }
  if(client_type == BOB)
  {
    mpz_set_ui(bias, 0);
  }

  s_log2(num, logtable, domain, client_type, lg);

  mpz_set(expe, lg);
  zn_to_z2(lg, bits, domain, client_type, v);
  z2_to_zn(v,  bits, domain_e, client_type, lg);
  s_plus(lg, bias, domain_e, exponent);

  ///*
  if(client_type == ALICE)
  {
    mpz_set_ui(disp, bits-1);
    mpz_sub(disp, disp, expe);
  }
  if(client_type == BOB)
  {
    mpz_sub(disp, zero, expe);
  }
  s_shift_left(num, disp, domain, client_type, tmp);

  zn_to_z2(tmp, bits, domain, client_type, v);
  for(i = 0; i < MANTISSA_FIELD-2; i++)
    mpz_set_ui(v2[i], 0);
  int count = MANTISSA_FIELD-2;
  for(i = bits-1; i >= 0; i--)
    mpz_set(v2[--count], v[i]);
  //*/

  z2_to_zn(v2, MANTISSA_FIELD-2, domain_m, client_type, mantissa);

  mpz_t n1, n2, m;
  mpz_init(n1);
  mpz_init(n2);
  mpz_init(m);

  //mantissa<0?
  s_negative(mantissa,domain_m,client_type,n1);
  mpz_neg(m,mantissa);
  //mantissa>0?
  s_negative(m,domain_m,client_type,n2);
  //mantissa!=0
  s_or(n1,n2,domain_2,client_type,m);
  //if mantissa=0, exponent=0.
  s_if_then_else(m,exponent,zero,domain_e,client_type,exponent);
  mpz_clear(n1);
  mpz_clear(n2);
  mpz_clear(m);

  SFloat s;
  s.set(sign, exponent, mantissa);

  mpz_clear(sign);
  mpz_clear(exponent);
  mpz_clear(mantissa);

  mpz_clear(d);
  mpz_clear(on);
  mpz_clear(bias);
  mpz_clear(lg);
  mpz_clear(expe);
  mpz_clear(disp);
  mpz_clear(tmp);
  mpz_clear(n);

  mpz_clear(zero);
  mpz_clear(domain_s);
  mpz_clear(domain_e);
  mpz_clear(domain_m);

  clear_vector(v, bits);
  clear_vector(v2, MANTISSA_FIELD-2);

  return s;
}

void s_log2(mpz_t v, mpz_t logtable[], mpz_t domain, int client_type, mpz_t result)
{
  int i, j, n1;
  int parts;
  int p;
  int bits = binary_size(domain);

  parts = 8 - 1;
  p = parts + 1;

  mpz_t temp[p];
  mpz_t n, c;
  mpz_t table_value;
  mpz_t r_tmp, d;
  mpz_t zero;
  mpz_t r;

  mpz_init(n);
  mpz_init(c);
  mpz_init(d);
  mpz_init(r);
  mpz_init(table_value);
  mpz_init(r_tmp);
  mpz_init_set_ui(zero, 0);
  init_vector(temp, parts+1);

  mpz_set_ui(n, 0);
  for(i = 0; i < (parts+1); i++)
  {
    if(client_type == ALICE)
      mpz_set_ui(n, (parts-i) * 8);

    s_shift_right(v, n, domain, client_type, temp[i]);
  }	

  mpz_set_ui(r, 0);

  //i = 4;
  i = 8;

  while(i > 0)
  {

    if(client_type == ALICE)
      mpz_set_ui(c, (8-i)*8);
    if(client_type == BOB)
      mpz_set_ui(c, 0);

    if(bits < 8)	n1 = mpz_get_ui(domain);
    else			n1 = 256;

    s_access_array(temp[i-1], n1, logtable, domain, client_type, table_value);

    s_plus(c, table_value, domain, r_tmp);
    s_equal(temp[i-1], zero, domain, client_type, d);	
    s_if_then_else(d, r, r_tmp, domain, client_type, result);
    mpz_set(r, result);

    i = i-1;
  }

  //mpz_set(result, r);

  mpz_clear(n);
  mpz_clear(c);
  mpz_clear(d);
  mpz_clear(r);
  mpz_clear(table_value);
  mpz_clear(r_tmp);
  clear_vector(temp, parts+1);
}

