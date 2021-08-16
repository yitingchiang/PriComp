#include "cSFloat.h"

SFloat::SFloat(int iinfinity = 0, int inan = 0)
{
	mpz_init_set_ui(s, 0);
	mpz_init_set_ui(e, 0);
	mpz_init_set_ui(m, 0);
	mpz_init_set_ui(infinity, iinfinity); //true: 1 , false: 0
	mpz_init_set_ui(NaN, inan);           //true: 1 , false: 0
	return;
}

SFloat::SFloat()
{
	mpz_init_set_ui(s, 0);
	mpz_init_set_ui(e, 0);
	mpz_init_set_ui(m, 0);
	mpz_init_set_ui(infinity, 0); 
	mpz_init_set_ui(NaN, 0);      
	return;
}

SFloat::SFloat(const SFloat &s)
{
	mpz_init_set(this->s, s.s);
	mpz_init_set(this->e, s.e);
	mpz_init_set(this->m, s.m);
	mpz_init_set(this->infinity, s.infinity); 
	mpz_init_set(this->NaN, s.NaN);      
	return;
}

SFloat& SFloat::operator=(const SFloat& v)
{
	if(this==&v) return *this;
        
	mpz_set(s,v.s);
	mpz_set(e,v.e);
	mpz_set(m,v.m);
	mpz_set(infinity,v.infinity);
	mpz_set(NaN,v.NaN);
	return *this;
}

void SFloat::setMax()
{
	mpz_set_ui(s, 0);
	mpz_set_ui(e, DOMAIN_E/2-1);
	mpz_set_ui(m, DOMAIN_M);
	mpz_init_set_ui(infinity, 0); //true: 1 , false: 0
	mpz_init_set_ui(NaN, 0);      //true: 1 , false: 0
}

void SFloat::set_ui(unsigned long int is, unsigned long int ie, unsigned long int im)
{
	mpz_set_ui(s, is);
	mpz_set_ui(e, ie);
	mpz_set_ui(m, im);
	mpz_init_set_ui(infinity, 0); //true: 1 , false: 0
	mpz_init_set_ui(NaN, 0);           //true: 1 , false: 0
}

void SFloat::set(mpz_t is, mpz_t ie, mpz_t im)
{
	mpz_set(s, is);
	mpz_set(e, ie);
	mpz_set(m, im);
	mpz_init_set_ui(infinity, 0); //true: 1 , false: 0
	mpz_init_set_ui(NaN, 0);           //true: 1 , false: 0
}

/**
 */
void SFloat::to_s()
{
	fprintf(stdout,"s: ");
	mpz_out_str(stdout, 10, s);
	fprintf(stdout,",e: ");
	mpz_out_str(stdout, 2, e);
	fprintf(stdout,",m: ");
	mpz_out_str(stdout, 2, m);
	fprintf(stdout,"\n");
}

void SFloat::negate()
{
	if(mpz_cmp_ui(s,0)==0) mpz_set_ui(s,1);
	else if(mpz_cmp_ui(s,1)==0) mpz_set_ui(s,0);
	else 
	{ 
		fprintf(stderr,"ERROR on sign bit:\n");
		mpz_out_str(stderr,10,s);
		fprintf(stderr,"\n");
	}
}

void SFloat::writeRaw(FILE* fp)
{
  mpz_out_raw(fp,infinity);
  mpz_out_raw(fp,NaN);
  mpz_out_raw(fp,s);
  mpz_out_raw(fp,e);
  mpz_out_raw(fp,m);
}

void SFloat::readRaw(FILE* fp)
{
  mpz_inp_raw(infinity,fp);
  mpz_inp_raw(NaN,fp);
  mpz_inp_raw(s,fp);
  mpz_inp_raw(e,fp);
  mpz_inp_raw(m,fp);
}

