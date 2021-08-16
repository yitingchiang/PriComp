#include "debug_util.h"
#include <cstring>
#include <gmp.h>

void print(SFloat s)
{
  double v=sFloat2float(s);
  printf("%f\n",v);
}

double AddSFloat(SFloat f1,SFloat f2)
{
  SFloat v;
  SFloat v2;

  mpz_add(v.m,f1.m,f2.m);
  mpz_mod_ui(v2.m,v.m,DOMAIN_M);
  mpz_add(v.e,f1.e,f2.e);
  mpz_mod_ui(v2.e,v.e,DOMAIN_E);
  mpz_add(v.s,f1.s,f2.s);
  mpz_mod_ui(v2.s,v.s,DOMAIN_S);

  return sFloat2float(v2);
}

void sendSFloat(int socket, SFloat v)
{
  ulong128 data[3];
  data[0]=mpz_get_ui(v.s);
  data[1]=mpz_get_ui(v.m);
  data[2]=mpz_get_ui(v.e);
  char* msg=(char*)(data);
  protocol_write(socket,msg,sizeof(ulong128)*3);
}

int getSFloat(int socket, SFloat* v)
{
  ulong128 data[3];
  char* msg=(char*)data;
  int len=protocol_read(socket,msg,sizeof(ulong128)*3);
  mpz_set_long128(v->s,data[0]);
  mpz_set_long128(v->m,data[1]);
  mpz_set_long128(v->e,data[2]);
  return len;
}

int sMerge(mpz_t v, mpz_t domain, int client_type)
{
  int local_ans=mpz_get_si(v);
  int remote_ans=0;
  int d=mpz_get_ui(domain);

  if(client_type==ALICE)
  {    
    protocol_readint(SP_p_socket,&remote_ans);
    protocol_writeint(SP_p_socket,local_ans);
  }
  else
  {
    protocol_writeint(SP_p_socket,local_ans);
    protocol_readint(SP_p_socket,&remote_ans);
  }
  return (local_ans+remote_ans)%d;
}

double fpMerge(SFloat v1,int client_type)
{
  double v;
  SFloat v2;
  if(client_type==ALICE) 
  {
    getSFloat(SP_p_socket,&v2);
    sendSFloat(SP_p_socket,v1);
  }
  else 
  {
    sendSFloat(SP_p_socket,v1);
    getSFloat(SP_p_socket,&v2);
  }

  return AddSFloat(v1,v2);
}
