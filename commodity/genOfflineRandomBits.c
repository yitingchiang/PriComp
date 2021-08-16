#include<stdio.h>
#include<stdlib.h>
#include<inttypes.h>
#include<sys/types.h>
#include<unistd.h>
#include<string.h>
#include<time.h>

#include "commodity.h"
#define ulong128 unsigned __int128
#define long128 __int128
void usage()
{
  printf("Generate alice.dat and bob.dat for offline commodity mode.\n");
  printf("usage: genData <dimension> <domain>\n");
}

long128 inner_product(ulong128* Ra,ulong128* Rb,int dim,long128 domain)
{
  ulong128 v=0;
  for(int i=0;i<dim;i++) v+=((Ra[i]*Rb[i])%domain);
  return v;
}

int main(int argc, char** argv)
{
  if(argc!=3) { usage(); return 1;}
  srandom(time(NULL));
  ulong128 domain=(ulong128)strtoll(argv[2],NULL,10);
  int dim=(int)strtol(argv[1],NULL,10);

  ulong128* Ra=(ulong128*)malloc(sizeof(long128)*dim);
  ulong128* Rb=(ulong128*)malloc(sizeof(long128)*dim);
  ulong128* ra=(ulong128*)malloc(sizeof(long128)*dim);
  ulong128* rb=(ulong128*)malloc(sizeof(long128)*dim);

  FILE* awp=fopen("alice.dat","wb");
  FILE* bwp=fopen("bob.dat","wb");

  for(int i=0;i<dim;i++)
  {
    Ra[i]=random()%domain;
    Rb[i]=random()%domain;
  }
  printf("\n");

  fwrite(&dim,sizeof(int),1,awp); 
  fwrite(&dim,sizeof(int),1,bwp); 
  fwrite(Ra,sizeof(long128)*dim,1,awp); 
  fwrite(Rb,sizeof(long128)*dim,1,bwp);

  for(int i=0;i<dim;i++)
  {
    if(Ra[i]*Rb[i]!=0)
    {
      ra[i]=random()%(Ra[i]*Rb[i]);
      rb[i]=Ra[i]*Rb[i]-ra[i];
    }
    else
    {
      ra[i]=random()%domain;
      rb[i]=domain-ra[i];
    }

    fwrite(ra+i,sizeof(long128),1,awp); 
    fwrite(rb+i,sizeof(long128),1,bwp);
  }

  fclose(awp); fclose(bwp);    
  free(Ra); free(Rb); free(ra); free(rb);

  return 0;
}
