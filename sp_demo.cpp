#define _WITH_GETLINE

#include <cstring>
#include <cstdlib>
#include "src/two_parties/sp_init.h"
#include "src/util.h"
#include "src/communicate.h"
#include "src/debug_util.h"
#include "time_eval.h"

#define DOMAIN 100

struct timeval Stamp[2];

void usage()
{
  fprintf(stderr,"Testing scalar product protocol.\n");
  fprintf(stderr,"usage: client <conf_file> <datafile> <client_type>\n");
  fprintf(stderr,"\t<conf_file> is the configuration file of Pricomp. Check the cfg files in the conf directory.");
  fprintf(stderr,"The format of <datafile>:\nIn the first row is a value N, followed by 2*N rows that the 2*N rows that one row is a value M and another row is M values, separated by space, comma, or tab.\n");
  fprintf(stderr,"\tFor example:\n");
  fprintf(stderr,"\tN\n");
  fprintf(stderr,"\tM_1\n");
  fprintf(stderr,"\tv_[1 1] v_[1 2] ... v_[1 M_1]\n");
  fprintf(stderr,"\tM_2\n");
  fprintf(stderr,"\tv_[2 1] v_[2 2] ... v_[2 M_2]\n");
  fprintf(stderr,"\t......\n");
  fprintf(stderr,"\tM_N\n");
  fprintf(stderr,"\tv_[N 1] v_[N 2] ... v_[N M_N]\n");
  fprintf(stderr,"\t<client_type>: 1 (Alice) or 2 (Bob).\n\n");
}

int main(int argc, char** argv)
{	
  if(argc!=4) { usage(); return 1; }
  int client_type = (int)strtol(argv[3], NULL, 10);
  FILE* rp=fopen(argv[2],"r");
  if(!rp) { printf("can not open file %s!\n",argv[2]); usage(); return 1;}

  char buf[65536];
  fgets(buf,65535,rp);//# of Vectors
  int nVectors=(int)strtol(buf,NULL,10);
  long int **vec=(long int**)malloc(sizeof(long int*)*nVectors);
  int* dims=(int*)malloc(sizeof(int)*nVectors);
  mpz_t* ans=(mpz_t*)malloc(sizeof(mpz_t)*nVectors);;
  int total_dims=0;

  //init
  for(int i=0;i<nVectors;i++)
  {
    mpz_init(ans[i]);

    fgets(buf,65535,rp);//dim
    dims[i]=(int)strtol(buf,NULL,10);
    total_dims+=dims[i];
    vec[i]=(long int*)malloc(sizeof(long int)*dims[i]);
    fgets(buf,65535,rp);//vector data
    char* tmps=strtok(buf,"\t\n, ");
    for(int j=0;tmps;j++)
    {
      vec[i][j]=(long int)strtol(tmps,NULL,10);
      tmps=strtok(NULL,"\t\n, ");
    }
  }

  mpz_t* all_x=(mpz_t*)malloc(sizeof(mpz_t)*total_dims);
  int idx=0;
  for(int i=0;i<nVectors;i++)
  {
    long int *ptr=vec[i];

    for(int j=0;j<dims[i];j++)
    {
      mpz_init(all_x[idx]);
      mpz_set_si(all_x[idx],ptr[j]);
      idx++;
    }
  }

  mpz_t domain;

  mpz_init_set_ui(domain,DOMAIN);

  // Call SP_init before performing any protocols
  SP_init(argv[1], client_type);

  int r;
  int offset=0;
  time_set(Stamp,0);
  for(int j=0;j<nVectors;j++)
  {
    //run scalar product protocol
    r=scalar_product(all_x+offset,client_type, dims[j], domain, ans[0]);
    if(r<0) printf("sp error no=%d\n",r);
    else {
      printf("ans:\n");
      mpz_out_str(stdout,10,ans[0]);
      printf("\n");
    }
    offset+=dims[j];
  }
  time_set(Stamp,1);
  time_print(Stamp,0,1);

  //call this function to perform multiple scalar product at once.
  r=multi_scalar_product(all_x, client_type, nVectors, dims, domain, ans);

  if(r<0) printf("multi sp error! error no=%d\n",r);
  else {
    for(int i=0;i<nVectors;i++)
    {
      printf("ans:\n");
      mpz_out_str(stdout,10,ans[i]);
      printf("\n");
    }
  }

  //free resources
  for(int i=0;i<nVectors;i++) { mpz_clear(ans[i]); free(vec[i]);}
  free(vec);
  free(ans);
  for(int i=0;i<total_dims;i++) mpz_clear(all_x[i]);
  free(all_x);
  mpz_clear(domain);

  //call SP_clear to clear all resource allocated in SP_init
  SP_clear(client_type);

  return 0;
}
