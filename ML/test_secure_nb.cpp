#define _WITH_GETLINE

#include "time_eval.h"
#include "ml_protocol/nb.h"
#include <sys/time.h>

#include <cstdio>
#include <ctime>

struct timeval Stamp[2];

void usage(const char* fname)
{
  fprintf(stderr,"Load a separated naive bayes model file and perform prediction.\n\n");
  fprintf(stderr,"usage %s <conf_file> <model_file> <test_file> <client_type>\n",fname);
  fprintf(stderr,"\t<conf_file> is the configuration file of Pricomp. Check the cfg files in the conf directory.");
  fprintf(stderr,"\t<model_file> is the separated model.\n");
  fprintf(stderr,"\t<test_file> is the testing dataset. Both the two parties have this file.\n");
  fprintf(stderr,"\t<client_type>: 1 (Alice) or 2 (Bob).\n\n");
}


int main(int argc, char** argv)
{
  if(argc!=5) { usage(argv[0]); return 1;}

  const char* test_file=argv[3];
  int client_type=(int)strtol(argv[4],NULL,10);
  struct timeval t_main[4];
  char mdlName[256];
  sprintf(mdlName,"%s",argv[2]);
  SecureNB* mdl=(SecureNB*)malloc(sizeof(SecureNB));

  //The following codes create NB model from a saved model file
  SP_init(argv[1], client_type);
  ReadSecureNB(mdlName,mdl); 
  printf("successfully read the secure NB model\n");

  time_set(t_main,0);
  double acc=eval_NB(test_file,mdl,client_type);
  time_set(t_main,1);
  printf("****NB Testing time:****\n");
  time_print(t_main,0,1);
  printf("\n");
  printf("acc=%f\n",acc);

  SP_clear(client_type);

  return 0;
}

