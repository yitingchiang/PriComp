#define _WITH_GETLINE

#include "util.h"
#include "debug_util.h"
#include "time_eval.h"
#include "two_parties/two_parties.h"
#include "two_parties/sp_init.h"
#include "ml_protocol/secure_ml.h"
#include "ml_protocol/kmeans.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>

void usage(const char* fname)
{
  fprintf(stderr,"Run kmeans prediction. A secure k means model file is required.\n\n");
  fprintf(stderr,"usage %s <conf_file> <secure_mdl_file> <test_file> <client_type>\n",fname);
  fprintf(stderr,"\t<conf_file> is the configuration file of Pricomp. Check the cfg files in the conf directory.\n");
  fprintf(stderr,"\t<secure_mdl_file> is the secure model file that are shared between the two parties.\n");
  fprintf(stderr,"\t<local_test_file> is the testing dataset. Both the two parties have this file. Try iris_cluster.test in the dataset directory\n");
  fprintf(stderr,"\t<client_type>: 1 (Alice) or 2 (Bob).\n\n");
}

int main(int argc, char** argv)
{	
  if(argc!=5) { usage(argv[0]); return 1;}

  const char *testfile=argv[3];
  int client_type=(int)strtol(argv[4],NULL,10);
  char mdlfile[256];

  strcpy(mdlfile,argv[2]);

  Secure_Sample* ss_center;
  int testsize=getRowSize(testfile);

  SP_init(argv[1], client_type);

  struct timeval t_main[2];
  time_set(t_main,0);

  int *result=load_and_predict(mdlfile,testfile,client_type);

  time_print(t_main,0,1);
  fprintf(stderr,"\n");

  SP_clear(client_type);

  for(int i=0;i<testsize;i++)
  {
    printf("%d\n",result[i]);
  }
  free(result);

  return 0;
}


