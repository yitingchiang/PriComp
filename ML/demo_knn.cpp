#define _WITH_GETLINE

#include "ml_protocol/knn.h"
#include "time_eval.h"
#include <sys/time.h>

#include <cstdio>
#include <ctime>

void usage(const char* fname)
{
  fprintf(stderr,"Demonstration of using Pricomp to build and test KNN model.\n");
  fprintf(stderr,"Find the detailed implementation under the directory ml_protocol.\n\n");
  fprintf(stderr,"usage: %s <sconf_file> <train_file> <party_dbsize> <local_test_file> <K> <client_type>\n",fname);
  fprintf(stderr,"\t<conf_file> is the configuration file of Pricomp. Check the cfg files in the conf directory.");
  fprintf(stderr,"\t<train_file> is the local training dataset file. Try iris_cluster.alice and iris_cluster.bob in the dataset directory.\n");
  fprintf(stderr,"\t<party_dbsize> is the number of samples hold by the PEER party.\n");
  fprintf(stderr,"\t<local_test_file> is the testing dataset. Both the two parties have this file. Try iris_cluster.test in the dataset directory\n");
  fprintf(stderr,"\t<K> is the number of class. This parameter is given here for convenience.\n");
  fprintf(stderr,"\t<client_type>: 1 (Alice) or 2 (Bob).\n\n");
  //fprintf(stderr,"This program will generate a secure model file with file extension \"sKM\". The two parties have to cooperate to use the model to perform prediction.\n");
}


int main(int argc, char** argv)
{
  if(argc!=7) { usage(argv[0]); return 1;}

  const char* train_file=argv[2];
  const char* test_file=argv[4];
  int party_datasize=(int)strtol(argv[3],NULL,10);
  int client_type=(int)strtol(argv[6],NULL,10);
  int K=(int)strtol(argv[5],NULL,10);

  struct timeval t_main[2];

  SP_init(argv[1], client_type);
  time_set(t_main,0);

  int* predict=secure_KNN(train_file,test_file,K,client_type);
  time_set(t_main,1);

  printf("****KNN time:****\n");
  time_print(t_main,0,1);
  printf("\n");
  delete[] predict;
  SP_clear(client_type);

  return 0;
}
