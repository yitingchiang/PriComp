#define _WITH_GETLINE

#include "time_eval.h"
#include "ml_protocol/nb.h"
#include <sys/time.h>

#include <cstdio>
#include <ctime>

void usage(const char* fname)
{
  fprintf(stderr,"Naive bayes model implementation using PriComp. It is assumed that the dataset is horizontally separated. That is, each party holds some rows (samples). In addition, each party have to know how many samples (rows) the peer party has. A naive bayes model will be trained using all the samples hold by the two parties.\n\n");
  fprintf(stderr,"usage %s <conf_file> <train_file> <party_dbsize> <local_test_file> <client_type>\n",fname);
  fprintf(stderr,"\t<conf_file> is the configuration file of Pricomp. Check the cfg files in the conf directory.");
  fprintf(stderr,"\t<train_file> is the local training dataset file. Try iris.alice and iris.bob in the dataset directory.\n");
  fprintf(stderr,"\t<party_dbsize> is the number of samples hold by the PEER party.\n");
  fprintf(stderr,"\t<local_test_file> is the testing dataset. Both the two parties have this file.\n");
  fprintf(stderr,"\t<client_type>: 1 (Alice) or 2 (Bob).\n\n");
  fprintf(stderr,"This program will generate a secure model file with file extension \"sNB\". The two parties have to cooperate to use the model to perform prediction.\n");
}


int main(int argc, char** argv)
{
  if(argc!=6) { usage(argv[0]); return 1;}

  const char* test_file=argv[4];
  int party_datasize=(int)strtol(argv[3],NULL,10);
  int client_type=(int)strtol(argv[5],NULL,10);
  struct timeval t_main[2];
  char mdlName[256];
  sprintf(mdlName,"%s.sNB",argv[2]);
  SecureNB* mdl=(SecureNB*)malloc(sizeof(SecureNB));


  //Train Secure NB model
  SP_init(argv[1], client_type);
  time_set(t_main,0);
  mdl=TrainSecureNB(argv[2],party_datasize,client_type);
  time_set(t_main,1);

  fprintf(stderr,"****NB Training time:****\n");
  time_print(t_main,0,1);
  fprintf(stderr,"\n");

  //write model
  WriteSecureNB(mdlName,mdl); 

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

