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
  fprintf(stderr,"KMeans model implementation using PriComp. It is assumed that the dataset is horizontally separated. That is, each party holds some rows (samples). In addition, each party have to know how many samples (rows) the peer party has. After a secure KMeans model being trained using all the samples hold by the two parties, party Alice will read the testing dataset and the two parties will cooperate to perform prediction.\n\n");
  fprintf(stderr,"usage %s <conf_file> <train_file> <party_dbsize> <local_test_file> <K> <client_type>\n",fname);
  fprintf(stderr,"\t<conf_file> is the configuration file of Pricomp. Check the cfg files in the conf directory.");
  fprintf(stderr,"\t<train_file> is the local training dataset file. Try iris_cluster.alice and iris_cluster.bob in the dataset directory.\n");
  fprintf(stderr,"\t<party_dbsize> is the number of samples hold by the PEER party.\n");
  fprintf(stderr,"\t<local_test_file> is the testing dataset. Both the two parties have this file. Try iris_cluster.test in the dataset directory\n");
  fprintf(stderr,"\t<K> is the number of clusters.\n");
  fprintf(stderr,"\t<client_type>: 1 (Alice) or 2 (Bob).\n\n");
  fprintf(stderr,"This program will generate a secure model file with file extension \"sKM\". The two parties have to cooperate to use the model to perform prediction.\n");
}

int main(int argc, char** argv)
{	
  if(argc!=7) { usage(argv[0]); return 1;}

  const char* test_file=argv[4];
  int partyDataSize=(int)strtol(argv[3],NULL,10);
  int client_type=(int)strtol(argv[6],NULL,10);
  int K=(int)strtol(argv[5],NULL,10);
  struct timeval t_main[2];
  char mdlFileName[256];
  sprintf(mdlFileName,"%s.sKM",argv[2]);
  int nIter=10;//default to run 100 iterations

  int nCols=getColSize(argv[2],"\n\t:;, ");
  SP_init(argv[1], client_type);
  time_set(t_main,0);
  Secure_Sample* centers=buildSecureKMeans(argv[2],partyDataSize,K,nIter,client_type);
  time_set(t_main,1);

  for(int i=0; i<K;i++)
  {
    printf("%d-th cluster:\n",i);
    for(int j=0; j<nCols; j++)
    {
      printf("\t%d-th column: %f\n",j,fpMerge(centers[i].getX(j),client_type));
    }
  }

  printf("training time:\n");
  time_print(t_main,0,1);

  printf("nCols=%d\n",nCols);
  writeKmeanModel(centers,K,nCols,mdlFileName);

  int testsize=getRowSize(test_file);
  Secure_Sample *sTestData=NULL;
  //perform prediction
  //only alice reads the testing data
  if(client_type==ALICE)
  {
    Sample *testData=load_all_samples(test_file,testsize,nCols,false);
    sTestData=setSecureData(testData,testsize,nCols);
  }
  else // Bob init null testing data
  {
    sTestData=new Secure_Sample[testsize];
    for(int i=0;i<testsize;i++)
    {
      sTestData[i].setColLen(nCols);
      for(int j=0;j<nCols;j++)
        sTestData[i].setX(float2SFloat(0),j);
    }
  }

  time_set(t_main,0);
  int *result=predict(centers,K,sTestData,testsize,nCols,client_type);
  time_set(t_main,1);
  printf("prediction time:\n");
  time_print(t_main,0,1);
  for(int i=0;i<testsize;i++)
  {
    printf("%d\n",result[i]);
  }
  free(result);    

  SP_clear(client_type);
  return 0;
}
