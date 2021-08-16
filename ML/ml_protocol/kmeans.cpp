#include "kmeans.h"
#include <cfloat>

//evenly divide each dimension in to K parts and set the center of each parts as each init centers' value of that dimension
Secure_Sample* init_KMeans_Centers(Sample* data,int K, int nCols, int nRows,Secure_Sample* ss_center)
{
  for(int i=0;i<nCols;i++)
  {
    double max=DBL_MIN;
    double min=DBL_MAX;
    //compute the max and min value of the i-th column (feature)
    for(int j=0;j<nRows;j++)
    {
      max = data[j].X[i] > max ? data[j].X[i] : max;
      min = data[j].X[i] < min ? data[j].X[i] : min;
    }
    double width=(max-min)/K;//divided into K parts
    double v=min+width/2;//starts from min+w/2

    //set the center of the k intervals as the value of the k centers's i-th dimension
    for(int j=0;j<K;j++)
    {
      ss_center[j].setX(float2SFloat(v),i);
      v+=width;
    }
  }
  return ss_center;
}
  
Secure_Sample* setCenter(Sample* data, int K, int nCols, int nRows, int client_type)
{
   //set centers
  Secure_Sample *ss_center=new Secure_Sample[K];
  for(int i=0;i<K;i++) ss_center[i].setColLen(nCols);
  
  init_KMeans_Centers(data,K,nCols,nRows,ss_center);
  return ss_center;

}

//center: current cluster centers
//data: the sample to decide which cluster it belongs to
// nCols: # of features of the sample
// result: the result cluster id
// client_type: Alice or Bob.
void getNearest(Secure_Sample *center, int K, Secure_Sample *data, int nCols, mpz_t domain, int client_type, mpz_t result) 
{
  mpz_t flag;
  mpz_t default_cluster;
  mpz_init(flag);
  mpz_init(default_cluster);

  //initially, set to cluster 0
  mpz_set_ui(result,0);
  SFloat minDistance=s_L2_distance(center+0,data,nCols,client_type);
   
  //check if center of cluster k is nearer
  for(int j=1;j<K;j++)
  {
    if(client_type==ALICE) mpz_set_ui(default_cluster,j);
    SFloat dist=s_L2_distance(center+j,data,nCols,client_type);
    //is dist < minDistance ?
    f_less_than(dist,minDistance,client_type,flag);

    //minDistance = (dist<minDistance) ? dist : minDistance
    minDistance=f_if_then_else(flag,dist,minDistance,client_type);
    //cluster= (dist<min_Distance) ? j : cluster
    s_if_then_else(flag,default_cluster,result,domain,client_type,result);
  }
  mpz_clear(flag);
  mpz_clear(default_cluster); 
}

// compute and set new cluster according to the distance
void updateCluster(Secure_Sample *data, int size, Secure_Sample *center, int nCols, mpz_t *cluster, int K, mpz_t domain, int client_type)
{
  //For all samples, find the nearest cluster centers
  for(int i=0;i<size;i++)
    getNearest(center, K, data+i,nCols,domain,client_type,cluster[i]);

}

//init secure samples from raw samples.
Secure_Sample* setSecureData(Sample* data,int nRows, int nCols)
{
  Secure_Sample* s_Data=new Secure_Sample[nRows];

  for(int i=0;i<nRows;i++)
  {
    s_Data[i].setColLen(nCols);
    for(int j=0;j<nCols;j++) 
      s_Data[i].setX(float2SFloat(data[i].X[j]),j);
  }
  return s_Data;
}

// sData: secured samples
// cluster: cluster[i] is the (added-shared) cluster id of sample i
  //sum: A K*nCols SFloat array that sum[i*nCols+j] (sum(i][j]) will be the sum of cluster i's j-th feature
 void computeSum(SFloat *sum,mpz_t *nC, Secure_Sample *sData,int nCols, mpz_t *cluster, int DataSize,int K,mpz_t domain,int client_type)
{
  mpz_t flag;
  mpz_t default_cluster;
  mpz_t nC_plus_one;
  mpz_init(flag);
  mpz_init(default_cluster);
  mpz_init(nC_plus_one);

  for(int i=0;i<DataSize;i++)//For each sample
  {
    for(int j=0;j<K;j++)//For each cluster
    {
       if(client_type==ALICE) mpz_set_si(default_cluster,j);
       else mpz_set_si(default_cluster,0);

       //is sample i belongs to cluster is j?
       s_equal(cluster[i],default_cluster,domain,client_type,flag);

       //set nC_plus_one = nC[j]+1
       //this way works because nC[j] is added shared
       if(client_type==ALICE) mpz_add_ui(nC_plus_one,nC[j],1);
       else mpz_set(nC_plus_one,nC[j]);
       
       //if flag == 1, then nC[i]=nC[i]+1
       s_if_then_else(flag,nC_plus_one,nC[j],domain,client_type,nC[j]);
       //j-th center's c-th column add s_Data[i]'s c-th column
       for(int c=0;c<nCols;c++)
       {
         //sum[j][c]=(cluster[i] == j) ? sum[j][c]+s_Data[i].X[c] : sum[j][c];
         sum[j*nCols+c]=f_if_then_else(flag,f_plus(sData[i].getX(c),sum[j*nCols+c],client_type),sum[j*nCols+c],client_type);
       }
    }
  }


  mpz_clear(nC_plus_one);
  mpz_clear(default_cluster);
  mpz_clear(flag);
}

//return value = K. Retuen value<0 indicates error
int readSecureKmeans(Secure_Sample **ss_center, const char *mdlfilename,int client_type)
{
  int K=-1;
  FILE *rp=fopen(mdlfilename,"rb");
  if(!rp) { fprintf(stderr,"Error open file %s to read\n", mdlfilename); return K;}
  int nCols=0;
  
  fread(&K,sizeof(int),1,rp);
  fread(&nCols,sizeof(int),1,rp);

  printf("K=%d, nCols=%d\n",K,nCols);

  //if ss_center is not initialized  
  if(!(*ss_center)) {
    printf("allocate memory for secure centers\n");
    *ss_center=new Secure_Sample[K];
    printf("set secure centers\n");
    for(int i=0; i<K;i++) ((*ss_center)[i]).setColLen(nCols);
  }
  else { //make sure all features' memory has been allocated
    printf("set secure centers\n");
    for(int i=0;i<K;i++)
      ((*ss_center)[i]).setColLen(nCols);
  }
  
  //read all centers
  for (int i=0; i<K;i++)
  {
     printf("%d-th center:\n",i);
     for(int j=0; j<nCols; j++) {
       SFloat v;
       v.readRaw(rp);
       printf("\t%d-th column has value: %f\n",j,fpMerge(v,client_type));
       ((*ss_center)[i]).setX(v,j);//set local shared X
     }
  }
  fclose(rp);
  return K;
}

//write the value K and the K centers
int writeKmeanModel(Secure_Sample *ss_center, int K, int nCols, const char *mdlfilename)
{
  FILE *wp=fopen(mdlfilename,"wb");
  if(!wp) { fprintf(stderr,"Error open file %s to write\n", mdlfilename); return -1;}

  //write K
  fwrite(&K,sizeof(int),1,wp);
  fwrite(&nCols,sizeof(int),1,wp);
  //write all centers
  for (int i=0; i<K;i++)
  {
     for(int j=0; j<nCols; j++)
     {
       ss_center[i].getX(j).writeRaw(wp);
     }
  }
  fclose(wp);
  return 1;
}

//return shared cluster ids
int *load_and_predict(const char *mdlfile, const char *testfile, int client_type)
{
  int nCols=getColSize(testfile,"\n\t:;, ");
  int testsize=getRowSize(testfile);
  Secure_Sample *ss_center=NULL;
  int K=readSecureKmeans(&ss_center,mdlfile,client_type);
  
  Secure_Sample *sTestData=NULL;
  //only alice reads the testing data
  if(client_type==ALICE)
  {
    Sample *testData=load_all_samples(testfile,testsize,nCols,false);
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
  printf("start to predict\n");
  fflush(stdout);
  int *result=predict(ss_center,K,sTestData,testsize,nCols,client_type);
  return result;
}

//ss_centers: the K centers
int *predict(Secure_Sample *ss_centers, int K, Secure_Sample *sData, int size, int nCols, int client_type)
{
  mpz_t cluster[size];
  mpz_t domain;
  for(int i=0;i<size;i++)//for each sample, init its cluster id
    mpz_init(cluster[i]);
  mpz_init(domain);
  mpz_set_ui(domain,MAX_CLUSTERS);

  printf("there are %d clusters\n",K);
  printf("there are %d samples, each with %d columns in the testing dataset\n",size, nCols);
 
  updateCluster(sData, size, ss_centers, nCols, cluster, K, domain, client_type);

  int* result=(int*)malloc(sizeof(int)*size);
  for(int i=0;i<size;i++)
    result[i]=mpz_get_ui(cluster[i]);
  
  mpz_clear(domain);
  for(int i=0;i<size;i++)//for each sample, init its cluster id
    mpz_clear(cluster[i]);
  return result;
}

// nIter: the number of iterations to run
// Since the secure Kmeans cannot choose to stop the procedure by checking if the schedule converges (otherwise there will be information leakage), specifying a number of iterations is necessary.
//Bob will generate the init K centers.
Secure_Sample* buildSecureKMeans(const char *train_file, int partyDataSize, int K, int nIter, int client_type)
{
  int localDataSize=getRowSize(train_file);
  int nCols=getColSize(train_file,"\n\t,;: ");
  Sample* trainData=load_all_samples(train_file,localDataSize,nCols,false);

  //init secure party data
  Secure_Sample* s_partyData=new Secure_Sample[partyDataSize];
  for(int i=0;i<partyDataSize;i++)
  {
    s_partyData[i].setColLen(nCols);
    for(int j=0;j<nCols;j++) 
      s_partyData[i].setX(float2SFloat(0),j);
  }

  Secure_Sample* s_Data=setSecureData(trainData,localDataSize,nCols);
  
  Secure_Sample* ss_center=NULL;
  //init centers
  //initial centers will be selected from Bob's data samples
  if(client_type==ALICE)
  {
    ss_center= new Secure_Sample[K];// the K centers
    for(int i=0;i<K;i++)
    {
      ss_center[i].setColLen(nCols);
      for(int j=0;j<nCols;j++)
        ss_center[i].setX(float2SFloat(0),j);
    }
  } 
  else //Bob
    ss_center=setCenter(trainData,K,nCols,localDataSize,client_type);

  //clear unsecure format data
  clearSample(trainData,localDataSize);
  
  int totalDataSize=localDataSize+partyDataSize;
  mpz_t* cluster=(mpz_t*)malloc(sizeof(mpz_t)*totalDataSize);
  mpz_t domain;
  mpz_init(domain);
  mpz_set_ui(domain,MAX_CLUSTERS);

  //init: set all samples will belong to cluster 0 
  for(int i=0;i<totalDataSize;i++) 
      mpz_init(cluster[i]);

  //run nIter iterations
  //In each iteration, first update all sample's cluster, then update the centers

  mpz_t nC[K]; // number of samples in clusters clusters
  for(int i=0;i<K;i++) { mpz_init(nC[i]);} //init
  //A K*nCols SFloat array that sum[i*nCols+j] (sum(i][j]) will be the sum of cluster i's j-th feature
  SFloat* sum=new SFloat[K*nCols];

  for(int iter=0;iter<nIter;iter++)
  {
    //step 1: compute new clusters
    if(client_type==ALICE)
    {
      //update using ALICE's data
      updateCluster(s_Data,localDataSize,ss_center,nCols,cluster,K,domain,client_type);
      //update using Bob's data
      updateCluster(s_partyData,partyDataSize,ss_center,nCols,cluster+localDataSize,K,domain,client_type);
    }
    else // in Bob's case
    {
      //update using ALICE's data
      updateCluster(s_partyData,partyDataSize,ss_center,nCols,cluster,K,domain,client_type);
      //update using Bob's data
      updateCluster(s_Data,localDataSize,ss_center,nCols,cluster+partyDataSize,K,domain,client_type);
    }

    //step 2: update cluster centers
    for(int i=0;i<K*nCols;i++) sum[i]=float2SFloat(0.0); //init
    for(int i=0;i<K;i++) { mpz_set_ui(nC[i],0);} //init # of samples in clusters
    if(client_type==ALICE)
    {
      //set according to Alice's samples
      computeSum(sum,nC,s_Data,nCols,cluster,localDataSize,K,domain,client_type);//use alice's data
      //set according to Bob's samples
      computeSum(sum,nC,s_partyData,nCols,cluster+localDataSize,partyDataSize,K,domain,client_type);//use bob's data
    }
    else //Bob
    {
      //set according to Alice's samples
      computeSum(sum,nC,s_partyData,nCols,cluster,partyDataSize,K,domain,client_type);//use Alice's data
      //set according to Bob's samples
      computeSum(sum,nC,s_Data,nCols,cluster+partyDataSize,localDataSize,K,domain,client_type);//use Bob's data
    }
     
    for(int i=0;i<K;i++)
    {
      double C = mpz_get_d (nC[i]);
      SFloat s_nC=fpShare(C,client_type);
      for(int j=0;j<nCols;j++)
      {
        ss_center[i].setX(f_division(sum[i*nCols+j],s_nC,client_type),j);
      }
    }
  }
  for(int i=0;i<K;i++) { mpz_clear(nC[i]);} //clear
  delete[] sum;

  delete[] s_partyData;
  delete[] s_Data;

  for(int i=0;i< totalDataSize;i++) { 
    mpz_clear(cluster[i]);
  }
  free(cluster);
  mpz_clear(domain);
  return ss_center;
}

