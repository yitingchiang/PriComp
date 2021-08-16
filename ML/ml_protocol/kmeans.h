/*
Header file for two-party secure KMeans
It is assumed that the label sets in the datasets of the two parties are the same. That is, their datasets have the same number and the same kinds of labels.
*/
#ifndef SECURE_KMEANS
#define SECURE_KMEANS

#include "ml_util.h"
#include "secure_ml.h"
#define MAX_CLUSTERS 512

//center: current cluster centers
//data: the sample to decide which cluster it belongs to
// nCols: # of features of the sample
// clusterID: clusterID[i]=i. Just for reducing the number of times we need to call mpz_init (to allocate memory)
// result: the result cluster id
// client_type: Alice or Bob.
void getNearest(Secure_Sample *center, int K, Secure_Sample *data, int nCols, mpz_t domain, int client_type, mpz_t result);
Secure_Sample* buildSecureKMeans(const char *train_file,int partyDataSize, int K, int nIter, int client_type);
int *predict(Secure_Sample *ss_centers, int K, Secure_Sample *sData, int size, int nCols, int client_type);
int writeKmeanModel(Secure_Sample *ss_center, int K, int nCols, const char *mdlname);
//set secure samples from raw samples
Secure_Sample* setSecureData(Sample* data,int nRows, int nCols);
//load model files and testing files and return the predicted cluster ids
int *load_and_predict(const char *mdlfile, const char *testfile, int client_type);
//return K
int readSecureKmeans(Secure_Sample *ss_center, const char *mdlfilename,int client_type);

#endif
