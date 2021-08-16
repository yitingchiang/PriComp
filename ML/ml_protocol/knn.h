/*
Header file for two-party secure KNN model
It is assumed that the label sets in the datasets of the two parties are the same. That is, their datasets have the same number and the same kinds of labels.
*/
#ifndef SECURE_KNN
#define SECURE_KNN
#include <map>
#include <cstring>

#include "ml_util.h"
#include "secure_ml.h"
#define MAX_LABEL_CHARS 512
#define MAX_FEATURE_CHARS 512

typedef struct
{
   Sample *s;
   int label;
   double distance;
} Candidate;

Sample* init_KNN(const char* train_file, int &nRow, int &nFeatures, int K);
int* secure_KNN(const char* train_file, const char* test_file, int K, int client_type);
#endif
