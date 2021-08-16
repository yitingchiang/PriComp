/*
Header file for two-party secure naive bayes model
It is assumed that the label sets in the datasets of the two parties are the same. That is, their datasets have the same number and the same kinds of labels.
*/
#ifndef SECURE_NB
#define SECURE_NB
#include <map>
#include <cstring>

#include "secure_ml.h"
#include "ml_util.h"
#define MAX_LABEL_CHARS 512
#define MAX_FEATURE_CHARS 512

typedef struct
{
    std::map<char*,int,cmpstr> *label_map;
    std::map<char*,int,cmpstr> *feature_map;
    SGaussian* condProb;
    SFloat* prior;
    int nFeature;
    int nClass;
} SecureNB;

void WriteSecureNB(char* file, SecureNB* s_NB);
void ReadSecureNB(char* file, SecureNB* s_NB);
mpz_t* nb_predict(SFloat* data, SecureNB* modeli,int client_type);
SecureNB* TrainSecureNB(const char* file,long party_datasize, int client_type);
double eval_NB(const char* file, SecureNB* mdl, int client_type);

#endif
