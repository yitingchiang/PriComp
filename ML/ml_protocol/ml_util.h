#ifndef MLUTIL
#define MLUTIL
#include <map>
#include <cstring>
#include <cstdbool>
struct cmpstr
{
  bool operator() (char* s1, char* s2) const
  {
    return strcmp(s1,s2)<0;
  }
};

typedef struct
{
  double* X;
  double sq_sum;
  char* label;
} Sample;

int getRowSize(const char* file);
int getColSize(const char* file,const char* sep);
long getLabels(const char* datafile,const char* sep,std::map<char*,int,cmpstr>* &ele_map);
void create_map(const char* datafile,const char* sep,std::map<char*,int,cmpstr> *ele_map, int nCol);
int collect_by_label(const char* datafile,const char* sep,std::map<char*,int,cmpstr> label_map, std::map<char*,int,cmpstr>* feature_map, int nf);
int* GetAndSet(std::map<char*,int,cmpstr> *map, int set);
double* getLine_double(FILE* rp, int nf, char* label);
Sample* load_all_samples(const char* file, int nRow, int nF, bool isLabeled);
std::map<char*,int,cmpstr>* setLabelIDMap(const char* file, const char* sep,int setID);
void clearSample(Sample* s,int size);
#endif
