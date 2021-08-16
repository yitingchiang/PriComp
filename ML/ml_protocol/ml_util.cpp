#include <cstdio>
#include <cstring>
#include <cstdbool>
#include "ml_util.h"

//count # of rows in a dataset file
//return -1 if there is error
int getRowSize(const char* file)
{
   char *line = NULL;
   size_t len = 0;
   FILE* rp=fopen(file,"r");
   int nRow=0;
   if(!rp) { fprintf(stderr,"Cannot open file %s\n",file); return -1;}
   while(getline(&line,&len,rp)!=-1)
       nRow++;

   fclose(rp);
   return nRow;
}

void clearSample(Sample* s,int size) 
{
  for(int i=0;i<size;i++)
  {
    if(s[i].label) free(s[i].label);
    if(s[i].X) free(s[i].X);
  }
  free(s);
}

//count the number of columns according to the delimiters 
//return -1 if there is error
int getColSize(const char* file,const char* sep)
{
    char *line = NULL;
    size_t len = 0;
    FILE* rp=fopen(file,"r");
    int nCol=0;
    char* tmps=NULL;
    if (!rp) 
    {
        fprintf(stderr,"Cannot open file %s!\n",file);
        return -1;
    }

    if(getline(&line,&len,rp)==-1)
    {
        fprintf(stderr,"Failed to read any data from file %s!\n",file);
        fclose(rp);
        return -1;
    }
    fclose(rp);
    tmps=strtok(line,sep);
    for(nCol=0;tmps!=NULL;)
    {
       nCol++;
       tmps=strtok(NULL,sep);
    }
    free(line);
    return nCol;
}

//create a map<char*,int> from a given column in a dataset file 
void create_map(const char* datafile,const char* sep,std::map<char*,int,cmpstr> *ele_map, int nCol)
{
   FILE* rp=fopen(datafile,"r");
   char *line = NULL;
   size_t len = 0;
   
   for(int j=0;getline(&line,&len,rp)!=-1;j++)
   {
       if(strlen(line)<nCol) continue;
       char* tmps=strtok(line,sep);
       for(int i=0;i<nCol;i++)
       {
           std::map<char*,int,cmpstr>::iterator it=ele_map[i].find(tmps);
           if(it==ele_map[i].end())//the key does not exist
           {
               char* data=new char[strlen(tmps)];
               strcpy(data,tmps);
               ele_map[i].insert(std::pair<char*,int>(data,1));
           }
           else it->second+=1;//the key exists
           tmps=strtok(NULL,sep);
       }
   }
   fclose(rp);
}


//get and store the labels in a map<str,int> container where the second component is the number of that label in the samples.
//assume the label is in the last column
//return -1 for error. otherwise, return # of samples
long getLabels(const char* datafile,const char* sep,std::map<char*,int,cmpstr>* &ele_map)
{
   char *line = NULL;
   size_t len = 0;
   int nCol=getColSize(datafile,sep);
   if(nCol<0) { fprintf(stderr,"Failed to get the number of columns from file %s\n",datafile); return -1;}
   FILE* rp=fopen(datafile,"r");
   long datasize=0;
   for(datasize=0;getline(&line,&len,rp)!=-1;datasize++)
   {
       if(strlen(line)<nCol) continue;
       char* tmps=strtok(line,sep);
       //skip non-label columns
       for(int i=0;i+1<nCol;i++) tmps=strtok(NULL,sep);
       
       std::map<char*,int,cmpstr>::iterator it=ele_map->find(tmps);
       if(it==ele_map->end())//the key does not exist
       {
           char* data=new char[strlen(tmps)];
           strcpy(data,tmps);
           ele_map->insert(std::pair<char*,int>(data,1));
       }
       else it->second+=1;//the key exists
       tmps=strtok(NULL,sep);
   }
   fclose(rp);
   return datasize;
}


//label_map: <char*,int> that the first component (char*) is the name of the label, and the second component (int) is the index of the label in feature_map;
//feature_map: a length-(n*m), where n is the number of labels and m is the number of features
int collect_by_label(const char* datafile,const char* sep,std::map<char*,int,cmpstr> label_map, std::map<char*,int,cmpstr>* feature_map, int nf)
{
   char *line = NULL;
   size_t len = 0;
   int nC=label_map.size();
   FILE* rp=fopen(datafile,"r");
   
   for(int j=0;getline(&line,&len,rp)!=-1;j++)
   {
       if(strlen(line)<nf) continue;
       char* data=(char*)malloc(sizeof(char)*strlen(line));
       strcpy(data,line);
       char* tmps=strtok(line,sep);
       
       //skip non-label columns
       for(int i=0;i<nf;i++) tmps=strtok(NULL,sep);
       
       std::map<char*,int,cmpstr>::iterator it=label_map.find(tmps);

       if(it==label_map.end())//the key does not exist
       {
           fprintf(stderr,"ERROR. label value %s not found!\n",tmps);
           fflush(stderr);
           fclose(rp);
           free(data);
           return -1;
       }
       else
       {
           int lab_idx=it->second;
           tmps=strtok(data,sep);
           for(int i=0;i<nf;i++) 
           {
               int feature_idx=lab_idx*nf+i;
               std::map<char*,int,cmpstr>::iterator it=feature_map[feature_idx].find(tmps);
               if(it==feature_map[feature_idx].end())//a new value
               {
                   char* data=new char[strlen(tmps)];
                   strcpy(data,tmps);
                   feature_map[feature_idx].insert(std::pair<char*,int>(data,1));
               }
               else it->second+=1;
               tmps=strtok(NULL,sep);
           }
       }
   }
   fclose(rp);
   return 1;
}


//if set != 0, the second element of the map will be set to a serial number
//as a result, the first element will be mapped to a int (serial number)
//Anyway, an int array count will be the second element of the map
int* GetAndSet(std::map<char*,int,cmpstr> *map, int set)
{
    int len=map->size();
    int* count=new int[len];

    for(int i=0;i<len;i++) count[i]=0;
    
    std::map<char*,int,cmpstr>::iterator it=map->begin();
    for(int i=0;i<len;i++,it++)
    {
        count[i]=it->second;//get
        if(set) it->second=i;//set
    }

    return count;
}

//get one row and convert it into double features.
//it is assumed that these are numeric features and the label is in the last column.
//Note: nf does not count the label, and the memory of the returned array is allocated automatically by using new 
double* getLine_double(FILE* rp, int nf, char* label)
{
    const char* sep="\n\r\t ,;";
    char *line = NULL;
    size_t len = 0;
    double* data=new double[nf];

    if(getline(&line,&len,rp)==-1) { delete[] data; return NULL;//EOF 
    }
    if(strlen(line)<nf) { fprintf(stderr,"get a null line!\n"); delete[] data; free(line); return NULL;}//skip null line

    char* tmps=strtok(line,sep);
    for(int i=0;i<nf;i++)
    {
        data[i]=strtod(tmps,NULL);
        tmps=strtok(NULL,sep);
    }
   
    if(tmps&&label) strcpy(label,tmps);
    free(line);
    return data; 
}

std::map<char*,int,cmpstr>* setLabelIDMap(const char* file, const char* sep,int setID)
{
   std::map<char*,int,cmpstr>* label_map=new std::map<char*,int,cmpstr>;

   if(getLabels(file,sep,label_map)<0) return NULL;

   int* count=GetAndSet(label_map,setID);//set the second element of the map to a serial number
   delete[] count;//we do not use the number of samples of each label

   return label_map;
}

// read all samples 
Sample* load_all_samples(const char* file, int nRow, int nF, bool isLabeled)
{
    int i=0;
    Sample* sample=new Sample[nRow];
    char* label;
    FILE* rp=fopen(file,"r");
     
    for(i=0;i<nRow;i++)
    {
        if(isLabeled) sample[i].label=new char[256];
        else sample[i].label=NULL;
        if((sample[i].X=getLine_double(rp,nF,sample[i].label))==NULL)
        {
          fprintf(stderr,"Error in loading data samples from file %s.\n",file);
          continue;
        }
    }
    fclose(rp);
    return sample;  
}

