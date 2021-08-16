#include "knn.h"
#include <cstring>
#include <cmath>
#include <limits>

//Euclidean distance
double L2_distance(double* p1, double* p2, int nF)
{
   int i=0;
   double sum=0;
   for(i=0;i<nF;i++)
   {
      double v=p1[i]-p2[i];
      sum+=(v*v);
   }

   return sqrt(sum);
}

//manhattan distance
double L1_distance(double* p1, double* p2, int nF)
{
   int i=0;
   double sum=0;
   for(i=0;i<nF;i++)
   {
       double v=abs(p1[i]-p2[i]);
       sum+=v;
   }
   return sum;
}


//find in a descendently sorted array
double binary_search(double* a, int l, int r, double v)
{
    while(l<r)
    {
        int mid=(l+r)/2;
        if(a[mid]<v) r=mid-1;
        else if(a[mid] > v) l=mid+1;
        else return mid;
    }
}

void swap_candidate(Candidate* c1, Candidate* c2)
{
    Sample* tmps=c1->s;
    double tmpd=c1->distance;
    int tmpi=c1->label;
    
    c1->s=c2->s;
    c1->distance=c2->distance;
    c1->label=c2->label;

    c2->s=tmps;
    c2->distance=tmpd;
    c2->label=tmpi;
}

int local_vote(Candidate* c,int K,std::map<char*,int,cmpstr> *label_map)
{
    int lsize=label_map->size();
    int* votes=new int[lsize];

    for(int i=0;i<lsize;i++) votes[i]=0;

    for(int i=0;i<K;i++)
    {
        std::map<char*,int,cmpstr>::iterator it=label_map->find(c[i].s->label);
        if(it==label_map->end())
        {
            fprintf(stderr,"In vote(): unknown label: \"%s\"",c[i].s->label);
            continue;
        }

        votes[it->second]++;
        c[i].label=it->second;
    }

    int MAX=0;
    int MAX_IDX=-1;
    for(int i=0;i<lsize;i++)
    {
        if(votes[i]>MAX) { MAX=votes[i]; MAX_IDX=i;}
    }
    delete[] votes;

    return MAX_IDX;
}

int* SecureVote(Candidate* c,int K,int nC, int client_type)
{
    mpz_t* alice_svotes=new mpz_t[K];
    mpz_t* bob_svotes=new mpz_t[K];
    mpz_t* merged_svotes=new mpz_t[K];
    mpz_t is_less;
    mpz_t domain;
    SFloat* alice_dist=new SFloat[K];
    SFloat* bob_dist=new SFloat[K];
    SFloat* merged_dist=new SFloat[K];
    SFloat max;
    max.setMax();

    mpz_init(is_less);
    mpz_init(domain);
    mpz_set_ui(domain,nC);

    int* result=new int[K];

    //init the distance 
    for(int i=0;i<K;i++)
    {
        mpz_init(alice_svotes[i]);
        mpz_init(bob_svotes[i]);
        mpz_init(merged_svotes[i]);
        if(client_type==ALICE)
        {
          mpz_set_ui(alice_svotes[i],c[i].label);
          mpz_set_ui(bob_svotes[i],0);
          alice_dist[i]=float2SFloat(c[i].distance);
          bob_dist[i]=float2SFloat(0.0);
        }
        else
        {
          mpz_set_ui(bob_svotes[i],c[i].label);
          mpz_set_ui(alice_svotes[i],0);
          bob_dist[i]=float2SFloat(c[i].distance);
          alice_dist[i]=float2SFloat(0.0);
        }

    }

    //from the two k-candidates, select the top-K smallest distance.
    for(int i=0;i<K;i++)
    {
        //comparing with alice's vote bob's vote 
        f_less_than(alice_dist[0],bob_dist[0],client_type,is_less);//is alice_dist[0]< bob_dist[0]?
        //save the smaller one
        merged_dist[i]=f_if_then_else(is_less,alice_dist[0],bob_dist[0],client_type);
        //get the i-th vote
        s_if_then_else(is_less,alice_svotes[0],bob_svotes[0],domain,client_type,merged_svotes[i]);

        //shift elements:
        //if alice_data[0]<bob_data[0], then get alice_data[0] and move alice_data[0]<-alice_data[1], alice_data[1]<-alice_data[2],...
        //the same procedure is also applied to bob
        //only move from idx 1 (1 to 0), 2 (2 to 1),..., to idx K-i (K-i to K-1-i)
        for(int j=0;j<K-1-i;j++)
        {
            alice_dist[j]=f_if_then_else(is_less,alice_dist[j+1],alice_dist[j],client_type);
            s_if_then_else(is_less,alice_svotes[j+1],alice_svotes[j],domain,client_type,alice_svotes[j]);
            bob_dist[j]=f_if_then_else(is_less,bob_dist[j],bob_dist[j+1],client_type);
            s_if_then_else(is_less,bob_svotes[j],bob_svotes[j+1],domain,client_type,bob_svotes[j]);
        }
  
    }  
   
    for(int i=0;i<K;i++)
    {
        result[i]=sMerge(merged_svotes[i],domain,client_type);
    }

    mpz_clear(is_less);
    mpz_clear(domain);
    for(int i=0;i<K;i++)
    {
        mpz_clear(merged_svotes[i]);
        mpz_clear(alice_svotes[i]);
        mpz_clear(bob_svotes[i]);
    }
    return result;
}

//load all samples
Sample* init_KNN(const char* train_file, int &nRow, int &nFeatures, int K)
{
   const char* sep="\n\r\t ,;";
   nFeatures=getColSize(train_file, sep)-1;//do not count the label column
   if(nFeatures<0) return NULL;
   nRow=getRowSize(train_file);
   if(nRow<0) return NULL;

   Sample* all_samples=load_all_samples(train_file,nRow,nFeatures, true);
   if(!all_samples) return NULL;

   return all_samples;
}

int vote2predict(int* vote, int K, int len)
{
    int max=0;
    int max_idx=-1;
    int* n=new int[len];

    for(int i=0;i<len;i++) n[i]=0;

    for(int i=0;i<K;i++) 
       n[vote[i]]++; 

    for(int i=0; i<len;i++)
    {
        if(n[i]>max) {max=n[i]; max_idx=i;}
    }
    delete[] n;
    return max_idx;
}

//return prediction
//the secure KNN method first run KNN locally, and then the two parties start a secure vote procedure to decide the output.
int* secure_KNN(const char* train_file, const char* test_file, int K, int client_type)
{
   const char* sep="\n\r\t ,;";
   int nF;
   int nRow;

   Sample* all_samples=init_KNN(train_file,nRow,nF,K);
   if(!all_samples) return NULL;
   
   std::map<char*,int,cmpstr>* label_map=setLabelIDMap(train_file,sep,true);
   if(!label_map) return NULL;
   
   FILE* rp=NULL;
   if((rp=fopen(test_file,"r"))==NULL)
   {
       fprintf(stderr,"Error! Cannot open test file %s.\n",test_file);
       return NULL;
   }

   Candidate* c=new Candidate[K];

   char label[256];
   int nTest=0, hit=0;
   double* test_sample;
  
   nTest=getRowSize(test_file);
   int* predict=new int[nTest];
   for(int idx=0;(test_sample=getLine_double(rp,nF,label))!=NULL;idx++)//get a new sample from the testing data to predict
   {
       //init k candidates
       for(int i=0;i<K;i++) 
       {
           c[i].distance=std::numeric_limits<double>::max(); 
           c[i].s=NULL; 
       }

       for(int i=0;i<nRow;i++)
       {
         if(all_samples[i].X==NULL) //corresponding to a null line
             continue;
         double d=L2_distance(test_sample,all_samples[i].X,nF);
         Candidate tmpc;
         std::map<char*,int,cmpstr>::iterator it=label_map->find(all_samples[i].label);
         tmpc.s=all_samples+i;
         tmpc.distance=d;
         tmpc.label=it->second;
         for(int j=0;j<K;j++)
         {
           if(d<c[j].distance)
           {
              if(c[j].s==NULL)//not yet found K candidates
              {
                  c[j].s=tmpc.s;
                  c[j].distance=d;
                  c[j].label=tmpc.label;
                  break;
              }
              else
              {
                  swap_candidate(c+j,&tmpc);
                  d=tmpc.distance;
              }
           }
         }
       }

       std::map<char*,int,cmpstr>::iterator it=label_map->find(label);
       if(it==label_map->end())//do not find a known label
       {
           fprintf(stderr,"unknown label %s!",label);
           continue;
       }
       int lid=it->second;

       //local KNN is ready. run corss-party KNN\n;
       int* svote=SecureVote(c,K,label_map->size(),client_type);
       
       predict[idx]=vote2predict(svote,K,label_map->size());
       if(lid==predict[idx]) hit++;
       delete[] svote;
       delete[] test_sample;
   }
   fclose(rp);

   printf("Accuracy = %f (%d/%d)\n",(double)hit/nTest,hit,nTest);

   for(int i=0;i<nRow;i++)
   {
       delete[] all_samples[i].X;
       delete[] all_samples[i].label;
   }
   delete[] all_samples;

   delete label_map;
 
   return predict;
}
