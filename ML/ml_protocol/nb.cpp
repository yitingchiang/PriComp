#include "nb.h"
#include <cstring>

int local_discrete_nb(const char* datafile, const char* mdl)
{
   int i=0, j=0;
   const char* sep="\n\r\t ,;";
   int nCol=getColSize(datafile,sep);
   int cnt=0;

   if(nCol<0) return -1;
   std::map<char*,int,cmpstr> *ele_map=new std::map<char*,int,cmpstr>[nCol];
  
   create_map(datafile,sep,ele_map,nCol);

   //clear maps
   for(i=0;i<nCol;i++)
       for(std::map<char*,int,cmpstr>::iterator it=ele_map[i].begin();it!=ele_map[i].end();it++)
       {
           delete[] it->first;
           ele_map[i].erase(it);
       }

   delete[] ele_map;
   return 0;
}

SFloat* computePrior(int* lcount, int nC, long N, int client_type)
{
    SFloat* p=new SFloat[nC];
    for(int i=0;i<nC;i++)
    {
        double v=(double)lcount[i]/N;
        p[i]=f_average_2(1,&v,1,client_type);
    }
    return p;
}

void nb_predict(SFloat* data, SecureNB* mdl, int client_type, mpz_t result)
{
    int nC=mdl->nClass;
    int nf=mdl->nFeature;
    std::map<char*,int,cmpstr>* lmap=mdl->label_map;
    std::map<char*,int,cmpstr> *fmap=mdl->feature_map;
    SGaussian* condProb=mdl->condProb;
    SFloat* prior=mdl->prior;
    SFloat* prob=new SFloat[nC];

    mpz_t candidate;
    mpz_t flag;
    mpz_t domain;
    mpz_init(flag);
    mpz_init(candidate);
    mpz_init(domain);
    mpz_set_ui(domain,nC);
    for(int i=0;i<nC;i++)
    {
        prob[i]=float2SFloat(0.0);
        for(int j=0; j<nf; j++)
        {
            //get P(x_j|C_i)
            SFloat v=mdl->condProb[i*nf+j].getLnProb(data[j],client_type);

            prob[i]=f_plus(prob[i],v,client_type);
        }
        //sum_j ln(P(X_j|C_i))+ ln(P(C_i))
        SFloat tmp=f_ln(prior[i],client_type);
        prob[i]=f_plus(prob[i],tmp,client_type);
    }

    SFloat max;
    max.set(prob[0].s,prob[0].e,prob[0].m);

    //initialize the final result to class 0
    mpz_set_ui(result,0);

    for(int i=1;i<nC;i++)
    {
       //set candidate winner
       if(client_type==ALICE) mpz_set_ui(candidate,i);
       else mpz_set_ui(candidate,0);


       SFloat tmp=f_minus(max,prob[i],client_type);//max-prior[i]
       f_is_negative(tmp,client_type,flag);//is max<prior[i]?
       max=f_if_then_else(flag,prob[i],max,client_type);
       s_if_then_else(flag,candidate,result,domain,client_type,result);

    }
    mpz_clear(flag);    
    mpz_clear(candidate);    
    mpz_clear(domain);    
    

    return ;
}

SecureNB* TrainSecureNB(const char* file,long party_datasize, int client_type)
{
    std::map<char*,int,cmpstr> *label_map;
    std::map<char*,int,cmpstr> *feature_map;
    const char* sep="\n\r\t ,;";
    int nf=getColSize(file,sep)-1;//# of features
    int nC=0;
    long local_datasize=0;

    label_map=new std::map<char*,int,cmpstr>;
    if((local_datasize=getLabels(file,sep,label_map))<0) return NULL;
    feature_map=new std::map<char*,int,cmpstr>[nf*label_map->size()];

    nC=label_map->size();
    
    int* label_count=GetAndSet(label_map,true);

    SecureNB* mdl=(SecureNB*)malloc(sizeof(SecureNB));
    mdl->nFeature=nf;
    mdl->nClass=nC;
    mdl->label_map=label_map;
    mdl->feature_map=feature_map;

    if(collect_by_label(file,sep,*label_map,feature_map,nf)<0) return NULL;

    mdl->prior=computePrior(label_count,nC,local_datasize+party_datasize,client_type);


    std::map<char*,int,cmpstr>::iterator it=label_map->begin();

    //generate nC secure gaussian distributions
    mdl->condProb=createSGaussian(*label_map,feature_map,label_count,nf,party_datasize,local_datasize,client_type);

    return mdl;
}


void ReadSecureNB(char* file, SecureNB* s_NB)
{
    std::map<char*,int,cmpstr> *label_map=NULL;
    std::map<char*,int,cmpstr> *feature_map=NULL;
    SGaussian* condProb=NULL;
    SFloat* prior=NULL;
    int i=0, j=0;
    FILE* rp=fopen(file,"rb");
    if(!rp) { fprintf(stderr,"ERROR! ccan not read local model file %s!\n",file); return ;}

    fread(&(s_NB->nFeature),sizeof(int),1,rp);
    fread(&(s_NB->nClass),sizeof(int),1,rp);

    int nF=s_NB->nFeature;
    int nC=s_NB->nClass;
    char label[MAX_LABEL_CHARS];
    char feature[MAX_FEATURE_CHARS];

    //allocate memory
    s_NB->feature_map=new std::map<char*,int,cmpstr>[nC*nF];
    s_NB->label_map=new std::map<char*,int,cmpstr>;
    s_NB->condProb=new SGaussian[nC*nF];
    s_NB->prior=new SFloat[nC];

    feature_map=s_NB->feature_map;
    label_map=s_NB->label_map;
    condProb=s_NB->condProb;
    prior=s_NB->prior;

    int value=0;
    char* tmps=NULL;
    for(i=0;i<nC;i++)
    {
        fread(label,sizeof(char)*MAX_LABEL_CHARS,1,rp);
        fread(&value,sizeof(int),1,rp);
        tmps=(char*)malloc(sizeof(char)*strlen(label));
        strcpy(tmps,label);
        label_map->insert(std::pair<char*,int>(tmps,value));
        prior[i].readRaw(rp);

        for(j=0;j<nF;j++)
        {
            int idx=i*nF+j;
            int len;
            fread(&len,sizeof(int),1,rp);

            for(int k=0;k<len;k++)
            {
                 fread(feature,sizeof(char)*MAX_FEATURE_CHARS,1,rp);
                 fread(&value,sizeof(int),1,rp);
                 tmps=(char*)malloc(sizeof(char)*strlen(feature));
                 strcpy(tmps,feature);
                 feature_map[idx].insert(std::pair<char*,int>(tmps,value));
            }
            condProb[idx].readRaw(rp);
        }
    }
    fclose(rp);
}

void WriteSecureNB(char* file, SecureNB* s_NB)
{
    std::map<char*,int,cmpstr> *label_map=s_NB->label_map;
    std::map<char*,int,cmpstr> *feature_map=s_NB->feature_map;
    SGaussian* condProb=s_NB->condProb;
    SFloat* prior=s_NB->prior;
    int nF=s_NB->nFeature;
    int nC=s_NB->nClass;
    int i=0, j=0;
   
    FILE* wp=fopen(file,"wb");
    if(!wp) { fprintf(stderr,"ERROR! ccan not create local model file %s!\n",file); return ;}

    fwrite(&nF,sizeof(int),1,wp);
    fwrite(&nC,sizeof(int),1,wp);
    std::map<char*,int,cmpstr>::iterator lit=label_map->begin();
    for(i=0;i<nC;i++,lit++)
    {
        fwrite(lit->first,sizeof(char)*MAX_LABEL_CHARS,1,wp);
        fwrite(&(lit->second),sizeof(int),1,wp);
        prior[i].writeRaw(wp);
        for(j=0;j<nF;j++)
        {
            int idx=i*nF+j;
            int nValue=feature_map[idx].size();
            fwrite(&nValue,sizeof(int),1,wp);
            std::map<char*,int,cmpstr>::iterator fit=feature_map[idx].begin();
            for(int k=0;k<nValue;k++,fit++)
            {
                 fwrite(fit->first,sizeof(char)*MAX_FEATURE_CHARS,1,wp);
                 fwrite(&(fit->second),sizeof(int),1,wp);
            }
            
            condProb[idx].writeRaw(wp);
        }
    }
    fclose(wp);

}

double eval_NB(const char* file, SecureNB* mdl, int client_type)
{
    const char* sep="\n\r\t ,;";
    int nf=getColSize(file,sep)-1;//# of features

    fprintf(stderr,"chk\n")   ;
    fflush(stderr);
    SFloat* data=NULL;
    SFloat* zeros=new SFloat[nf];
    std::map<char*,int,cmpstr> *label_map=mdl->label_map;
    std::map<char*,int,cmpstr> *feature_map=mdl->feature_map;

    char label[256];
    int nC=mdl->nClass;
    FILE* rp=fopen(file,"r");
    if(!rp) { fprintf(stderr,"cannot open file %s\n",file); return -1.0;}

    mpz_t predict;
    mpz_t domain;
    mpz_init(predict);
    mpz_init(domain);
    mpz_set_ui(domain,nC);
    double acc;
    int t_size=0;
    int hit=0; 
    for(int i=0;i<nf;i++) zeros[i]=float2SFloat(0.0);

    for(t_size=0;(data=getLine(rp,nf,label))!=NULL;t_size++)
    {
        if(!data) { break;}
      
        //fflush(stdout);
        if(client_type==ALICE)
            nb_predict(data,mdl,client_type,predict);
        else
            nb_predict(zeros,mdl,client_type,predict);
        std::map<char*,int,cmpstr>::iterator it=label_map->find(label);

        if(it==label_map->end())//label in testing data is not exists in the label set, which means something wrong!
        {
            fprintf(stderr,"ERROR in eval: get an unknown label %s from the testing data!\n",label);
            t_size--;//decrease the counter
        }
        else
        {
            //parties get the true label
            int predict_label=sMerge(predict,domain,client_type);
            if(it->second==predict_label) 
            {
                hit++;
            }
        }
    }
    mpz_clear(predict);
    mpz_clear(domain);
    delete[] zeros;
    return (double)hit/t_size;
}

