#include "secure_ml.h"

Secure_Sample::Secure_Sample(int nCols,bool maintain_sq) {
  this->X.resize(nCols);
  this->maintain_square=maintain_sq;
  if(this->maintain_square)
  {
    this->square_sumX=float2SFloat(0.0);
    this->square_X.resize(nCols);
  }
}

void Secure_Sample::setMaintainSquare(bool b)
{
    maintain_square=b;
}

bool Secure_Sample::isMainTainSquare()
{
    return maintain_square;
}

int Secure_Sample::getColLen()
{
  return X.size();
}

Secure_Sample::~Secure_Sample() {
  this->X.clear();
  if(this->maintain_square) 
    this->square_X.clear();
}

void Secure_Sample::setColLen(int n) {
  this->X.resize(n);
  this->square_X.resize(n);
}

void Secure_Sample::setY(SFloat v) {
  this->Y=v;
}

void Secure_Sample::setX(SFloat v, int idx) {
  this->X[idx]=v;//update X[idx]
  //compute the new square sum 
  if(maintain_square) {
    double sq_sumX=sFloat2float(this->square_sumX);
    double sq_X=sFloat2float(this->square_X[idx]);
    double sq_v=sFloat2float(v);
    sq_v*=sq_v;
    this->square_X[idx]=float2SFloat(sq_v);
    this->square_sumX=float2SFloat(sq_sumX-sq_X+sq_v);
  }
}

void Secure_Sample::setX(SFloat v, int idx, int client_type) {
  this->X[idx]=v;//update X[idx]
  if(maintain_square) {
    //compute the new square sum 
    this->square_sumX=f_minus(this->square_sumX,this->square_X[idx],client_type);//square_sumX-=old X[idx]^2
    this->square_X[idx]=f_product(v,v,client_type);//update X[idx]^2

    this->square_sumX=f_plus(this->square_sumX,this->square_X[idx],client_type);//new square_sumX=square_sumX + new X[idx]^2
  }
}

//generate secure gaussian using the label map (lmap), the feature map (fmap).
//lcount[i] is the number of label i appeared in the dataset
SGaussian* createSGaussian(std::map<char*,int,cmpstr> lmap, std::map<char*,int,cmpstr> *fmap, int* lcount, int nf,int party_datasize,int local_datasize,int client_type)
{
    int nC=lmap.size();
    SGaussian* sg=new SGaussian[nC*nf];
    SFloat zero=float2SFloat(0.0);
    std::map<char*,int,cmpstr>::iterator lit=lmap.begin();
    SFloat SCSize;
    for(int i=0;lit!=lmap.end();lit++,i++)
    {
        int fidx=lit->second;
        SCSize=float2SFloat((double)lcount[i]);

        if(client_type==ALICE) 
        {
            SCSize=f_plus(SCSize,zero,client_type);//# of samples of thi label
        }
        else
        {
            SCSize=f_plus(zero,SCSize,client_type);
        }

        for(int j=0;j<nf;j++)
        {
            std::map<char*,int,cmpstr> fm=fmap[fidx*nf+j];
            int len=fm.size();
            double sum=0;
            double sq_sum=0;

            for(std::map<char*,int,cmpstr>::iterator fit=fm.begin();fit!=fm.end();fit++)
            {

              double v=(double)strtod(fit->first,NULL);
              sum+=v*fit->second;
              sq_sum+=(v*v*fit->second);//local square sum

            }

            SFloat EX2;
            SFloat E2X=float2SFloat(sum);
            if(client_type==ALICE)
            {
               
              EX2=f_plus(float2SFloat(sq_sum),zero,client_type);//total sum(X^2)
              E2X=f_plus(E2X,zero,client_type);
            }
            else
            {
              EX2=f_plus(zero,float2SFloat(sq_sum),client_type);
              E2X=f_plus(zero,E2X,client_type);
            }
            SGaussian *g=sg+(i*nf+j);//point to current secure gaussian
            
            
            E2X=f_division(E2X,SCSize,client_type);//E(X)
            
  
            g->setMean(E2X);
            EX2=f_division(EX2,SCSize,client_type);//secure sum(X^2)/N
            E2X=f_product(E2X,E2X,client_type);//E^2(X)
            g->setVar(f_minus(EX2,E2X,client_type));//E(X^2)-E^2(X)

        }
    }
    return sg;
} 

//get one row and convert it into SFloat datatype.
//it is assumed that these are numeric features and the label is in the last column.
//Note: nf does not count the label, and the memory of the returned array is allocated automatically by using new
SFloat* getLine(FILE* rp, int nf, char* label)
{
    const char* sep="\n\r\t ,;";
    char *line = NULL;
    size_t len = 0;
    SFloat* data=new SFloat[nf];
    if(!label)  { fprintf(stderr,"getLine() cannot accept a null string parameter\n"); return NULL;}

    if(getline(&line,&len,rp)==-1) { delete[] data; return NULL;//EOF 
    }
    if(strlen(line)<nf) { delete[] data; free(line); return NULL;}//skip null line

    char* tmps=strtok(line,sep);
    for(int i=0;i<nf;i++)
    {
        data[i]=float2SFloat(strtod(tmps,NULL));
        tmps=strtok(NULL,sep);
    }
    
    strcpy(label,tmps);
    free(line);
    return data; 
}

// secure sigmoid function
SFloat s_sigmoid(SFloat v, int client_type)
{
    SFloat ONE = float2SFloat(1.0);
    SFloat ZERO = float2SFloat(0.0);
    SFloat exp_x=f_exp(v,client_type); // exp(x)
    SFloat tmp=(client_type==ALICE) ? f_plus(ONE,exp_x,client_type): f_plus(ZERO,exp_x,client_type); //1+exp(x)
    SFloat result=f_division(exp_x,tmp,client_type);
    return result;
}

//f(x)=sign(x). That is, return 1 if x >0. oterwise 0
SFloat s_sign(SFloat v, int client_type)
{
    return f_is_positive(v,client_type);
}

SFloat s_relu(SFloat v, int client_type)
{
    return f_relu(v,client_type);
}

//s1, s2: two samples in unscured form
//nF: the number of features
//return L2^2
SFloat s_L2_distance(Sample *s1, Sample *s2, int nF, int client_type)
{
    double sum=0;
    //compute 2*sum xi*yi
    for(int i=0;i<nF;i++)
    {
      double v=2*s1->X[i]*s2->X[i];
      sum+=v;  
    }
    
    //result=\sum (xi^2+xi*yi+yi^2)
    double r=s1->sq_sum+sum+s2->sq_sum;
    SFloat result=fpShare(sum,client_type);
   
    return result;
}

//s1, s2: two samples in unscured form
//nF: the number of features
//return L2^2
SFloat s_L2_distance(Secure_Sample *s1, Secure_Sample *s2, int nF, int client_type)
{
    SFloat TWO=(client_type==ALICE) ? float2SFloat(2.0) : float2SFloat(0);
    SFloat sum=float2SFloat(0.0);
    for(int i=0;i<nF;i++)
    {
      //xi*yi
      SFloat v=f_product(s1->getX(i),s2->getX(i),client_type);
      //2*xi*yi
      v=f_product(TWO,v,client_type);
      //sum=sum-2*xi*yi
      sum=f_minus(sum,v,client_type);
      //xi*xi
      v=f_product(s1->getX(i),s1->getX(i),client_type);
      //sum=sum+xi*xi
      sum=f_plus(v,sum,client_type);

      //yi*yi
      v=f_product(s2->getX(i),s2->getX(i),client_type);
      //sum=sum+yi*yi
      sum=f_plus(v,sum,client_type);
    }
    
    return sum;
}
