#ifndef SECURE_ML
#define SECURE_ML
#include "ml_util.h"
#include "two_parties/two_parties.h"
#include "smc.h"
#include "debug_util.h"
#include "util.h"
#include <cstdbool>
#include <vector>
//Gaussian distribution with secure-type parameters
class SGaussian 
{
    private:
      SFloat mean;
      SFloat var;    
      double PI=3.14159265;
    public:
      void setMean(SFloat s)
      {
          mean=s;
      }
      void setVar(SFloat s)
      {
          var=s;
      }
      SFloat getMean()
      {
          return mean;
      }
      SFloat getVar()
      {
          return var;
      }
      void write(FILE* wp)
      {
          mean.write(wp);
          var.write(wp);
      }
      void read(FILE* rp)
      {
          mean.read(rp);
          var.read(rp);
      }

      void writeRaw(FILE* wp)
      {
          mean.writeRaw(wp);
          var.writeRaw(wp);
      }
      void readRaw(FILE* rp)
      {
          mean.readRaw(rp);
          var.readRaw(rp);
      }


      SFloat getLnProb(SFloat f,int client_type)
      {
          SFloat v=f_minus(f,mean,client_type);//(x-mean)

          SFloat minus_two=float2SFloat(-2.0);
          SFloat zero=float2SFloat(0.0);
          SFloat one=float2SFloat(1.0);
          SFloat f2;
          SFloat factor;
         
          v=f_product(v,v,client_type);//(x-mean)^2
 
          //compute -2*var
          if(client_type==ALICE)
            f2=f_product(minus_two,var,client_type);
          else 
            f2=f_product(zero,var,client_type);
        
          v=f_division(v,f2,client_type);//get the exp part
          
          if(client_type==ALICE)
              factor=float2SFloat(2*PI);
          else factor=zero;
  
          factor=f_product(factor,var,client_type);//2*pi*var

          factor=f_ln(factor,client_type);
          
          if(client_type==ALICE)
              factor=f_division(factor,minus_two,client_type);
          else 
              factor=f_division(factor,zero,client_type);

          v=f_plus(factor,v,client_type);
          return v;
  
      }
} ;

class Secure_Sample {
  private:
    std::vector<SFloat> X;
    std::vector<SFloat> square_X; 
    SFloat square_sumX;//\sum X[i]^2
    SFloat Y;
    bool maintain_square;
  public:
    Secure_Sample(int nCols,bool maintain_square);
    Secure_Sample(){};
    ~Secure_Sample();
    void setColLen(int n);//set the length of feature
    void setX(SFloat v, int idx, int client_type);
    void setX(SFloat v, int idx);
    void setY(SFloat v);
    void setMaintainSquare(bool);
    bool isMainTainSquare();
    SFloat getX(int idx) { return X[idx]; }
    SFloat getY() { return Y; }
    int getColLen();
};

SFloat* getLine(FILE* rp, int nf, char* label);
SGaussian* createSGaussian(std::map<char*,int,cmpstr> lmap, std::map<char*,int,cmpstr> *fmap, int* lcount, int nf,int party_datasize,int local_datasize,int client_type);

SFloat s_sigmoid(SFloat v, int client_type);
/**
 Decide if a shared SFloat is positive or negative
 The function f_is_positive(SFloat,client_type,mpz), which requires gmp.h, does the same thing.
 */
SFloat s_sign(SFloat v, int client_type);
SFloat s_relu(SFloat v, int client_type);
//secure L2 distance computation. return L2^2
SFloat s_L2_distance(Sample *s1, Sample *s2, int nF, int client_type);
//another secure L2 distance computation. return L2^2
SFloat s_L2_distance(Secure_Sample *s1, Secure_Sample *s2, int nF, int client_type);
//SFloat sharedInt2SFloat(mpz_t s, int client_type);

#endif
