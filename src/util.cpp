#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <bitset>
#include <iostream>
#include <sys/time.h>
#include "util.h"
#include "gmp.h"
#include "communicate.h"
#define MAX_LINE_LEN 1024

/** 
 print a mpz vector to the specified device.
 */
void print_vec(FILE* dev,mpz_t vector[], long128 dim)
{
  int i;
  fprintf(dev,"[");
  for(i = 0; i<dim; i++){
    mpz_out_str(stderr, 10, vector[i]);
    if((i % 10 == 0) && (i!=0)) fprintf(dev,",\n");
    else if(i+1<dim) fprintf(dev,", ");
  }
  fprintf(dev,"]\n");
}

/**
 * print a mpz value to a specified device. 
 */
void print_val(FILE* dev,mpz_t val)
{
  mpz_out_str(dev, 10, val);
}

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
    fprintf(stderr,"Failed to read data from file %s!\n",file);
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

int binary_size(mpz_t num)
{
  int bits = 0;
  mpz_sub_ui(num, num, 1);
  bits = mpz_sizeinbase(num, 2);
  mpz_add_ui(num, num, 1);

  return bits;
}


void mpz_binnegate(mpz_t result,mpz_t b)
{
  mpz_add_ui(result,b,1);
  mpz_mod_ui(result,result,2);
}

void clear_vector(mpz_t vec[], int dim)
{
  int i;
  for(i = 0; i < dim; i++) mpz_clear(vec[i]);
}
void init_vector(mpz_t vec[], int dim)
{
  int i;
  for(i = 0; i < dim; i++)
    mpz_init(vec[i]);
}

// result = vec1 * vec2 (* is scalar product)
void inner_product(mpz_t result, mpz_t vector1[], mpz_t vector2[], long128 dim)
{
  int i; 

  mpz_set_ui(result, 0);
  for(i = 0; i < dim; i++)
    mpz_addmul(result, vector1[i], vector2[i]);

} 

// add = vec1 + vec2
void vec_add(mpz_t sum[], mpz_t vector1[], mpz_t vector2[], long128 dim)
{
  int i;
  for(i = 0; i < dim; i++)
    mpz_add(sum[i], vector1[i], vector2[i]);
}

// dif = vec1 - vec2
void vec_sub(mpz_t dif[], mpz_t vector1[], mpz_t vector2[], long128 dim)
{
  int i;
  for(i = 0; i < dim; i++)
    mpz_sub(dif[i], vector1[i], vector2[i]);
}

ulong128 inner_product(ulong128* x,ulong128* y,int dim,mpz_t domain)
{
  mpz_t* px=(mpz_t*)malloc(sizeof(mpz_t)*dim);
  mpz_t* py=(mpz_t*)malloc(sizeof(mpz_t)*dim);
  mpz_t result;
  ulong128 s=0;

  init_vector(px,dim);
  init_vector(py,dim);
  mpz_init(result);
  for(int i=0;i<dim;i++) 
  {
    mpz_set_long128(px[i],x[i]);
    mpz_set_long128(py[i],y[i]);
  }

  inner_product(result,px,py,dim);

  mpz_mod(result,result,domain);

  s=mpz_to_long128(result);

  mpz_clear(result);
  clear_vector(px,dim);
  clear_vector(py,dim);
  free(px);
  free(py);
  return s;
}

void parse_err(char* line)
{
  fprintf(stderr,"Syntax error or unrecognized line in configuration file: \"%s\"\n",line);
}

bool parse_floating_point(char* args,int64_t& param,bool& setting)
{
  if(!args) setting=false;
  else param=(int64_t)strtoll(args,NULL,10);

  return setting;
}

char* rmNewLine(char* line)
{
  char* tmps=strchr(line,0x0a);
  if(tmps) *tmps='\0';
  tmps=strchr(line,0x0d);
  if(tmps) *tmps='\0';
  return line;
}


/** 
 Append a mpz_t value to a string, return address is at the end position (that is, the location of '\0') of the string.
 */
char* append_mpz(char* msg,mpz_t r)
{
  char buf[512]; //bzero(buf,512);
  char* tmp=mpz_get_str(buf,10,r);
  char* ptr=msg;
  int len=strlen(buf);
  for(int j=0;j<len;j++,ptr++) *ptr=*(tmp+j);//copy char by char
  *ptr='\0';
  return ptr;
}

/** 
 Append a mpz_t array to a string that elements are separated by comma
 */
char* append_mpz_vec(char* msg,mpz_t R[], long128 dim)
{
  char* ptr=msg;
  for(int i=0;i<dim;i++)
  {
    ptr=append_mpz(ptr,R[i]);
    *ptr=',';ptr++;
  }
  ptr--; //remove last comma
  *ptr='\0';//remove the last comma
  return ptr;
}

char* long128_to_str(char ptr[],long128 x)
{
  char digit[10]={'0','1','2','3','4','5','6','7','8','9'};
  ptr[255]='\0';
  int idx=254;
  bool isnegative=false;

  long128 v=x;
  v>>=127;

  char tmps[129];
  tmps[128]='\0';

  if(v!=0)
  {
    isnegative=true;
    x*=-1;
  }

  if(x==0) {*(ptr+254)='0'; return ptr+idx;}

  //x is a positive value now
  for(idx=254;idx>0&&x!=0;idx--) 
  {
    int d=x%10;

    *(ptr+idx)=digit[d];

    x/=10;

  }
  if(isnegative) { *(ptr+idx)='-'; idx--;}
  idx++;

  return ptr+idx;
}

void mpz_set_long128(mpz_t op,long128 x)
{
  char buf[256];
  char* ptr=long128_to_str(buf,x);
  mpz_set_str(op,ptr,10);
}

long128 mpz_to_long128(mpz_t n)
{
  mpz_t q,r;
  mpz_init(q); 
  mpz_init(r);
  long128 v=0;
  long128 d=1;
  bool isnegative=false;

  mpz_set(q,n);//copy
  if(mpz_sgn(n)<0)//is negative value
  {
    isnegative=true;
    mpz_neg(q,n);//let q>0
  }

  while(mpz_cmp_ui(q,0)>0)
  {
    v+=d*mpz_mod_ui(r,q,10);
    mpz_fdiv_q_ui(r,q,10);
    mpz_set(q,r);//copy
    d*=10;
  }

  if(isnegative) v*=-1;

  mpz_clear(q); mpz_clear(r);

  return v;
}

ulong128 MOD(long128 x, mpz_t y)
{
  int idx=0;
  mpz_t v;
  mpz_init(v);

  mpz_set_long128(v,x);

  mpz_mod(v,v,y);

  x=mpz_to_long128(v);
  mpz_clear(v);

  return x;
}

void vec_add(ulong128* result,ulong128* x,ulong128* y, int dim)
{
  for(int i=0;i<dim;i++) 
  {
    result[i]=x[i]+y[i];
  }
}

//The original dimension is D, and 
//The dimension of sub matrices M1, M2 is d.
//D and d are both assumed to be 2^n for some int n
int* MatrixSub(int d, int D, int* M1, int* M2)
{
  int size=d*d;
  int* M=(int*)malloc(sizeof(int)*size);
  bzero(M,sizeof(int)*size);
  for(int i=0;i<d;i++)
    for(int j=0;j<d;j++)
      M[i*d+j]=M1[i*D+j]-M2[i*D+j];
  return M;
}

//The original dimension is D, and 
//The dimension of sub matrices M1, M2 is d.
//D and d are both assumed to be 2^n for some int n
int* MatrixAdd(int d, int D, int* M1, int* M2)
{
  int size=d*d;
  int* M=(int*)malloc(sizeof(int)*size);
  bzero(M,sizeof(int)*size);
  for(int i=0;i<d;i++)
    for(int j=0;j<d;j++)
      M[i*d+j]=M1[i*D+j]+M2[i*D+j];
  return M;
}

//do matrix row operation
//Let the row of matrix M is M1,...,Mn, the resulting Mr=a*Ms+Mr
int* RowOperation(int r, int s, int a, int nCol, int * M)
{
  for(int i=0;i<nCol;i++)
    M[r*nCol+i]+=a*M[s*nCol+i];
  return M;
}

//generate a dim by dim symmetrix matrix M that M=M^-1
//and for all i,j, M_ij is integer.
//Method: generate A with matrix row operations from I,
//and let M=A*A^T. Since |A|=1, M=M^-1.
//parameters:
//  dim: dimension
//  rowop: (dim-1)*(dim-1) integers used to generate A
int* genMatrix(int d,int *rowop)
{
  int size=d*d;
  int* A=(int*)malloc(sizeof(int)*size);

  bzero(A,sizeof(int)*size);

  for(int i=0;i<d;i++) 
  { 
    A[i*d+i]=1;//create I
  }

  //do (dim-1)*(dim-1) row operations
  int c=0;
  int i=0;
  int j=0;
  for(;i<d;i++)
  {
    for (j=0;j<d;j++)
    {
      if(i==j) continue;
      A=RowOperation(j,i,rowop[c],d,A);
      c++;
    }
  }

  return A;
}

int parse_ProtocolARG(Protocol_ARG* arg,char* filename)
{
  FILE* rp=fopen(filename,"r");
  if(!rp) { fprintf(stderr,"Error in reading configuration file %s!\n",filename); return -1;}

  char* buf=(char*)malloc(sizeof(char)*MAX_LINE_LEN);
  if(!buf) { fprintf(stderr,"Failed to allocate memory of reading buffer!\n"); return -1;}
  char* line=(char*)malloc(sizeof(char)*MAX_LINE_LEN);
  if(!line) { fprintf(stderr,"Failed to allocate memory of backup reading buffer!\n"); return -1;}

  //init
  bool islocal=false;
  for(int i=0;i<3;i++)
  {
    strcpy(arg->host[i],"localhost"); arg->port[i]=0; 
  }
  arg->Sign=-1; arg->Exponent=-1; arg->Mantissa=-1;
  bool setting=true;

  //read from the setting file
  while(fgets(buf,MAX_LINE_LEN-1,rp))
  {
    if(strlen(buf)<=3) continue;//skip null lines
    if(buf[0]=='#') continue; //skip comment lines

    const char* sep="\n\t,= ";
    strcpy(line,buf);
    char* tmps=strtok(buf,sep);

    //check if successfully get the setting
    if(!tmps) { parse_err(line); }
    //run_type
    if(strcmp(tmps,"local_run")==0)//local run 
    {
      char* local_run=strtok(NULL,sep);

      if(local_run)
      {
        if(strstr(local_run,"true"))
        {
          islocal=true;
        }
        else if(strcasecmp(local_run,"false")!=0)
        {
          parse_err(line);
        }
      }

    }
    else if(strcmp(tmps,"adder")==0)//set adder
    {
      tmps=strtok(NULL,sep);
      arg->adderIDX=(int)strtol(tmps,NULL,10);

    }
    else if(strcmp(tmps,"$host0")==0)
    {
      if(islocal==false)//not run on local host
      {
        char* tmps=NULL;
        strcpy(buf,line);//reset buf from backup
        tmps=strtok(buf,":");
        tmps=strtok(NULL,":");
        char *host_str=strtok(tmps,"[\"], \n");
        for(int i=0;i<3;i++)
        {
          if(!host_str) { setting=false; break; }
          strcpy(arg->host[i],host_str);
          host_str=strtok(NULL,"[\"], \n");
        }
      }
    }
    else if(strcmp(tmps,"$port0")==0)
    {
      strcpy(buf,line);//reset buf from backup
      char* ptr=strchr(buf,'='); ptr++;
      tmps=strtok(ptr,",\n");//start to get port no.
      for(int i=0;tmps;i++)
      {
        if(!tmps) { setting=false; break; }
        arg->port[i]=(int)strtol(tmps,NULL,10);
        tmps=strtok(NULL,",\n");
      }
    }
    else if ((strcasecmp(tmps,"offline_commodity")==0))
    {
      tmps=strtok(NULL,sep);
      if(strstr(tmps,"true")) arg->Offline_Commodity=true;
      else arg->Offline_Commodity=false;
    }
    else if(strcasecmp(tmps,"offlineRNDDIR")==0)
    {
      tmps=strtok(NULL,sep);
      strcpy(arg->offlineRNDDIR,tmps);
    }
    else if ((strcasecmp(tmps,"Sign")==0))
      parse_floating_point(strtok(NULL,sep),arg->Sign,setting);
    else if(strcasecmp(tmps,"Exponent")==0)
      parse_floating_point(strtok(NULL,sep),arg->Exponent,setting);
    else if(strcasecmp(tmps,"Mantissa")==0)
      parse_floating_point(strtok(NULL,sep),arg->Mantissa,setting);
    else if(strcasecmp(tmps,"Dimension")==0)
      parse_floating_point(strtok(NULL,sep),arg->Dimension,setting);
    else if(strcasecmp(tmps,"Domain")==0)
      parse_floating_point(strtok(NULL,sep),arg->Domain,setting);
    else { parse_err(line);}
  }

  free(line);
  free(buf);
  fclose(rp);

  return 0;

}

void init_net_params(int client_type,char* &c_server,int& c_port,char* &peer_host, int& peer_port,int& myport,Protocol_ARG &arg)
{
  c_server=arg.host[0];//commodity server hostname
  c_port=arg.port[0];//commodity server port
  myport=arg.port[client_type];//this party will listen on this port

  if(client_type==1) //Alice
  {  
    peer_host=arg.host[2]; 
    peer_port=arg.port[2]; 
  }
  else  // Bob
  {  
    peer_host=arg.host[1]; 
    peer_port=arg.port[1]; 
  }
}

int set_random_data(ulong128* r, ulong128* vector,int c_socket,int dim)
{
  //read random data, long version
  int len=0;
  int v=0;
  int i=0;
  char* tmpl=(char*)vector;

  len=(dim+1)*sizeof(long128);

  v=protocol_read(c_socket,tmpl,len);
  if(v<0) { fprintf(stderr,"Error on getting random bits from the commodity server.\n"); return -1; }
  *r=vector[dim];//the last element is ra/rb 

  return v;
}

