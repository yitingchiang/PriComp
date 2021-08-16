#define _WITH_GETLINE
//do time evaluation and print debug message
#include "src/util.h"
#include "src/two_parties/two_parties.h"
#include "src/two_parties/sp_init.h"
#include "src/debug_util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "src/time_eval.h"

void usage(const char* fname)
{
    fprintf(stderr,"A demonstration of using Pricomp to compute some statistics.\n");
    fprintf(stderr,"This program reads data file a.test (ALICE) or b.test (Bob) from current directory and performs the specified operation\n");
    fprintf(stderr,"The format of a.test and b.test is one integer N followed by N values.\n");
    fprintf(stderr,"usage: %s <conf_file> <data_size> <op> <client_type>\n",fname);
    fprintf(stderr,"\t<conf_file> is the configuration file of Pricomp. Check the cfg files in the conf directory.");
    fprintf(stderr,"\t<op>:\n");
    fprintf(stderr,"\t\t 1: average\n");
    fprintf(stderr,"\t\t 2: variance\n");
    fprintf(stderr,"\t\t 3: median\n");
    fprintf(stderr,"\t<conf_file> is the setting file. You can fild some examples in the directory conf.\n");
    fprintf(stderr,"\t<data_size> is the value N in the first line of PEER's input file. For example, in Alice's case, this value is the first line in b.test.\n");
    fprintf(stderr,"\t<client_type>: 1 (Alice) or 2 (Bob).\n\n");
  fprintf(stderr,"Notice that ALICE HAVE TO RUN BEFORE BOB.\n");

}

int main(int argc, char** argv)
{	

    if(argc!=5) {usage(argv[0]); return 0;}
    int party_dim=(int)strtol(argv[2],NULL,10);
    int op = (int)strtol(argv[3], NULL, 10);
    int client_type = (int)strtol(argv[4], NULL, 10);
#ifdef EVAL_TIME    
    struct timeval t_main[6];
#endif
    SP_init(argv[1], client_type);
    char fname[16];
    if(client_type==ALICE) strcpy(fname,"a.test");
    else strcpy(fname, "b.test");
    
    FILE *fp=fopen(fname,"r");
    if(!fp) {fprintf(stderr,"cannot open file %s!\n",fname); return 1;}

    int dim;

    fscanf(fp,"%d",&dim);

    double* data=new double[dim];

    SFloat* S_data=new SFloat[dim+party_dim];
    SFloat result;
    float f;

    if(client_type==ALICE)
    {
      for(int i=0;i<dim;i++) { fscanf(fp,"%lf",data+i); S_data[i]=float2SFloat(data[i]);}
      for(int i=dim;i<dim+party_dim;i++) S_data[i]=float2SFloat(0.0);
    }
    else if(client_type==BOB)
    {
      for(int i=0;i<party_dim;i++) {  S_data[i]=float2SFloat(0.0);}
      for(int i=party_dim;i<dim+party_dim;i++) { fscanf(fp,"%lf",data+i-party_dim); S_data[i]=float2SFloat(data[i-party_dim]); }
    }
    else return 1;

    fclose(fp);

    const char ops[][8]={"average","var","median"};
    printf("compute %s\n",ops[op-1]);
#ifdef EVAL_TIME    
    time_set(t_main,0);
#endif
    switch(op)
    {
      case 1: //mean
        result = f_average_2(party_dim,data,dim,client_type);
        break;
      case 2: //var
        result=f_variance_2(party_dim,data,dim,client_type);
        break;
      case 3: //median
        if(client_type==ALICE)
           result=f_median(S_data,S_data+dim,dim,party_dim,client_type);
        else 
           result=f_median(S_data,S_data+party_dim,party_dim,dim,client_type);
        break;
        default:
        break;
    }
#ifdef EVAL_TIME    
    time_set(t_main,1);
#endif
    double v_recovered=fpMerge(result,client_type);
    
    printf("result: %f\n",v_recovered);
    delete[] S_data;
    delete[] data;

    printf("\n");

#ifdef EVAL_TIME    
    time_print(t_main,0,1);
#endif
    SP_clear(client_type);

}


