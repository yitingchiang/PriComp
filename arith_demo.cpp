#define _WITH_GETLINE

#include "src/two_parties/sp_init.h"
#include "src/util.h"
#include "src/time_eval.h"
#include "src/debug_util.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#ifdef EVAL_TIME
struct timeval Stamp[2];
#endif

void usage(const char* fname)
{
	fprintf(stderr,"Demo how to use PriComp to perform the four arithmetic operations.\n");
	fprintf(stderr,"usage: %s <conf_file> <function> <value> <round> <client_type>\n", fname);
	fprintf(stderr,"\t<conf_file> is the configuration file of Pricomp. Check the cfg files in the conf directory.");
	fprintf(stderr,"\t<function>:\n");
	fprintf(stderr,"\t\t0: addition\n");
	fprintf(stderr,"\t\t1: subtraction\n");
	fprintf(stderr,"\t\t2: multiplication\n");
	fprintf(stderr,"\t\t3: division\n");
	fprintf(stderr,"\t<value> is the private inout of the party.\n");
	fprintf(stderr,"\t<round> is how many times the specified operation are performed. This sis for performance evaluation.\n");
	fprintf(stderr,"\t<client_type>: 1 (Alice) or 2 (Bob).\n\n");
	fprintf(stderr,"Notice that ALICE HAVE TO RUN BEFORE BOB.\n");
}

int binary_s(int num)
{
	int bits = 1;
	int n = 2;
	while(num / n)
	{
		bits++;
		n *= 2;
	}
	return bits;
}

int main(int argc, char** argv)
{	
	if(argc!=6) { usage(argv[0]); return 1; }
	int client_type = (int)strtol(argv[5], NULL, 10);
	int K = (int)strtol(argv[4], NULL, 10);

	mpz_t domain;
	mpz_t result;
	mpz_init(domain);
	mpz_init(result);
	mpz_set_ui(domain,10);

	printf("Start to init secure computation\n");
	SP_init(argv[1], client_type);
	printf("successfully init secure computation\n");

	srandom(time(NULL));

	int fun_idx=(int)strtol(argv[2],NULL,10);
	double v=strtod(argv[3],NULL);

	SFloat(*ari_func[4])(SFloat,SFloat,int);
	ari_func[0]=f_plus;
	ari_func[1]=f_minus;
	ari_func[2]=f_product;
	ari_func[3]=f_division;
	const char* fun_name[4]={"addition","subtraction","multiplication","division"};

	SFloat d;
	SFloat f1;
	SFloat f2;
#ifdef EVAL_TIME
	struct timeval t_sum;
	struct timeval t_sub;
	t_sum.tv_sec=0;
	t_sum.tv_usec=0;
#endif
	if(client_type==ALICE) 
	{
		printf("Run as Alice with local input %f\n",v);
		f1=float2SFloat(v);
	}
	else
	{
		printf("Run as Bob with local input %f\n",v);
		f2=float2SFloat(v);
	}
	printf("Start to run %s protocols for %d times:\n",fun_name[fun_idx],K);
	for(int i=0;i<K;i++)
	{
		printf("%d-th iteration.\n",i+1);
#ifdef EVAL_TIME
		time_set(Stamp,0);
#endif
		d=ari_func[fun_idx](f1,f2,client_type);
#ifdef EVAL_TIME
		time_set(Stamp,1);
		time_print(Stamp,0,1);
#endif
	}
#ifdef DEBUG
	printf("ans: %f\n",fpMerge(d,client_type));
#else
	printf("Local shared SFloat:\n");
        d.to_s();
#endif
	mpz_clear(domain);
	mpz_clear(result);
	SP_clear(client_type);
	return 0;
}


