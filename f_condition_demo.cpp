#define _WITH_GETLINE

#include "src/two_parties/sp_init.h"
#include "src/util.h"
#include "src/time_eval.h"
#include "src/debug_util.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

void usage(const char* fname)
{
	fprintf(stderr,"Demo secure conditional expressions on floating-point values in PriComp. Similar functions for integer values can be found in s_fixnum.h\n");
  fprintf(stderr,"This program demonstrate how to check the (in)equality relationships between two values hold by Alice (v1) and Bob (v2) without revealing them to the peer party. \n");
  printf("To make the output easy to read, the result is v1 if the condition is true, otherwise v2.\n");
	fprintf(stderr,"usage: %s <conf_file> <value> <client_type>\n", fname);
	fprintf(stderr,"\t<conf_file> is the configuration file of Pricomp. Check the cfg files in the conf directory.");
	fprintf(stderr,"\t<value> is the private inout of the party.\n");
	fprintf(stderr,"\t<client_type>: 1 (Alice) or 2 (Bob).\n\n");
	fprintf(stderr,"Notice that ALICE HAVE TO RUN BEFORE BOB.\n");
}

void print_result(const char *msg, mpz_t v)
{
  printf("%s\n", msg);
  printf("Additively shared binary output: ");
  mpz_out_str(stdout,10,v);
  printf("\n");
}

int main(int argc, char** argv)
{	
	if(argc!=4) { usage(argv[0]); return 1; }

	int client_type = (int)strtol(argv[3], NULL, 10);
	double value = strtod(argv[2], NULL);

	mpz_t domain;
	mpz_t result;
	mpz_init(domain);
	mpz_init(result);
	mpz_set_ui(domain,10);

	printf("Start to init secure computation\n");
	SP_init(argv[1], client_type);
	printf("successfully init secure computation\n");

#ifdef EVAL_TIME
	struct timeval t_sum;
	struct timeval t_sub;
	t_sum.tv_sec=0;
	t_sum.tv_usec=0;
#endif

	SFloat v1;
	SFloat v2;
  SFloat output;

  if(client_type==ALICE) v1=float2SFloat(value);//Alice's case: v1=value, and the value v2's components (sign bit, exponent, and mantissa) are all zero.
  else v2=float2SFloat(value);

  /* you can also use fpShare to convert an additively-shared floating point value to SFloat as follows:
   */
  /*
  if(client_type==ALICE) 
  {
    v1=fpShare(value,client_type); 
    v2=fpShare(0,client_type); 
  }
  else //Bob
  {
    v1=fpShare(0,client_type);  
    v2=fpShare(value,client_type); 
  }
  */
 
  //Check the relationships between v1 and v2 and print results
  //Notice that the value of the binary variable result (1: true, 0: false) is additively shared between the two parties.
  //Function fpMerge is usually used to debug or to get the final result only.
  //Only use fpMerge to merge the component-wisedly shared value if revealing that value does not leak other information.
  //v1 > 0?
  f_is_positive(v1, client_type, result);
  print_result("The value hold by Alice is positive?",result);
  //output= (v1 > 0)? v1: v2
  output=f_if_then_else(result,v1,v2,client_type);
  printf("result is %f\n",fpMerge(output,client_type));

  //v2 < 0?
  f_is_negative(v2, client_type, result);
  print_result("The value hold by Bob is negative?",result);
  //output= (v2 < 0)? v1: v2
  output=f_if_then_else(result,v1,v2,client_type);
  printf("result is %f\n",fpMerge(output,client_type));


  // v1 is zero?
  f_is_zero(v1,client_type,result);
  print_result("The value hold by Alice is zero?",result);
  //output= (v1 == 0)? v1: v2
  output=f_if_then_else(result,v1,v2,client_type);
  printf("result is %f\n",fpMerge(output,client_type));


  //v1 > v2?
  f_great_than(v1,v2,client_type,result);
  print_result("The value hold by Alice is larger than that hold by Bob?",result);
  //output= (v1 > v2)? v1: v2
  output=f_if_then_else(result,v1,v2,client_type);
  printf("result is %f\n",fpMerge(output,client_type));

  //v1 < v2?
  f_less_than(v1,v2,client_type,result);
  print_result("The value hold by Alice is less than that hold by Bob?",result);
  //output= (v1 < v2)? v1: v2
  output=f_if_then_else(result,v1,v2,client_type);
  printf("result is %f\n",fpMerge(output,client_type));

  //v1 == v2?
  f_equal(v1,v2,client_type,result);
  print_result("The value hold by Alice is equal to that hold by Bob?",result);
  // output = (v1 == v2) ? v1 : v2
  output=f_if_then_else(result,v1,v2,client_type);
  printf("result is %f\n",fpMerge(output,client_type));

	mpz_clear(domain);
	mpz_clear(result);
	SP_clear(client_type);
	return 0;
}


