#include "src/smc.h"
#include <cstdlib>

int main(int argc, char** argv)
{	
  //init 
  int client_type = (int)strtol(argv[argc-1], NULL, 10);
  SP_init(argv[1], client_type);
  mpz_t alice_wins;
  mpz_init(alice_wins);

  //get the local input
  double v=strtod(argv[2],NULL);
  //Secure floating-point vars
  SFloat alice_value, bob_value;

  //convert a floating-point var to secure floating-point var
  if(client_type==ALICE) {
    alice_value=fpShare(v,client_type);
    bob_value=fpShare(0,client_type);
  }
  else {
    alice_value=fpShare(0,client_type);
    bob_value=fpShare(v,client_type);
  }

  //check if alice_value > bob_value. The result is shared
  f_great_than(alice_value,bob_value,client_type,alice_wins);
  printf("Alice is richer: ");   
  mpz_out_str(stdout,10,alice_wins);
  printf("\n");

  //release
  mpz_clear(alice_wins);
  SP_clear(client_type);
  return 0;
}


