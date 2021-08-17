/**
 This header file defines procedures to initialize and to end PriComp. The elemental procedure of PriComp, scalar product, is also defined here.
 */
#include <gmp.h>
#ifndef SP_INIT
#define SP_INIT
/**
 Defined in client.cpp, this socket is for the scalar product protocol that requires a commodity server. 
 */
extern int SP_TCP_c_socket;
/**
 Defined in client.cpp, this socket is for send/recv between the two parties.
 */
extern int SP_p_socket;
/**
 Defined in lcient.cpp. This variable is to count how many times the scalar product protocol is performed.
 This variable is for debugging only.
 */
extern long long int t_count;

/**
 Free the resources allocated in SP_init()
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 */
void SP_clear(int client_type);
/**
 Initialize network sockets for scalar product protocol. In addition, tables to compute exp and log are also initialized here.
 For offline commodity server mode, this function also init the random bits.
 To do: increase the robustness of the offline random bits.
 @param file The setting file.
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @return An integer to indicate if the initialization successes.
   - 0: Success
   - Otherwise: Fail
 */
int  SP_init(char* file,int& client_type);
/**
 */
int  scalar_product(mpz_t *px, int client_type, int dim, mpz_t dom,mpz_t result);

int multi_scalar_product(mpz_t *px, int client_type, int nSC, int* dim, mpz_t dom,mpz_t* result);
#endif
