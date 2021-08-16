#include <gmp.h>
#include "smc.h"

#ifndef UTIL
#define UTIL
/**
 PriComp Settings.
 Note that the number of bits in Sign, Mantissa, Dimension, and Domain is not used in current version. These settings is defined in two_parties/protocol/cSFloat.h
 */
typedef struct
{
  char host[3][256];
  char offlineRNDDIR[256];
  int port[3];
  int64_t Sign; /**< The number of bits of the sign bit in secure float (SFloat)*/ 
  int64_t Exponent; /**< The number of bits of the sign bit in secure float (SFloat)*/ 
  int64_t Mantissa; /**< The number of bits of the sign bit in secure float (SFloat)*/ 
  int64_t Dimension;/**< The number of bits of the sign bit in secure float (SFloat)*/ 
  int64_t Domain;/**< The domain that the ooperations on the above three components in SFloat are performed on.*/ 
  int adderIDX; /**< Which adder is used. See the setting file.*/
  bool Offline_Commodity; /**< Whether the three-party (one commodity server) scalar product protocol use random bits generated offline or not.*/
} Protocol_ARG;

/**
 This structure is for the three-party (one commodity server) scalar product protocol.
 The suructure is for conducting multiple times of scalar products by calling the function one time, so some variables, for example ra and rb, are arrays.
 */
typedef struct
{
    ulong128* vec_pa; /**< Alice's private data, or Alice's encoded data Bob gets */
    ulong128* vec_pb; /**< Bob's private data, or Bob's encoded data Alice gets */
    ulong128* vec_Ra; /**< The random vector for Ra */
    ulong128* vec_Rb; /**< The random vector for Rb */
    ulong128* sp; /**< Auxiliary space for scalar product protocol. */
    ulong128* ra;/**< For ra */
    ulong128* rb;/**< For rb */
    ulong128* u;/**< For the result of one of the parties. */
    ulong128* v;/**< For the result of one of the parties. */
    mpz_t *domain;/**< The domain the (integer) protocol performed on */
} SPData;

/**
 Initialize an mpz array.
 @param vec The mpz array.
 @param dim The length if vec.
 */
void init_vector(mpz_t vec[], int dim);

/**
 Free the memory allocated by the mpz array.
 @param vec The mpz array.
 @param dim The length if vec.
 */
void clear_vector(mpz_t vec[], int dim);

/**
 Compute the number of bits required to store the value num.
 @param num The value.
 @return the number of bits.
 */
int binary_size(mpz_t num);

/**
 Count the number of columns in a dataset file.
 @param fname The dataset filename.
 @param sep The delimiter 
*/
int getColSize(const char* fname, const char* sep);

/** 
 Binary domain negation. That is, make 0 => 1 and 1 => 0.
 @param b The mpz value to be set.
 @param result The result.
*/
void mpz_binnegate(mpz_t b,mpz_t result);

/**
 Conduct an inner product between vec1 and vec2. This is not a secure protocol.
 @param result The result.
 @param vec1 vec1
 @param vec2 vec2
 @param len The length of vec1 (and vec2)
 */
void inner_product(mpz_t result, mpz_t* vec1, mpz_t* vec2, long128 len);//result, vec1, vec2, vec length


/**
 Conduct vector addition between vec1 and vec2. This is not a secure protocol.
 @param result The result vector.
 @param vec1 vec1
 @param vec2 vec2
 @param len The length of vec1 (and vec2)
 */
void vec_add(mpz_t result[], mpz_t vec1[], mpz_t vec2[], long128 len);

/**
 Conduct vector subtraction between vec1 and vec2. This is not a secure protocol.
 @param result The result vector.
 @param vec1 vec1
 @param vec2 vec2
 @param len The length of vec1 (and vec2)
 */
void vec_sub(mpz_t result[], mpz_t vec1[], mpz_t vec2[], long128 len);

/**
 Remove 0x0d and 0x0a from a string.
 @param s The string.
 @return The result string.
 */
char* rmNewLine(char* s);

/**
 Read the setting of PriComp from a configuration file.
 @param arg The structure variable to save the settings
 @param filename The filename of the configuration file.
 @return -1 if there is any error. Otherwise, return 0.
 */
int parse_ProtocolARG(Protocol_ARG* arg,char* filename);

/**
 Convert a 128-bit integer x to mpz_t data type.
 The max length of x in decimal is 255 if x>0, and is 254 if x<0.
 @param array An char array with length at least 256.
 @param x The 128-bit integer value
 @return A string s that represents x. The address of s is in the range of the array.
 */
char* long128_to_str(char array[],long128 x);

/**
 Use an 128-bit integer to set a mpz_t variable.
 @param op The mpz_t variable to be set
 @param x The 128-bit value
 */
void mpz_set_long128(mpz_t op,long128 x);

/**
 Get an 128-bit integer from an mpz_t variable.
 @param n The mpz_t variable.
 @return An 128-bit integer.
 */
long128 mpz_to_long128(mpz_t n);

/**
 Compute x mod y. The sign of the divisor is ignored so ithe return value is always non-negative. Note that x is an 128-bit integer and y is a mpz_t variable.
 @param x The dividend.
 @param y The divisor.
 @return An 128-bit integer z s.t. z = x%y
 */
ulong128 MOD(long128 x, mpz_t y);

/**
 Compute the inner product of two unsigned 128-bit integer arrays. This is not a secure protocol.
 @param x One of the two vectors
 @param y One of the two vectors
 @param dim The length of vector x (and y)
 @param domain The domain this operation is performed on.
 @return The result
 */
ulong128 inner_product(ulong128* x,ulong128* y,int dim,mpz_t domain);

/**
 Compute vector addition z=x+y s.t. z[i]=x[i]+y[i]. This is not a secure protocol.
 @param z The result
 @param x An unsigned 128-bit integer array.
 @param y An unsigned 128-bit integer array.
 @param dim The length of array x (and y)
 */
void vec_add(ulong128* z,ulong128* x,ulong128* y, int dim);

/**
 Set Protocol_ARG arg according to the arguments.
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param c_server A string of commodity server. The value should be "localhost" if the option "local_run" is et to true. Otherwise, give the ip of the commodity server.
 @param c_port The port the commodity server uses to communicate with the two parties.
 @param peer_host ip or "localhost" of the peer
 @param peer_host The port the peer host uses to perform the sevure protocols.
 @param myport The port the local host uses to perform the secure protocols.
 @param arg The Protocol_ARG variable.
 */
void init_net_params(int client_type,char* &c_server,int& c_port,char* &peer_host, int& peer_port,int& myport,Protocol_ARG &arg);

/**
 Read random bits from the commodity server.
 @param r A pointer of unsigned 128-bit integer, which points to ra (Alice) or rb (Bob).
 @param vector A point of unsigned 128-bit integer, whic points to Ra (Alice) or Rb (Bob).
 @param c_socket The socket used to read data from commodity server.
 @param dim The length of Ra (or Rb).
 @return The number of bytes read. The result should be (dim+1)*16 (128-bit integers includingtes
 */
int set_random_data(ulong128* r, ulong128* vector,int c_socket,int dim);

#endif
