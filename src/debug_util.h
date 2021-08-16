/**
 This header file defines some functions that may be useful to debug your protocols.
 * */
#ifndef DEBUG_UTIL
#define DEBUG_UTIL
#include "smc.h"
#include "communicate.h"
#include "two_parties/two_parties.h"
#include "two_parties/cSFloat.h"
#include "util.h"

/**
 * Print the value of a given SFloat.
 * Notice that SFloat is usually a component-wisedshared floating-point value, so the value may be meaningless
 @param s The SFloat variable to be printed.
 */
void print(SFloat s);

/**
 * Add two SFloat variables f1 and f2 compoment by component. If a floating-point variable f is shared among two parties using f1 and f2, the result is f.
 @param f1 One of the two SFloat variables
 @param f2 One of the two SFloat variables
 @return A floating-variable f that is component-wised shared by f1 and f2.
 */
double AddSFloat(SFloat f1,SFloat f2);

/**
 Write a SFloat variable using the socket.
 @param socket The socket to write the data.
 @param v The SFloat variable to be sent using the socket.
 */
void sendSFloat(int socket, SFloat v);

/**
 Get a SFloat variable using the socket.
 @param socket The socket to read the data from.
 @param v The pointer of a SFloat variable to get the data.
 @return The number of bytes read from the socket.
 */
int getSFloat(int socket, SFloat *v);
/**
 Merge two additively shared integers v1 and v2. The result will be (v1+v2)%domain.
 @param v The value hold by this party.
@param domain The domain of the integer.
 @param client_type The type of this party:
  - ALICE (1)
  - BOB (2)
 */
int sMerge(mpz_t v, mpz_t domain, int client_type);
/**
 Merge a shared SFloat. This function can be used to debug.
 @param v The component-wised shared SFloat hold by the party.
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @return Both the two parties will get the value of the shared value.
 */
double fpMerge(SFloat v,int client_type);

#endif
