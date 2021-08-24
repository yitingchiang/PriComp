/**
 s_fixnum.h defines secure functions for two-party secure integer computation. The input and output are additively shared and should be in a given domain. For example, to compute  5 times -3 under domain 7, the two inputs should be 5 and 4 (i.e. 7 + -3). Using negative values results in wrong answer.
 */

#ifndef _SFIXNUM_H_
#define _SFIXNUM_H_

#include "gmp.h"

/**
 Compute x+y where x and y are both additive shared integers. The result is also additive shared. If the domain is K, then the result each party gets will smaller than K. But the sum of then may be equal (if the result is zero) or larger than K.
 @param x an additive shared integer.
 @param y an additive shared integer.
 @param domain the domain this operation is performed on.
 @param result result = (x+y) mod domain. The result is shared between the two parties.
 */
void s_plus(mpz_t x, mpz_t y, mpz_t domain, mpz_t result);
/**
 Compute x-y where x and y are both additive shared integers. The result is also additive shared. If the domain is K, then the result each party gets will smaller than K. But the sum of then may be equal (if the result is zero) or larger than K.
 @param x an additive shared integer.
 @param y an additive shared integer.
 @param domain the domain this operation is performed on.
 @param result result = (x-y) mod domain. The result is shared between the two parties.
 */
void s_minus(mpz_t x, mpz_t y, mpz_t domain, mpz_t result);
/**
 Compute x*y where x and y are both additive shared integers. The result is also additive shared. If the domain is K, then the result each party gets will smaller than K. But the sum of then may be equal (if the result is zero) or larger than K.
 @param x an additive shared integer.
 @param y an additive shared integer.
 @param domain the domain this operation is performed on.
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result result = (x*y) mod domain. The result is shared between the two parties.
 */
void s_product(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result);
/**
 Compute x^2 where x is an additive shared integer. The result is also additive shared. If the domain is K, then the result each party gets will smaller than K. But the sum of then may be equal (if the result is zero) or larger than K.
 @param x an additive shared integer.
 @param domain the domain this operation is performed on.
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result result = (x^2) mod domain. The result is shared between the two parties.
 */
void s_square(mpz_t x, mpz_t domain, int client_type, mpz_t result);
/**
 Set the result according to a binary flag. That is
 - result = x if flag = 1 (true)
 - result = y if flag = 0 (false)

 Only binary flags can get the meaningful result.
 @param flag a binary additive-shared integer.
 @param x the first candidate result which is additive shared integer.
 @param y the second candidate result which is additive shared integer.
 @param domain the domain this operation is performed on.
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result an additive shared integer. 
 */
void s_if_then_else(mpz_t flag, mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result);
/**
 Set multiple results according to one binary flag under the same domain. That is, for all i,
 - result[i] = x[i] if flag = 1 (true)
 - result[i] = y[i] if flag = 0 (false)
 
 This function sets result[0], result[1],...,result[nResults-1] to x[0](or y[0]), x[1] (or y[1]), ..., x[nResults-1](or y[Results-1]. Only binary flags can get meaningful results.
 @param flag a binary additive-shared integer.
 @param x the first candidate result which is additive shared integer.
 @param y the second candidate result which is additive shared integer.
 @param domain the domain this operation is performed on.
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result an additive shared integer. 
 @param nResults how many results are there to be set. 
 */
void s_multiple_if_then_else(mpz_t b, mpz_t* x, mpz_t* y, mpz_t domain, int client_type, mpz_t* result, int nResults);

/**
 Check if x is a negative value under the domain. 
 @param x an additive shared integer.
 @param domain the domain this operation is performed on.
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result result = 1 (x is negative) or 0 (x is not negative). The result is additive shared between the two parties.
 */
void s_negative(mpz_t x, mpz_t domain, int client_type, mpz_t result);
/**
 Check if x < y is satisfied under the domain. 
 @param x an additive shared integer.
 @param y an additive shared integer.
 @param domain the domain this operation is performed on.
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result result = 1 (x < y) or 0 (x >= y). the result is additive shared between the two parties.
 */
void s_less_than(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result);
/**
 Check if x > y is satisfied under the domain. 
 @param x an additive shared integer.
 @param y an additive shared integer.
 @param domain the domain this operation is performed on.
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result result = 1 (x > y) or 0 (x <= y). The result is additive shared between the two parties.
 */
void s_great_than(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result);
/**
 Compute -x under the given domain. 
 @param x an additive shared integer.
 @param domain the domain this operation is performed on. The domain should be 2 (as a mpz_t int).
 @param result result = -x (or -x + domain). The result is additive shared between the two parties.
 */
void s_neg(mpz_t x, mpz_t domain, mpz_t result);
/**
 Compute logic x or y under the given domain. Value 1 represents true; 0 represents false.
 @param x an additive shared integer.
 @param y an additive shared integer.
 @param domain the domain this operation is performed on. The domain should be 2 (as a mpz_t int).
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result result = 1 (x OR y is true) or 0 (x OR y is false). The result is additive shared between the two parties.
 */
void s_or(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result);
/**
 Compute logic x and y under the given domain. Value 1 represents true; 0 represents false.
 @param x an additive shared integer.
 @param y an additive shared integer.
 @param domain the domain this operation is performed on. The domain should be 2 (as a mpz_t int).
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result result = 1 (x AND y is true) or 0 (x AND y is false). The result is additive shared between the two parties.
 */
void s_and(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result);
/**
 Compute logic not x under the given domain. Value 1 represents true; 0 represents false. If the input x is not 1 or 0, the re
 @param x an additive shared integer.
 @param domain the domain this operation is performed on. The domain should be 2 (as a mpz_t int).
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result result = 1 (x is true) or 0 (x is false). The result is additive shared between the two parties.
 */
void s_not(mpz_t x, mpz_t domain, int client_type, mpz_t result);
/**
 Check if x is equal to y under the domain. 
 @param x an additive shared integer.
 @param y an additive shared integer.
 @param domain the domain this operation is performed on.
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result result = 1 (x is equal to y) or 0 (x is not equal to y). The result is additive shared between the two parties.
 */
void s_equal(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result);
/**
 Compute integer division of x/y under the given domain.
 @param x an additive shared integer.
 @param y an additive shared integer.
 @param domain the domain this operation is performed on. 
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result result = 1 (x AND y is true) or 0 (x AND y is false). The result is additive shared between the two parties.
 */
void s_division(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result);
/**
 Compute remainder of integer x/y under the given domain. That is, x mod y. 
 @param x an additive shared integer.
 @param y a positive additive-shared integer.
 @param domain the domain this operation is performed on. The domain should be 2 (as a mpz_t int).
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result result = 1 (x AND y is true) or 0 (x AND y is false). The result is additive shared between the two parties.
 */
void s_remainder(mpz_t x, mpz_t y, mpz_t domain, int client_type, mpz_t result);
/**
 Compute base^exp under the given domain.
 @param base the base. an additive shared integer.
 @param exp the exponent. exp must be a positive additive-shared integer.
 @param domain the domain this operation is performed on. 
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result result = base^exp. The result is additive shared between the two parties.
 */
void s_power(mpz_t base, mpz_t exp, mpz_t domain, int client_type, mpz_t result);
/**
 Compute x << n under the given domain.
 @param x an additive shared integer.
 @param n must be a positive additive-shared integer.
 @param domain the domain this operation is performed on. 
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result result = x * 2^n. The result is additive shared between the two parties.
 */
void s_shift_left(mpz_t x, mpz_t n, mpz_t domain, int client_type, mpz_t result);
/**
 Compute x >> n under the given domain.
 @param x an additive shared integer.
 @param n an additive shared integer.
 @param domain the domain this operation is performed on. 
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result result = x * 2^n. The result is additive shared between the two parties.
 */
void s_shift_right(mpz_t x, mpz_t n, mpz_t domain, int client_type, mpz_t result);
/**
 A ripple-cary adder that converts an additive shared x into the same value in binary domain. The output of the two parties get can perform XOR to get the result.
 @param x the input in the given domain
 @param bits the length of array v
 @param domain the domain of x (possible value of x is 0,1,2,...,domain-1)
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param v the result which is an mpz_t array
 */
void zn_to_z2(mpz_t x, int bits, mpz_t domain, int client_type, mpz_t v[]);
/**
 Convert two shared binary values to a shared value in a given domain. The true value can be computed by performing XORon the two shared binary values and by adding the two values the parties get. For example, value five in domain 7 may be shared in binary forms of
   - 110
   - 011

 and the two parties may get outputs like 
   - 6
   - 6
 or 
   - 1
   - 4

 As 1+4 = 6+6 = 5 mod 7. The two outputs are both correct
 
 @param v the length-bits binary array
 @param bits the length of array v
 @param domain the domain of x (possible value of x is 0,1,2,...,domain-1)
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result the shared result
 */
void z2_to_zn(mpz_t v[], int bits, mpz_t domain, int client_type, mpz_t result);

/**
 Access the i-th elements in array[].
 @param i a additive shared index value
 @param dim dimension of array
 @domain the domain of the return value
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param item the result (additive shared) 
*/
void s_access_array(mpz_t i, int dim, mpz_t array[], mpz_t domain, int client_type, mpz_t item);

#endif   //_SFIXNUM_H_
