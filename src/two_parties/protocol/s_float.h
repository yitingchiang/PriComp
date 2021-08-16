/**
 @file
 This header file defines most of the operations on secure floating-point variables (SFloat). Note that functions with the parameter client_type must be performed by both parties. In this case, the SFloat parameter are shared; Functions without the client_type parameter can be performed by one party. In this case, the SFloat is not shared.
*/

#ifndef S_FLOAT_H
#define S_FLOAT_H

#include <gmp.h>
#include "../cSFloat.h"
#include "s_fixnum.h"
#include "../../util.h"
//#include "../../debug.h"
#include "../tables.h"

/**
 * Convert a double to a SFloat variable
 @param f the variable to be converted
 @return f in SFloat format
 */
SFloat float2SFloat(double f);
/**
 Convert a SFloat variable to a double variable. Note that SFloat is usually shared by two parties. Therefore, this function is usually used for debugging only.
 @param f The variable to be converted
 @return The SFloat value in double data type.
 */
double sFloat2float(SFloat f);
/**
 Perfoem secure addition on two SFloat variables.
 @param f1 a shared secure floating-point variable
 @param f2 a shared secure floating-point variable
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @result a shared result of f1+f2
 */
SFloat f_plus(SFloat f1, SFloat f2, int client_type);
/**
 Given a additively shared variable v, the two parties get v1 and v2 that are component-wise shared of v. So to share a value k, make the sum of the two parties' input be k. For example, if k=10, then the two parties' input can be (0,10),(1,9),...,(9,1),(10,0). You can also use negative values (for example, (-1,11)).
 @param v the value to be shared.
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @result a shared result of v, which is the sum of the two parties' input.
 */
SFloat fpShare(double v, int client_type);
/**
 Perfoem secure subtraction on two SFloat variables.
 @param f1 a shared secure floating-point variable
 @param f2 a shared secure floating-point variable
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @result a shared result of f1-f2
 */
SFloat f_minus(SFloat f1, SFloat f2, int client_type);
/**
 Perfoem secure multiplication on two SFloat variables.
 @param f1 a shared secure floating-point variable
 @param f2 a shared secure floating-point variable
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @result a shared result of f1*f2
 */
SFloat f_product(SFloat f1, SFloat f2, int client_type);
/**
 Perfoem secure division on two SFloat variables.
 @param f1 a shared secure floating-point variable
 @param f2 a shared secure floating-point variable
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @result a shared result of f1/f2
 */
SFloat f_division(SFloat f1, SFloat f2, int client_type);
/**
 Compute relu(v).
 @param v a shared secure floating-point variable
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @result a shared result of relu(v)
 */
SFloat f_relu(SFloat v, int client_type);
/**
 Test if v is positive or negative.
 @param v a shared secure floating-point variable
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @result a shared result
    - 1 if v > 0
    - 0 otherwise
 */
SFloat f_is_positive(SFloat v,int client_type);
/**
 Test if v is positive or negative.
 @param v a shared secure floating-point variable
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result The shared result in gmp mpz_t format
    - 1 if v > 0
    - 0 otherwise
 */
void f_is_positive(SFloat f, int client_type, mpz_t result);

/**
 Test if f1 is greater than f2.
 @param f1 a shared secure floating-point variable
 @param f2 a shared secure floating-point variable
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result The shared result in gmp mpz_t format
    - 1 if f1 > f2
    - 0 otherwise
*/
void   f_great_than(SFloat f1, SFloat f2, int client_type,  mpz_t result);

/**
 Test if f is zero.
 @param f a shared secure floating-point variable
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result The shared result in gmp mpz_t format
    - 1 if f is zero
    - 0 otherwise
*/
void f_is_zero(SFloat f,int client_type,mpz_t result);
/**
 Test if f is a negative value.
 @param f a shared secure floating-point variable
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result The shared result in gmp mpz_t format
    - 1 if f is negative
    - 0 otherwise
*/
void f_is_negative(SFloat f,int client_type,mpz_t result);
/**
 Test if f1 < f2.
 @param f a shared secure floating-point variable
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result The shared result in gmp mpz_t format
    - 1 if f1 < f2
    - 0 otherwise
*/
void f_less_than(SFloat f1, SFloat f2, int client_type, mpz_t result);

/**
 Check if two SFloat values are equal.
 @param f1 The first SFloat value
 @param f1 The second SFloat value
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @param result The shared result in gmp mpz_t format
    - 1 if f1 is equal to f2
    - 0 otherwise
*/
void f_equal(SFloat f1,SFloat f2,int client_type,mpz_t result);

/**
 Return a result according to the flag.
 @param flag a mpz_t integer which should be false (0) or true (1)
 @param f1 the first candidate value to be returned
 @param f2 the second candidate value to be returned
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @return f1 if flag is true; f2 if flag is false.
 */
SFloat f_if_then_else(mpz_t flag, SFloat f1, SFloat f2, int client_type);
/**
 Compute e^x, where e is the Euler number (Napier's contant) 2.71828...
 @param x the power
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @return e^x 
 */
SFloat f_exp(SFloat x, int client_type);

/**
 An approximation of ln(x).

 ln(y)=ln(2^n*x)=n*ln(2)+ln(x), for x in [1.0,2.0), ln(x) can be approximated as:
    - ln(x)= -1.7417939 + (2.8212026 + (-1.4699568 + (0.44717955 - 0.056570851 * x) * x) * x) * x. 

 (From https://stackoverflow.com/questions/9799041/efficient-implementation-of-natural-logarithm-ln-and-exponentiation).

 This function use another approximation.
    - ln(x) = -1.9410661+(3.5293114+(-2.4612314+(1.1306327+(-.28874220+0.31104322e-1*x)*x)*x)*x)*x

 @param x an SFloat
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @return ln(x) 
*/
SFloat f_ln(SFloat x, int client_type);

#endif
