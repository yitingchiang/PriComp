/**
 cSFloat.h defines the data structure of SFloat, which is a component-wisedly shared floating-point data type. Specifically, the sign bit, exponent, and mantissa are additively shared between two parties. Some methods are also defined.
 */
#ifndef SFLOAT_H
#define SFLOAT_H

#include <gmp.h>

/**
 Sign bit length
 */
#define SIGN_FIELD 1
/**
 Number of bits of the exponent in SFloat
 */
#define EXPONENT_FIELD 11
/**
 Number of bits of the mantissa in SFloat
 */
#define MANTISSA_FIELD 55
/**
 The domain of operations on the sign bit.
 2^SIGN_FIELD
 */
#define DOMAIN_S 2
/**
 The domain of operations on the exponent.
 2^EXPONENT_FIELD
 */
#define DOMAIN_E 2048
/**
 The domain of operations on the mantissa.
  2^MANTISSA_FIELD
 */
#define DOMAIN_M 36028797018963968   

/**
 * A secure floating-point class that each components are additive sharing.
 * The length (bits) of exponent and mantissa are defined in s_float.h
 */
class SFloat
{
  public:
    mpz_t infinity; /**< If the value is infinite */ 
    mpz_t NaN;      /**< If the value is NaN */
    mpz_t s;        /**< sign bit */
    mpz_t e;        /**< exponent */
    mpz_t m;        /**< mantissa */
    /**
     * Default constructor. Set all elements to 0.
     */
    SFloat(); 

    /**
     The constructor of SFloat. Copy the value from another SFloat variable.
     @param s The SFloat variable to init a new SFloat variable.
     */
    SFloat(const SFloat &s);

    /**
      Initialize the SFloat object.
      @param iinfinity
        - 1 if the value is infinity
        - 0 (default) otherwise
      @param inan
        - 1 if is NaN
        - 0 (default) not NaN 
     */	
    SFloat(int iinfinity, int inan);
    /**
     * Override operator =
       @param v The SFloat value to set from.
     */
    SFloat& operator=(const SFloat& v);
    /**
       * Set this value to the Maximum vale
     */
    void setMax();
    /**
     * Set this value by giving the sign bit, exponent, and mantissa.
     * @param is sign bit
     * @param ie exponent
     * @param im mantissa
     */
    void set(mpz_t is, mpz_t ie, mpz_t im);

    /**
     * Set this value by the given sign bit, exponent, and mantissa which are all non-negative values.
     * @param is an unsigned long value as the sign bit
     * @param ie an unsigned long value as the exponent
     * @param im an unsigned long value as the mantissa
     */
    void set_ui(unsigned long int is, unsigned long int ie, unsigned long int im);
    
    /**
     * To negate this value. That is, make this value s to be -s
       Notice that a SFloat variable is component-wised sharing.
       So only one party needs to call this function.
     */
    void negate();
    /**
     * Write out a floating-pount value.
     * @param fp A c-style file pointer to write the data
     */
    void write(FILE* fp);
    /**
     * Read a floating-pount value and set.
     * @param fp A c-style file pointer to read data from
     */
    void read(FILE* fp);
    /** 
     * Write elements using gmp functions.
     * @param fp A c-style file pointer to write the data to.
     */
    void writeRaw(FILE* fp);
    /**
     * Read a SFloat using gmp functions.
     * @param fp A c-style file pointer to read data from.
     */
    void readRaw(FILE* fp);
    /**
       Print the sign bit, exponent, and mantissa to stderr.
       This function is mostly for debugging.
     */
    void to_s();
};

#endif
