#ifndef SMC
#define SMC
#include <stdint.h>
#define ulong128 unsigned __int128
#define long128 __int128
#include <gmp.h>
#include "two_parties/two_parties.h"
#include "two_parties/cSFloat.h"

/** 
  define ALICE=1, BOB=2
 */
enum { ALICE=1,BOB};

#endif

