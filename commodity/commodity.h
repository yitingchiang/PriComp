#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#include <unistd.h>
#include <string.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <netinet/in.h>
#include "smc.h"
#define ELEMENT_SIZE sizeof(long128)
/**
  For sending/receiving random bits to/from the parties/commodity server. The random bits are represented as 128-bit integers.
 */
typedef struct
{
    ulong128* msg_to[2];/**< For Ra or Rb */
    ulong128* r[2]; /**< For ra or rb */
} MSG;
/**
 The random bits pool.
 */
typedef struct
{
    ulong128* data[4];
} POOL;


#define Q (1<<16)
#define AIDX 0//alice's index
#define BIDX 1//

#define READIDX 0
#define WRITEIDX 1

