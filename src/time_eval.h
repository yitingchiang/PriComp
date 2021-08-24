/**
 This file provides function to do time performance evaluation.
 */
#ifndef TIME_EVAL
#define TIME_EVAL
#include <sys/time.h>

/**
 Call gettimeofday to set a timestamp in an array
 @param timeval An timestamp array
 @@param t Which timestamp to set
 */
void time_set(struct timeval*, int t);

/**
 Print the time of timestamp t[startIDX] and t[endIDX] and compute the time between these two timestamps.
 @param t The timestamp array
 @param startIDX startIDX
 @param endIDX endIDX
 */
void time_print(struct timeval* t, int startIDX,int endIDX);

/**
 Compute the time between t1 and t2 (t1 < t2). The result is put in the first parameter. The time is measure in milliseconds 
 @param result The result.
 @param t2 The end time
 @param t1 The start time
 */
int timeval_sub(struct timeval *result, struct timeval *t2, struct timeval *t1);
/**
 Add two time stamps and save the result to the first one.
 @param result The first time stamp and tje result.
 @param t2 The second time stamp.
 */
int timeval_add(struct timeval *result, struct timeval *t2);
/**
 Print a time stamp in format MM-DD-YY H:M:S where H in [0,23], M in [0,59], and S in [0,60].
 @param tv The time stamp to print.
 */
void timeval_print(struct timeval *tv);

/**
 Print the time between two timestamps Stamp[startIDX] and Stamp[endIDX].
 @param *Stamp The time stamp array.
 @param startIDX The start point.
 @param endIDX The end ppoint.
 */
void time_print(struct timeval *Stamp, int startIDX, int endIDX);
#endif
