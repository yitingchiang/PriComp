#include "time_eval.h"
#include <ctime>
#include <cstdio>
#include <cstring>

void time_set(struct timeval* t,int i)
{
  gettimeofday(t+i, NULL);
}
/*
void time_print(struct timeval *Stamp, int startIDX, int endIDX)
{
  struct timeval tvDiff;
  timeval_sub(&tvDiff, Stamp+endIDX, Stamp+startIDX);
  fprintf(stderr,"Time between timestamp %d and %d :\n%ld.%06ld\n",startIDX,endIDX,tvDiff.tv_sec, tvDiff.tv_usec);
}
*/
void time_print(struct timeval* t, int startIDX,int endIDX)
{
  struct timeval tvDiff;
  char buf_start[32];
  char buf_end[32];
  char* tmps=NULL;

  timeval_sub(&tvDiff, t+endIDX, t+startIDX);
  //convert time_t to string format
  ctime_r(&(t[endIDX].tv_sec),buf_end);
  ctime_r(&(t[startIDX].tv_sec),buf_start);
  //cut the year part
  tmps=strrchr(buf_start,' ');
  *tmps='\0';
  tmps=strrchr(buf_end,' ');
  *tmps='\0';

  fprintf(stdout,"Timestamp %d: %s.%ld\n",startIDX,buf_start,t[startIDX].tv_usec);
  fprintf(stdout,"Timestamp %d: %s.%ld\n",endIDX,buf_end,t[endIDX].tv_usec);
  fprintf(stdout,"Time between timestamp %d and %d :\n%ld.%06ld\n", startIDX,endIDX,tvDiff.tv_sec, tvDiff.tv_usec);
}

/* time measure related (milisec) */
int timeval_sub(struct timeval *result, struct timeval *t2, struct timeval *t1)
{
  long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec);
  result->tv_sec = diff / 1000000;
  result->tv_usec = diff % 1000000;

  return (diff<0);
}

int timeval_add(struct timeval *result, struct timeval *t2)
{
  long int diff = (t2->tv_usec + 1000000 * t2->tv_sec);
  diff += result->tv_usec+1000000*result->tv_sec;

  result->tv_sec = diff / 1000000;
  result->tv_usec = diff % 1000000;
  return (diff<0);
}

void timeval_print(struct timeval *tv)
{
  char buffer[30];
  time_t curtime;

  printf("%ld.%06ld", tv->tv_sec, tv->tv_usec);
  curtime = tv->tv_sec;
  strftime(buffer, 30, "%m-%d-%Y  %T", localtime(&curtime));
  printf(" = %s.%06ld\n", buffer, tv->tv_usec);
}

