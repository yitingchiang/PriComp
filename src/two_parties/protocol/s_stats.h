
#include "../cSFloat.h"
#ifndef S_STATS_H
#define S_STATS_H
/**
 Compute the summation of a SFloat array.
 @param operands The SFloat array
 @param dim The length of the input array
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @return Summation of the elements in the input array   
 */
SFloat f_summation(SFloat *operands, int dim, int client_type);
/**
 Compute the average of the values in a SFloat array.
 @param operands The SFloat array
 @param dim The length of the input array
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @return Average of the elements in the input array   
 */
SFloat f_average(SFloat *operands, int dim, int client_type);
/**
 Compute the average of the values in two arrays hold by the two parties.
 The scenario of using this function is that one array is separated in to two arrays, hold by the two parties. This function can compute the average of the elements in the two arrays such that the two parties cannot know the array hold by another party.
 @param party_dim The length of array hold by another party.
 @param local_data The array hold by this party.
 @param local_dim The length of the array hold by this party
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @return The average of the elements in the two arrays
 */
SFloat f_average_2(int party_dim, double* local_data, int local_dim, int client_type);
/**
 Compute the variance of the values in a SFloat array.
 @param operands The SFloat array
 @param dim The length of the input array
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @return The variance of the elements in the input array   
 */
SFloat f_variance(SFloat *operands, int dim, int client_type);
/** 
  Compute the median of the values in two SFloat arrays.
  Party A and Party B holds the two SFloat array. NOTICE: This function requires that data1 and data2 are respectively sorted in ascendant order.
  @param data1: Party A's data. Party B has to give an Zero array of langth dim1 
  @param data2: Party B's data. Party A has to give an Zero array of langth dim2
  @param dim1: The length of array in Party A
  @param dim2: The length of array in Party B
  @return The median of all the elements in the two arrays.
*/
SFloat f_median(SFloat *data1, SFloat* data2, int dim1, int dim2,int client_type);

/**
 Compute the variance of the values in two arrays hold by the two parties.
 The scenario of using this function is that one array is separated in to two arrays, hold by the two parties. This function can compute the average of the elements in the two arrays such that the two parties cannot know the array hold by another party.
 @param party_dim The length of array hold by another party.
 @param local_data The array hold by this party.
 @param local_dim The length of the array hold by this party
 @param client_type The type of this party:
    - ALICE (1)
    - BOB (2)
 @return The variance of the elements in the two arrays
 */
SFloat f_variance_2(int party_dim, double* local_data, int local_dim, int client_type);
#endif
