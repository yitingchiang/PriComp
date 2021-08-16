#define _WITH_GETLINE
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <ctime>
#include <sys/resource.h>
#include <cerrno>
#include <climits>
#include <cinttypes>
#include <sys/types.h>
#include <cstdint>
#include <cstring>

#include "../smc.h"
#include "sp_init.h"
#include "../communicate.h"
#include "../util.h"
#include "tables.h"

/**
 The table to compute exp(x).
 Defined and initialized in table.cpp
 */
extern SFloat **exp_table;

/**
 The table to compute log(x)
 Defined and initialized in table.cpp
 */
extern mpz_t *logtable;

/**
  for three party scalar protocol which has a commodity server. This socket is to communicate with the commodity server.
 */
int SP_TCP_c_socket;

/**
  For the scalar protocol which has a commodity server. This socket is for the parties to communicate with the commodity server.
 */
int SP_p_socket;//sockets to send recv data from peer

/**
 The socket for Alice to play as a server to wait for Bob to connect.
 */
int SP_mysocket;

/**
 The index to access the R vector.
 This vector is for offline commodity server.
 That is, is used to generate Ra and Rb in the theoretically secure scalar-product protocol.
 */
unsigned int R_index;
/**
 The index to access the r vector.
 This vector is for offline commodity server.
 That is, is used to generate ra and rb in the theoretically secure scalar-product protocol.
 */
unsigned int r_index;

/**
 Which adder is used as the zn_to_z2 operation.
 */
int adderIDX=0;

#ifdef EVAL_TIME
/**
 Count how many times the scalar-product protocol is performed.
 For performance evaluation.
 */
long long int t_count = 0;
#endif

/**
 The number of pairs of (Ra,Rb) and (ra,ab) in the offline random table.
 This variable is set in SP_init.
 */
unsigned long MAX_RND_LEN=0;

/**
 The offline random bits for Ra (or Rb).
 */
long128* R=NULL;

/**
 The offline random bits for ra (or rb).
 */
long128* r=NULL;

/**
 Write a length-dim unsigned long vector using the socket.
 @param vec A vector with dim elements
 @param socket The socket to write the data
 @param dim The length of vec
 @return The return value is the same as the send(). That is, it return the number of bytes sent, or -1 if error occurs  
 */
int send_vector(ulong128* vec,int socket,int dim)
{
  return protocol_write(socket,(char*)vec,(sizeof(long128)*(dim)));
}

/**
 Read a length-dim unsigned long vector using the socket.
 @param vec A vector with dim elements
 @param socket The socket to read the data
 @param dim The length of vec
 @return The return value is the same as the recv(). That is, it return the number of bytes sent, or -1 if error occurs  
 */
int get_vector(ulong128* vec,int socket, int dim)
{
  char* tmps=(char*)vec;
  return protocol_read(socket,tmps,sizeof(long128)*dim);
}

/**
 Copy the value in px to a field of sarg structure.
 @param sarg The sarg structure to copy data to.
 @param client_type The type of this party:
  - ALICE (1)
  - BOB (2)
 @param px The vector of data to copy from
 @param dim The length of the vector
 */
void set_spdata_multi(SPData *sarg,int client_type,mpz_t* px,int dim)
{
  if(client_type==ALICE)//ALICE
  {
    for(int i=0;i<dim;i++)
      sarg->vec_pa[i]=mpz_to_long128(px[i]);
  }
  else // BOB
  {
    for(int i=0;i<dim;i++) 
      sarg->vec_pb[i]=mpz_to_long128(px[i]);
  }
}

/**
 Compute the addition of nV vectors with some random vectors get from R. That is, compute multiple V+R where V is the vector hold by the party and R is the random vector Ra or Rb. The length of the i-th vector is dim[i] (i=0,1,...,nV-1). The summation of the length of these vectors is dim_len
 @param result The result
 @param get_vec To store the random vectors get from R
 @param local_vec The vectors of the party
 @param dim The length of the vectors.
 @param dim_len The sum of the length of the vectors
 @param nV The number of vectors in this operation
 */
void local_multi_get_and_add(ulong128* result, ulong128* get_vec,ulong128* local_vec, int* dim,int dim_len,int nV)
{
  int v=0;//count how many bytes have been get

  for(int i=0;i<nV;i++)
  {
    int d=dim[i]-1;//index in R

    for(int j=0;j<=d;j++)
    {
      get_vec[v+j]=R[(R_index)%MAX_RND_LEN];
      R_index++;
    }
    v+=(d+1);

  }
  R_index%=MAX_RND_LEN;//set R_index after getting data in this round;

  //compute result[i] = get_vec[i]+local_vec[i]
  vec_add(result,get_vec,local_vec,dim_len);
}


/**
 Compute result=X+Y. Y is given by the party, and X is used to store the data getting from the another party. In the scalar product protocol, X is used to get Ra (for Alice) or Rb (for Bob) from the commodity server.You have to allocate X's memory before calling this function.
 @param result The result vector
 @param socket The socket used to communicate with the peer party.
 @param X A vector to store the data get from another party. 
 @param Y The local vector to do vector addition.
 @pram total_dim The length of X (Y)
 @param domain The domain this operation performs on
 @return 1 if success. Otherwise -1.
 */
int multi_get_and_add(ulong128* result, int socket,ulong128* X,ulong128* Y, int total_dim,mpz_t domain)
{
  char* tmps=(char*)X;
  int total_len=sizeof(long128)*(total_dim);
  unsigned int vec_idx=0;//the index of the vector to do inner product
  unsigned int offset=0;	
  unsigned int fraction=0;//for small fractions that only a part of a long128 byte(s) is get. 
  int unget_len=total_len;

  int counter=0;
  int v=0;//count how many bytes have been get

  while(unget_len>0)
  {
    //get data from the peer party
    int s=recv(socket,tmps+v,unget_len,0);
    if(s<0) { perror("Error when read from socket!"); return -1; }
    unget_len-=s;
    v+=s;
    s+=fraction;

    int len=s/sizeof(long128);
    fraction=s%sizeof(long128);

    vec_add(result+offset,X+offset,Y+offset,len);
    offset+=len;
    counter++;
  }
  return 1;
}

/**
 Compute result=X*Y. Y is given by the party, and X is used to store the data getting from the another party. In the scalar product protocol, X is used to get Ra (for Alice) or Rb (for Bob) from the commodity server.You have to allocate X's memory before calling this function.
 Note that Y may contains more than one vectors. That is, There are nSC vectors Y_1,...,Y_{nSC} that the length of Y_i is dim[i]. This function gets X_i from another party, computing inner product X[i]*Y[i], and then put the result in result[i].
 @param result An 128-bit unsigned int vector with length nSC.
 @param socket The socket used to get another vector from the commodity party.
 @param X A vector to store the data get from another party. 
 @param Y The local vector to do vector addition.
 @param dim The length of the vectors. That is, dim[i] is the length of the i-th vector.
 @pram total_dim The length of X (Y)
 @param domain The domain this operation performs on
 @return 1 if success. Otherwise -1.
*/
int multi_get_and_product(ulong128* result, int socket,ulong128* X,ulong128* Y, int* dim,int total_dim,int nSC,mpz_t domain)
{
  int v=0;//count how many bytes have been get
  char* tmps=(char*)X;
  int total_len=sizeof(long128)*(total_dim);
  unsigned int start_idx=0;
  unsigned int vec_idx=0;//the index of the vector to do inner product
  ulong128 local_result=0;
  unsigned int dim_offset=0;
  unsigned int fraction=0;//for small fractions that only a part of a long128 byte(s) is get. 

  int unget_len=total_len;

  int remained_len=dim[0];
  while(v<total_len)
  {
    int s=recv(socket,tmps+v,unget_len,0);
    if(s<0) { perror("Error when read from socket!"); return -1; }

    unget_len-=s;
    v+=s;
    s+=fraction;
    int len=s/sizeof(long128);
    fraction=s%sizeof(long128);

    while(len>0)
    {
      bool has_complete_vector=(remained_len-len<=0)||(v>=total_len) ? true : false;
      int local_len=has_complete_vector==true ? remained_len: len;
      remained_len-=len;

      ulong128 new_local_result=inner_product(X+dim_offset,Y+dim_offset,local_len,domain);

      local_result+=new_local_result;
      dim_offset+=local_len;
      len-=local_len;

      //get a complete vector;
      //if local_len>0, it means that thereveived data contains more than one (may be less than two) vectors
      //if local_len=0, it means exactly (or with a small fraction of one long128) one vector is get
      if(has_complete_vector)
      {
        result[vec_idx]=local_result; 
        local_result=0; 
        vec_idx++;
        if(vec_idx<nSC) 
        {
          remained_len=dim[vec_idx];
        }
      }
    }
  }
  return 1;
}

/**
 Compute n times of inner products on two vectors X and Y. Note that X and Y are two-dimentioanl matrix represented as one-dimentional vector. That is, there are n vectors in X and Y. X and Y are the flattened form of the vectors in them. The i-th vector in X (and Y) has length dim[i]. So the length of X (and) Y is dim[0]+dim[1]+...+dim[n-1]. After this function is done, the i-th result is in result[i].
 @param result The results which is a length-n vector
 @param X The vector X
 @param Y The vector Y
 @param dim A vector that the length of the vectors conducting the i-th inner product is dim[i]
 @param n How many vectors are there in X (and Y). And also the number of inner product computations this function will conduct.
 @param domain The domain these computations are performed on.
 */
void multi_inner_product(ulong128* result, ulong128* X,ulong128* Y, int* dim,int n,mpz_t domain)
{
  result[0]=inner_product(X,Y,dim[0],domain);
  
  unsigned int dim_offset=dim[0];

  for(int i=1;i<n;i++)
  {
    result[i]=inner_product(X+dim_offset,Y+dim_offset,dim[i],domain);
    dim_offset+=dim[i];
  }
}

/**
  An implementation scalar product protocol that performs scalar product on multiple pairs of vectors.
  In addition, this implementation can support performing the protocol using offline random bits. That is, reading random bits from local files instead of from a commodity server.
  @param sarg some parameters and variables used for performing this protocol. The structure is defined in sp/client.
  @param socket The socket to send/recv data from the peer party
  @param nSC The number of scalar-product protocols to perform. This is the length of dim.
  @param dim The length the vectors to perform the scalar-product protocol. dim[k] is the k-th vector's length. There are totally nSC vectors, so |dim|=nSC.
  @param total_dim Summation of dim[i] where i=0,1,2,...,nSC-1.
  @param domain The domain these scalar product protocols perform on.
  //@param tbBegin For time evaluation. This param will be removed later.
*/
int run_multi_bob(SPData* sarg,int socket,int nSC, int* dim, int total_dim, mpz_t domain)
{
  ulong128* Xa=sarg->vec_pa;//Alice's private vector Xa, or the Xa' that Bob gets
  ulong128* Xb=sarg->vec_pb;//Bob's private vector, or the Xb' that Alice gets
  ulong128* sp=sarg->sp;//temp vector
  ulong128* ra=sarg->ra;//random variable from commodity server
  ulong128* rb=sarg->rb;//random variable from commodity server
  ulong128* Ra=sarg->vec_Ra; //random vector from commodity server
  ulong128* Rb=sarg->vec_Rb; //random vector from commodity server
  ulong128* u=sarg->u;
  ulong128* v=sarg->v;

  int len=0; 

  //get Rb and compute Ra=Rb+Xb
  //c_socket>0 means using commodity server
  if(SP_TCP_c_socket>0)
  {
    multi_get_and_add(Ra,SP_TCP_c_socket,Rb,Xb,total_dim,domain);
  }
  else // use local random bits
  {
    local_multi_get_and_add(Ra,Rb,Xb,dim,total_dim,nSC);
  }

  //get Xa' from Alice and compute ra=Xa'*Xb.
  multi_get_and_product(ra,socket,Xa,Xb,dim,total_dim,nSC,domain);
  send_vector(Ra,socket,total_dim);

  if(SP_TCP_c_socket>0) 
	  get_vector(rb,SP_TCP_c_socket,nSC);
  else //use offline random bits
  {
    for(int i=0;i<nSC;i++)
    {
      rb[i]=0;
      int d=dim[i];
      for(int j=0;j<d;j++)
      {
        rb[i]+=r[(r_index)%MAX_RND_LEN];
        r_index++;
      }
    }
    r_index%=MAX_RND_LEN;//set r_index after getting random data in this round
  }

  for(int i=0;i<nSC;i++)
  {
    u[i]=MOD(random(),domain);

    //compute s
    v[i]=MOD((rb[i]-u[i]+ra[i]),domain);
  }

  send_vector(v,socket,nSC);
#ifdef EVAL_TIME
  t_count++;
#endif
  return 0;
}

/**
  An implementation scalar product protocol that performs scalar product on multiple pairs of vectors.
  In addition, this implementation can support performing the protocol using offline random bits. That is, reading random bits from local files instead of from a commodity server.
  @param sarg some parameters and variables used for performing this protocol. The structure is defined in sp/client.
  @param socket The socket to send/recv data from the peer party
  @param nSC The number of scalar-product protocols to perform. This is the length of dim.
  @param dim The length the vectors to perform the scalar-product protocol. dim[k] is the k-th vector's length. There are totally nSC vectors, so |dim|=nSC.
  @param total_dim Summation of dim[i] where i=0,1,2,...,nSC-1.
  @param domain The domain these scalar product protocols perform on.
  //@param tbBegin For time evaluation. This param will be removed later.
 */
int run_multi_alice(SPData* sarg,int socket,int nSC, int* dim,int total_dim,mpz_t domain)
{
  ulong128* Xa=sarg->vec_pa;//Alice's private vector Xa, or the Xa' that Bob gets
  ulong128* Xb=sarg->vec_pb;//Bob's private vector, or the Xb' that Alice gets
  ulong128* sp=sarg->sp;//temp vector
  ulong128* ra=sarg->ra;//random variable from commodity server
  ulong128* rb=sarg->rb;//random variable from commodity server
  ulong128* Ra=sarg->vec_Ra; //random vector from commodity server
  ulong128* Rb=sarg->vec_Rb; //random vector from commodity server
  ulong128* u=sarg->u;
  ulong128* v=sarg->v;

  int len=0;
#ifdef EVAL_TIME
  struct timeval tvEnd, tvDiff, tvBegin0,tvBegin1,tvBegin2, tvBegin3;
#endif

  //use commodity server
  if(SP_TCP_c_socket>0)
    multi_get_and_add(Rb,SP_TCP_c_socket,Ra,Xa,total_dim,domain);
  else //use local random bits
    local_multi_get_and_add(Rb,Ra,Xa,dim,total_dim,nSC);
  send_vector(Rb,socket,total_dim);

  multi_get_and_product(rb,socket,Xb,Ra,dim,total_dim,nSC,domain);

  get_vector(sp,socket,nSC);//get s 
  if(SP_TCP_c_socket>0) // use commodity server
    get_vector(ra,SP_TCP_c_socket,nSC);
  else //use offline random bits
  {
    for(int i=0;i<nSC;i++)
    {
      ra[i]=0;	
      int d=dim[i];
      for(int j=0;j<d;j++)
      {
        ra[i]+=r[(r_index)%MAX_RND_LEN];
        r_index++;
      }
    }
    r_index%=MAX_RND_LEN;
  }

  for(int i=0;i<nSC;i++)
  {
    v[i]=MOD(((sp[i])+(ra[i])-(rb[i])),domain);
  }

#ifdef EVAL_TIME
  t_count++;
#endif
  return 0;
}

/**
 Initialize a mpz_t matrix (an array which represents a two-dimential matrix) using a 128-bit matrix. M will be used to perform the scalar product which using a non-singular Matrix.
 @param M the mpz_t matrix to be initialized
 @param 1M The 128-bit int matrix used to set M
 @param len The length of M and lM
 */
void initM(mpz_t M[],long128 lM[],int len)
{
  init_vector(M,len);
  for(int i=0;i<len;i++)
    mpz_set_long128(M[i],lM[i]);
}

/**
 Compute M*A, where M is a square matrix (represented as a vector of length dim*dim) and A is a length-dim vector. The result is set to V, an length-dim vector.
 @param V The result vector
 @param A The vector A in M*A
 @param M The square matrix M in M*A
 @param dim the length of A and V
 @param domain The domain this computation is performed on.
 @return The vector V
 */
long128* ArrayMatrixProduct(long128* V,ulong128* A, const long128* M, int dim,mpz_t domain)
{
  bzero(V,sizeof(ulong128)*dim);

  mpz_t px[dim];
  mpz_t py[dim];
  init_vector(py,dim);
  init_vector(px,dim);
  mpz_t result;
  mpz_init(result);

  for(int i=0;i<dim;i++) 
  {
    mpz_set_long128(px[i],A[i]);

  }

  for(int i=0;i<dim;i++)
  {
    for(int j=0;j<dim;j++)
    {
      mpz_set_long128(py[j],M[i*dim+j]);
    }
    inner_product(result,px,py,dim);
    V[i]=mpz_to_long128(result);
  }


  mpz_clear(result);
  clear_vector(px,dim);
  clear_vector(py,dim);

  return V;
}

/**
 Initialize the SPData struct, including memory allocation. For the definition of struct SPData, check util.h.
 @param px The party's private input.
 @param sarg The SPData structure
 */
int init_sarg(mpz_t *px, SPData *sarg, int* dim, int nSC, int client_type)
{
  int total_dim=0;
  for(int i=0;i<nSC;i++) { 
    total_dim+=dim[i];
  }
  sarg->vec_pa=(ulong128*)malloc(sizeof(ulong128)*total_dim*4);
  if(!(sarg->vec_pa)) { fprintf(stderr,"Failed to allocate memory for running scalar product!\n"); return -1; }
  sarg->vec_pb=sarg->vec_pa+total_dim;
  sarg->vec_Ra=sarg->vec_pa+total_dim*2;
  sarg->vec_Rb=sarg->vec_pa+total_dim*3;

  sarg->sp=(ulong128*)malloc(sizeof(ulong128)*nSC*5);

  if(!(sarg->sp)) { fprintf(stderr,"Failed to allocate memory for running scalar product!\n"); free(sarg->vec_pa); return -1; }
  //set memory address
  sarg->ra=sarg->sp+nSC;
  sarg->rb=sarg->sp+nSC*2;
  sarg->u=sarg->sp+nSC*3;
  sarg->v=sarg->sp+nSC*4;
  set_spdata_multi(sarg,client_type,px,total_dim);
  return total_dim;
}

/**
 Free the memory used by a SPData struct. The memory is allocated in function init_sarg. For the definition of struct SPData, check util.h.
 @param sarg The SPData struct.
 */
void free_sarg(SPData *sarg)
{
  free(sarg->vec_pa);
  free(sarg->sp);
}

/**
 The implementation of secure scalar-product protocol which leaks some amount of information. This is for Bob. This protocol only supports vectors with dimension 2,3,11,53,55, or 108. In addition, if the dimension is 1, 256, or 32, the three-party (one commodity server) scalar-product protocol that using offline random bits (so the commodity server is not necessary) will be used.
 Check the following papers:
 Du, W. & Zhan, Z., A Practical Approach to Solve Secure Multi-party Computation Problems, Proceedings of New Security Paradigms Workshop, 2002.
 Chiang, Y.-T.; Wang, D.-W.; Liau, C.-J. & Hsu, T.-s., Secrecy of Two-Party Secure Computation, Data and Applications Security XIX, Springer Berlin Heidelberg, 2005, 114-123.
 @param sarg The data structure to perform scalar-product protocol. Check util.h to see the definition of struct SPData.
 @param socket The socet to connect with the peer party
 @param dim Dimension of the input vector
 @param domain The domain this operation is performed on
 
*/
int run_bob_matrix(SPData *sarg,int socket,int dim, mpz_t domain)
{
  if(dim==1||dim==256||dim==32)
  {
    return  run_multi_bob(sarg,socket,1,&dim,dim,domain);
  }

  const long128* M[6]={invM55,invM11,invM2,invM3,invM53,invM108};
  int midx=(dim==53) ? 4: dim%5;
  midx=(dim==108) ? 5 : midx;

  int send_len=dim/2;
  int get_len=dim-send_len;

  ulong128* Xa=sarg->vec_pb;//Alice's private vector Xa
  long128 Ra[dim];
  //=sarg->vec_Ra;//to send
  long128 Rb[dim];
  //=sarg->vec_Rb;//to get
  ulong128* u=sarg->u;

  ArrayMatrixProduct(Ra,Xa,M[midx],dim,domain);

  long128* ptr=Rb;
  if(protocol_read(socket,(char*)ptr,sizeof(long128)*get_len)<0)
  {
    fprintf(stderr,"Failed to get vector from Alice!\n");
  }
  ptr=Ra+get_len;
  if(protocol_write(socket,(char*)ptr,sizeof(long128)*send_len)<0)
  {
    fprintf(stderr,"Failed to send vector from Alice!\n");
  }

  long128 local_u=0;

  for(int i=0;i<get_len;i++)
  {
    local_u+=Ra[i]*Rb[i];
  }
  *u=MOD(local_u,domain);
  free_sarg(sarg);
#ifdef EVAL_TIME
  t_count++;
#endif
  return 0;
}

/**
 The implementation of secure scalar-product protocol which leaks some amount of information. This is for Alice. This protocol only supports vectors with dimension 2,3,11,53,55, or 108. In addition, if the dimension is 1, 256, or 32, the three-party (one commodity server) scalar-product protocol that using offline random bits (so the commodity server is not necessary) will be used.
 Check the following papers:
 Du, W. & Zhan, Z., A Practical Approach to Solve Secure Multi-party Computation Problems, Proceedings of New Security Paradigms Workshop, 2002.
 Chiang, Y.-T.; Wang, D.-W.; Liau, C.-J. & Hsu, T.-s., Secrecy of Two-Party Secure Computation, Data and Applications Security XIX, Springer Berlin Heidelberg, 2005, 114-123.
 @param sarg The data structure to perform scalar-product protocol. Check util.h to see the definition of struct SPData.
 @param socket The socet to connect with the peer party
 @param dim Dimension of the input vector
 @param domain The domain this operation is performed on
 
*/
int run_alice_matrix(SPData* sarg,int socket,int dim,mpz_t domain)
{
  //call 
  if(dim==1||dim==256||dim==32)
  {
    return run_multi_alice(sarg,socket,1,&dim,dim,domain);
  }

  const long128* M[6]={M55,M11,M2,M3,M53,M108};
  int midx=(dim==53) ? 4: dim%5; 
  midx=(dim==108) ? 5 : midx;

  int get_len=dim/2;
  int send_len=dim-get_len;
  ulong128* Xa=sarg->vec_pa;//Bob's private vector Xb
  long128 Ra[dim];
  long128 Rb[dim];
  ulong128* v=sarg->v;
  ArrayMatrixProduct(Ra,Xa,M[midx],dim,domain);

  long128* ptr=Ra;

  if(protocol_write(socket,(char*)ptr,sizeof(long128)*send_len)<0)
  {
    fprintf(stderr,"Failed to send vector to Bob!\n");
  }

  ptr=Rb;
  if(protocol_read(socket,(char*)ptr,sizeof(long128)*get_len)<0)
  {
    fprintf(stderr,"Failed to get vector from Bob!\n");
  }
  long128 local_v=0;

  for(int i=0;i<get_len;i++)
  {
    local_v+=Ra[send_len+i]*Rb[i];
  }
  *v=MOD(local_v,domain);
  free_sarg(sarg);
#ifdef EVAL_TIME  
  t_count++;
#endif
  return 0;
}

/** 
 The procedure for secure product protocol. This procedure conduct secure inner product on two vector.
 Each party has to provide their secure vector.
 @param px Private input of the party with datatype mpz_t, which is a (pointer to) integer.
 @param client_type The type of this party:
  - ALICE (1)
  - BOB (2)
 @param dim Dimension (Length) of the secure vector px
 @param dom The domain
 @param result The (additively shared) result
 @return return value: 
  - 0: successfully performed the computation. The result is stored in result
  - negative value: Failed to perform the computation.
 */
int scalar_product(mpz_t *px, int client_type, int dim, mpz_t dom,mpz_t result)
{
  //int status=0;//result status of this procedure
  int r=0;//status of run_bob or run_alice
  mpz_t single_r[1];
  mpz_init(single_r[0]);

  int m_return=multi_scalar_product(px,client_type,1,&dim,dom,single_r);
  mpz_set(result,single_r[0]);
  mpz_clear(single_r[0]);

  return m_return;

}

/**
  The main procedure to run multiple scalar products.
  This function can conduct multiple scalar products at once.
  @param px The private input of alice or bob.
  @param client_type The type of this party:
  - ALICE (1)
  - BOB (2)
  @param dim Dimensionality of px
  @param domain Domain of px and this operation
  @param result Result of these scalar products
  @param nSC The number of vectors to conduct scalar products
  @return An integer to indicat if this function terminated normally.
  - 0 successfully exectued the computation. The result is stored in u
  - negative value: Failed to execute the computation.
 */
int multi_scalar_product(mpz_t *px, int client_type, int nSC, int* dim, mpz_t dom,mpz_t* result)
{
  int status=0;//result status of this procedure
  int r=0;//status of run_bob or run_alice

  struct timeval tvBegin,tvEnd, tvDiff;

  SPData sarg;
  int total_dim=init_sarg(px,&sarg,dim,nSC,client_type);
  if(total_dim<0) return -1;

  //set socket id, and put px to smag->vec_pa (ALICE) or vec_pb (BOB)

  //use commodity party
  if(SP_TCP_c_socket>0)
  {
    if(client_type==ALICE)
    {
      protocol_writeint(SP_TCP_c_socket,nSC);

      protocol_write(SP_TCP_c_socket,(char*)dim,sizeof(int)*nSC);
    }
  }

  if(client_type==1)
  {
    r=run_multi_alice(&sarg,SP_p_socket,nSC,dim,total_dim,dom);
    if(r<0)//If there is something wrong 
    {
      fprintf(stderr,"ERROR when running run_multi_alice()\n");
      status=-3;
    }
    else 
    {
      char buf[256];

      for(int i=0;i<nSC;i++)
      {
        char* ptr=long128_to_str(buf,sarg.v[i]);
        mpz_set_str(result[i],ptr,10);
      }
    }
  }
  else
  {
    r=run_multi_bob(&sarg,SP_p_socket,nSC,dim,total_dim,dom);
    if(r<0)//Something wromg
    {
      fprintf(stderr,"ERROR when running run_multi_bob()\n");
      fflush(stderr);
      status=-3;
    }
    else 
    {
      char buf[256];
      for(int i=0;i<nSC;i++)
      {
        char* ptr=long128_to_str(buf,sarg.u[i]);
        mpz_set_str(result[i],ptr,10);
      }
    }
  }

  free_sarg(&sarg);

  return status;
}

/** 
  Init sockets for scalar product protocol
  The scalar product protocol will fail if SP_init returns a negative value
  @param file The setting file.
  @param client_type The type of this party:
  - ALICE (1)
  - BOB (2)
**/
int SP_init(char* file,int& client_type)
{

  Protocol_ARG arg;
  char* peer_host;//peer party's host
  int 	peer_port;//peer party's port
  char* c_server=NULL;//commodity server's host
  int 	c_port=0; //commodity server's port
  int 	myport; 
  int 	i=0;
  int 	ADD_PORTNO=10;

  if(parse_ProtocolARG(&arg,file)<0) return -1;
  adderIDX=arg.adderIDX;
  printf("init parameters from %s\n",file);
  init_net_params(client_type,c_server,c_port,peer_host,peer_port,myport,arg);

  printf("successfully init parameters from file %s\n",file);

  if(client_type==ALICE) 
  {
    SP_mysocket=setServer(myport);
    printf("Alice is waiting for connection on port %d\n",myport);
    SP_p_socket=wait_connection(SP_mysocket);
    if(SP_p_socket<0)
    { fprintf(stderr,"Failed to set socket"); perror(""); return -2;}

  }
  else//bob
  {
    printf("Bob start to connect to %s:%d...\n",peer_host,peer_port);
    SP_p_socket=setClient(peer_host,peer_port);
    if(SP_p_socket<0) { fprintf(stderr,"Can not connect to alice!\n"); return -2;}
  }

  //Online commodity server mode
  if(arg.Offline_Commodity==false)
  {
    if(client_type==BOB)//sync. wait for ithe peer to create connection with the commodity server first.		
    {
      int v=0;
      protocol_readint(SP_p_socket,&v);//sync. wait Alice to connect to the commodity server.

      //After Alice connected to the commodity, Bob starts to connect
      SP_TCP_c_socket=setClient(c_server,c_port);

      if(SP_TCP_c_socket<0) { fprintf(stderr,"Party %d failed to set connection to the commodity server\n",client_type); return -2;}

    }
    else//Alice's case
    {
      //Alice conncts to the commodity server first
      SP_TCP_c_socket=setClient(c_server,c_port);

      if(SP_TCP_c_socket<0) { fprintf(stderr,"Party %d failed to set connection to the commodity server\n",client_type); return -2;}
    protocol_writeint(SP_p_socket,1);//sync with Bob. Telling Bob that Alice's connection is done.

    }

  }
  else //using offline random bits
  {
    FILE* rp=NULL;
    int rnd_len=0;
    SP_TCP_c_socket=-1;//indicates that we don;t use commodity server.
    char rndfile[256];
    //open the file to read offline random bits
    strcpy(rndfile,arg.offlineRNDDIR);
    if(client_type==ALICE)
      strcat(rndfile,"/alice.dat");
    else 
      strcat(rndfile,"/bob.dat");
    rp=fopen(rndfile,"rb");
    if(!rp)
    {
      char tmpbuf[256];
      fprintf(stderr,"ERROR! Failed to init local random bits from file %s.\n", rndfile);
      fprintf(stderr,"current working dir: %s\n",getcwd(tmpbuf,256));
      return 1;
    }
    fread(&rnd_len,sizeof(int),1,rp);
    MAX_RND_LEN=rnd_len;
    R=(long128*)malloc(sizeof(long128)*MAX_RND_LEN);
    r=(long128*)malloc(sizeof(long128)*MAX_RND_LEN);
    size_t size=sizeof(long128)*rnd_len;
    fread(R,size,1,rp);
    fread(r,size,1,rp);

    R_index=0;
    r_index=0;
    fclose(rp);
  }

  init_ExpTable(client_type);
  init_logTable(client_type);

  return 0;
}

/**
 Close sockets
  @param client_type The type of this party:
  - ALICE (1)
  - BOB (2)
 */
void SP_clear(int client_type)
{
#ifdef EVAL_TIME	
  //number of scalar-protocol calls. For performance checking
  printf("Conducted scalar-product protocol for %lld times\n", t_count);
#endif
  //a party should call for end
  if(SP_TCP_c_socket>0)
  {
    if(client_type==ALICE)
    {
      //tell the commodity serevr to end this session
      protocol_writeint(SP_TCP_c_socket,0);
    }
    close(SP_TCP_c_socket);
  }
  else //offline commodity server mode
  {
    free(R);
    free(r);
  }
  close(SP_p_socket);

  if(client_type==ALICE) 
  {
    close(SP_mysocket);
  }
  free_ExpTable();
  free_logTable();
}

