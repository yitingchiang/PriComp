/*
 A simple commodity-server implementation to generate random bots online.
 * */
#include <stdbool.h>
#include "commodity.h"
#include "util.h"
#include "communicate.h"

//#define DEBUG

void usage(const char *fname)
{
  fprintf(stderr,"The commodity server\n");
  fprintf(stderr,"usage: %s <conf_file>\n", fname);
  fprintf(stderr,"\t<conf_file> is the configuration file of Pricomp. Check the cfg files in the conf directory.");
  fprintf(stderr,"Notice that the commodity server have to be run BEFORE Alice and Bob.\n");
}

/**
 Free allocated memory for messages
 @param m The memory address of the message to free
 */
void free_msg(MSG* m)
{
  free(m->msg_to[0]);
  free(m->msg_to[1]);
  for(int i=0;i<2;i++)
    if(m->r[i]) free(m->r[i]);
  free(m);
}

/**
 Free allocated memory for the message pool
 @param m The memory address of the message pool
 */
void free_pool(POOL* p)
{
  int i=0;
  for(i=0;i<4;i++) free(p->data[i]);
  free(p);
}
/** 
 * Init Commodity Server. Read all data from the files.
 */
POOL* init_commodity_server(const char* file)
{
  POOL* p=NULL;
  FILE* rp=NULL;
  char* buf=NULL;
  char* ptr1=NULL;
  char* ptr2=NULL;
  int dim=0;
  size_t size=0;
  int i=0;
  int j=0;

  rp=fopen(file,"r");
  if(!rp)
  {
    fprintf(stderr,"Failed to open file %s to prepare random bits.\n",file);
    return NULL;
  }

  p=(POOL*)malloc(sizeof(POOL));
  if(!p)
  {
    fprintf(stderr,"Failed to generate MSG pool!\n");
    fclose(rp);
    return NULL;
  }

  for(i=0;i<4;i++)
  {
    size=getline(&buf,&size,rp);

    dim=size/21;
    p->data[i]=(ulong128*)malloc(dim*sizeof(long128));//"Ra,ra\0". 

    if(!p->data[i])
    {
      fprintf(stderr,"Failed to allocate %ld memory for %d-th pool data!\n",size,i);
      for(j=0;j<i;j++) free(p->data[j]);
      free(buf);
      fclose(rp);
      free(p);
      return NULL;
    }

    ptr1=buf;

    char tmps[256];

    for(j=0;ptr1-buf<size;j++)
    {
      p->data[i][j]=strtol(ptr1,&ptr2,10);
      ptr1=ptr2;
      ptr1++;
    }
  }

  free(buf);
  fclose(rp);

  return p;
}

//nSC: how many scalar products is here to perform
//the length of dim array=nSC
MSG* set_msg(POOL* p,int* dim,int nSC)
{
  MSG* m=NULL;
  char* ptr1=NULL;
  char* ptr2=NULL;
  int i=0;
  int j=0;
  ulong128 q[2];
  ulong128 rnd[2];

  m=(MSG*)malloc(sizeof(MSG));

  if(!m)
  {
    fprintf(stderr,"Failed to generate MSG pool!\n");
    return NULL;
  }

  int d_sum=0;

  for(i=0;i<nSC;i++) d_sum+=dim[i];

  //the total bytes to send is the sum of length of the vectors times sizeof(ulong128), which is for Ra (length=d_sum), and nSC long128 integers, which is for ra
  size_t size=sizeof(ulong128)*(d_sum+nSC);
  //allocate memory for msg
  for(i=0;i<2;i++)
  {
    m->msg_to[i]=(ulong128*)malloc(size);//"Ra,ra\0". to alice
    if(!m->msg_to[i])
    {
      fprintf(stderr,"Failed to allocate %ld memory for msg.\n",size);
      for(j=0;j<i;j++) free(m->msg_to[j]);
      return NULL;
    }
  }

  //preparing msgs. there are nSC sequences to be set.
  int offset=0;

  for(int k=0;k<nSC;k++)
  {
    //set random data for the k-th vector
    int local_dim=dim[k];//the length of the k-th vector
    //copy Ra and Rb to msg
    for(i=0;i<2;i++)
    {
      for(j=0;j<local_dim;j++)
      {
        m->msg_to[i][offset+j]=p->data[i][j];
      }
    }

    for(i=2;i<4;i++)
    {
      q[i-2]=p->data[i][local_dim-1];
      if(q[i-2] == 0) rnd[i-2] = 0;
      else
        rnd[i-2]=random()%q[i-2];//smaller than q[i]
    }

    m->msg_to[0][offset+local_dim]=(q[0]-rnd[0])*(Q)+(q[1]-rnd[1]);//ra
    m->msg_to[1][offset+local_dim]=rnd[0]*(Q)+rnd[1];//rb
    offset+=(local_dim+1);
  }
  return m;
}


//nSC: how many scalar products is here to perform
//the length of dim array=nSC
MSG* set_multi_msg(POOL* p,int* dim,int nSC)
{
  MSG* m=NULL;
  char* ptr1=NULL;
  char* ptr2=NULL;
  int i=0;
  int j=0;
  ulong128 q[2];
  ulong128 rnd[2];

  m=(MSG*)malloc(sizeof(MSG));

  if(!m)
  {
    fprintf(stderr,"Failed to generate MSG pool!\n");
    return NULL;
  }
  int d_sum=0;

  for(i=0;i<nSC;i++) d_sum+=dim[i];


  //the total size in sending ulong128 is the sum of length of the vectors, which is for Ra, and nSC long128 integers, which is for ra
  size_t R_size=sizeof(ulong128)*(d_sum);
  size_t r_size=sizeof(ulong128)*(nSC);
  //allocate memory for msg
  for(i=0;i<2;i++)
  {
    m->msg_to[i]=(ulong128*)malloc(R_size);//"R"
    m->r[i]=(ulong128*)malloc(r_size);//"r".
    if(!m->msg_to[i])
    {
      fprintf(stderr,"Failed to allocate %ld memory for Ra/Rb.\n",R_size);
      for(j=0;j<i;j++) { free(m->msg_to[j]); free(m->r[j]);}
      return NULL;
    }
    else if(!m->r[i])
    {
      fprintf(stderr,"Failed to allocate %ld memory for ra/rb.\n",r_size);
      free(m->msg_to[i]);
      for(j=0;j<i;j++) { free(m->msg_to[j]); free(m->r[j]);}
      return NULL;
    }

  }

  int offset=0;

  for(int k=0;k<nSC;k++)
  {
    int local_dim=dim[k];//the length of the k-th vector
    int local_size=dim[k]*sizeof(ulong128);//the length of the k-th vector

    //copy Ra and Rb to msg
    ulong128* p1=p->data[0];
    ulong128* p2=p->data[1];
    for(int i=0;i<local_dim;i++)
    {
      m->msg_to[0][i+offset]=p1[i];
      m->msg_to[1][i+offset]=p2[i];
    }
    offset+=(local_dim);
  }

  //reset offset to set ra/rb
  offset=0;
  for(int k=0;k<nSC;k++)
  {
    int local_dim=dim[k];//the length of the k-th vector
    for(i=2;i<4;i++)
    {
      q[i-2]=p->data[i][local_dim-1];
      if(q[i-2] == 0) rnd[i-2] = 0;
      else
        rnd[i-2]=random()%q[i-2];//smaller than q[i]
    }

    m->r[0][k]=(q[0]-rnd[0])*(Q)+(q[1]-rnd[1]);//ra
    m->r[1][k]=rnd[0]*(Q)+rnd[1];//rb
  }

  return m;
}

void run_server(int party_socket[2],POOL* p)
{
  bool getconnection=true;

  while(getconnection)
  {
    //get dimension sequence, which is separated by commas
    int nSC=0;
    //get the number of vectors
    if(protocol_readint(party_socket[AIDX],&nSC)<=0) { printf("failed to read how many times the scalar product protocol is going to conduct.\n"); break;}
#ifdef DEBUG
    printf("Run scalar prodduct protocol for %d times\n",nSC);
#endif
    if(nSC==0)
    {
      printf("End of service request. close connection!\n");
      getconnection=false;
      break;
    }

    int* dim=(int*)malloc(sizeof(int)*nSC);
    if(!dim) { fprintf(stderr,"cannot init dim array to get dimensions\n"); return;}

    if(protocol_read(party_socket[AIDX],(char*)dim,sizeof(int)*nSC)<0) break;
#ifdef DEBUG    
    printf("Read # of dims from Alice:\n");
    for(int n=0;n<nSC;n++) printf("%d ",dim[n]);
    printf("\n");
#endif
    int d_sum=0;//The sum of dimensions of all vectors
    for(int i=0;i<nSC;i++) d_sum+=dim[i];

    //generate and check if msgs are successfully generated
    MSG* m=set_multi_msg(p,dim,nSC);
    free(dim);
    if(!m)
    {
      fprintf(stderr,"Failed to generate random data!\n");
      for(int j=0;j<2;j++) 
      {
        protocol_writeint(party_socket[j],0);
        close(party_socket[j]);
      }
      fprintf(stderr,"can not init data pool!\n"); 
      getconnection=false;
      break;
    }

    //Bob first
    //the total size in sending ulong128 is the sum of length of the vectors, which is for Ra, and nSC long128 integers, which is for ra (or rb)
    //int d=sizeof(ulong128)*(d_sum+nSC);
    int d1=sizeof(ulong128)*(d_sum);
    int d2=sizeof(ulong128)*(nSC);

    //write Rb to Bob using socket party_socket[BIDX]
    protocol_write(party_socket[BIDX],(char*)(m->msg_to[1]),d1);
#ifdef DEBUG
    //debug
    fprintf(stderr,"write to Bob: \n");
    for(int i=0;i<d_sum;i++)
    {
      fprintf(stderr,"%lld",m->msg_to[1][i]);
      if(i+1<d_sum) fprintf(stderr,",");
    }
    fprintf(stderr,"\n");
    //end debug
#endif
    //write Ra to Alice using socket party_socket[AIDX]
    protocol_write(party_socket[AIDX],(char*)(m->msg_to[0]),d1);
#ifdef DEBUG
    //debug
    fprintf(stderr,"write to Alice:\n");
    for(int i=0;i<d_sum;i++)
    {
      fprintf(stderr,"%lld",m->msg_to[0][i]);
      if(i+1<d_sum) fprintf(stderr,",");
    }
    fprintf(stderr,"\n");
    //end debug
#endif
    //write rb to Bob using socket party_socket[BIDX]		
    protocol_write(party_socket[BIDX],(char*)(m->r[1]),d2);

#ifdef DEBUG
    //debug
    fprintf(stderr,"rb\n");
    for(int i=0;i<nSC;i++)
    {
      fprintf(stderr,"%lld",m->r[1][i]);
      if(i+1<nSC) fprintf(stderr,",");
    }
    fprintf(stderr,"\n");
    //end debug
#endif

    //write ra to Alice using socket party_socket[AIDX]
    protocol_write(party_socket[AIDX],(char*)(m->r[0]),d2);

#ifdef DEBUG
    //debug
    fprintf(stderr,"ra\n");
    for(int i=0;i<nSC;i++)
    {
      fprintf(stderr,"%lld",m->r[0][i]);
      if(i+1<nSC) fprintf(stderr,",");
    }
    fprintf(stderr,"\n");
    //end debug
#endif
    free_msg(m);

  }//end while
  /*
   */
  printf("End run_server on commodity server\n");
}

int main(int argc, char** argv)
{
  Protocol_ARG arg;
  int socket;
  int sid=0;
  int party_socket[2];
  int dim=0;
  int go=255;
  int idx=0;
  int i=0;
  int j=0;
  int partyid[2];

  void* t_status;

  MSG* m=NULL;
  POOL* p=NULL;
  //TH_RND_MSG t_rnd_data[2];
  if(argc!=2) { usage(argv[0]); return 1; }
  if(parse_ProtocolARG(&arg,argv[1])<0) return -1;
  srandom(time(NULL));
  socket=setServer(arg.port[0]);
  //if(socket<0) return -1;

  printf("listen on port %d...\n",arg.port[0]);

  p=init_commodity_server("commodity.dat");
  if(!p) { fprintf(stderr,"Failed to init data pool!\n"); close(socket);return -2;}

  printf("successfully init commodity server\n");

  while(1)
  {
    int pid=-1;//process id;

    printf("waiting for new connection request\n");

    for(int i=0;i<2;i++) 
    {
      party_socket[i]=wait_connection(socket);//waiting for the first party
      if(party_socket[i]<0)
      { 
        fprintf(stderr,"socket id=%d. Can not create connection with the %d-th party.\n",party_socket[i],i);
        for(int j=0;j<i;j++) close(party_socket[j]);
        break;
      }
      else 
      {
        fprintf(stderr,"get connection request from the %d-th child, socket id=%d\n",i,party_socket[i]);
      }
    }

    pid=fork();
    if(pid<0)//failed to fork
    { fprintf(stderr,"failed to fork!\n"); return -1;}
    if(pid==0)//child. start a new service
    {
      run_server(party_socket,p);
      for(int i=0;i<2;i++) close(party_socket[i]);
      close(socket);//close the socket used to accept connection
      free_pool(p);
      return 0;
    }

    for(int i=0;i<2;i++) close(party_socket[i]);

  }//end while(1)

  return 0;
}

