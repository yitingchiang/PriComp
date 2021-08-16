#include <cstdio>
#include <sys/socket.h>
#include <netinet/tcp.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <cerrno>
#include <cstring>
#include <unistd.h>
#include <sys/uio.h>
#include <cstdlib>
#include "communicate.h"

int setServer(int portno)
{
  struct sockaddr_in serv_addr;

  int sockfd = socket(AF_INET, SOCK_STREAM, 0);
  if (sockfd < 0) { perror("ERROR opening socket"); return -1; }
  int v=1;
  socklen_t optlen=sizeof(v);

  //set reusable socket
  if(setsockopt(sockfd,SOL_SOCKET,SO_REUSEADDR,(void*)&v,optlen)<0) { perror("Failed to set reusable socket!"); return -1; }
  //no TCP delay
  if(setsockopt(sockfd,6,TCP_NODELAY,(void*)&v,optlen) < 0)
    printf("Cannot set TCP_NODELAY option on listen socket (%s)\n", strerror(errno));

  bzero((char *) &serv_addr, sizeof(serv_addr));
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = INADDR_ANY;//accept all
  serv_addr.sin_port = htons(portno);
  if (bind(sockfd, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0)
  {
    perror("ERROR on bind()!!");
    return -1;
  }

  if(listen(sockfd,5)<0)
  {
    perror("ERROR on listen()!");
    return -1;
  }

  return sockfd;
}

int wait_connection(int socket)
{
  socklen_t clilen;
  struct sockaddr_in cli_addr;
  clilen = sizeof(cli_addr);
  int newsocket = accept(socket,(struct sockaddr *) &cli_addr, &clilen);//wait for the cooodinator to call for start
  if (newsocket < 0)
  {
    perror("Error on accept coordinator's requirement.");
    return -1;
  }

  int v=1;	
  socklen_t optlen=sizeof(v);
  if(setsockopt(newsocket,6,TCP_NODELAY,(void*)&v,optlen) < 0)
    printf("Cannot set TCP_NODELAY option on listen socket (%s)\n", strerror(errno));

  return newsocket;
}

//for debugging
#ifdef DEBUG
void printBits(size_t const size, void const* const ptr)
{
  unsigned char *b = (unsigned char*) ptr;
  unsigned char byte;
  int i, j;

  for (i=size-1;i>=0;i--)
  {
    for (j=7;j>=0;j--)
    {
      byte = b[i] & (1<<j);
      byte >>= j;
      fprintf(stderr,"%u", byte);
    }
  }
  fprintf(stderr,"\n");
  //puts("");
}
#endif

int setClient(char* s_host,int s_port)
{
  int c_socket=socket(AF_INET,SOCK_STREAM, 0);
  if(c_socket<0) { perror("Error in setting client socket!"); return -1; }

  int v=1;
  socklen_t optlen=sizeof(v);
  if(setsockopt(c_socket,SOL_SOCKET,SO_REUSEADDR,(void*)&v,optlen)<0) { perror("Failed to set reusable socket!"); return -1; }

  sockaddr_in serverAddr;
  bzero((char *)&serverAddr,sizeof(serverAddr));
  serverAddr.sin_family=AF_INET;
  serverAddr.sin_port=htons(s_port);
  serverAddr.sin_addr.s_addr=inet_addr(s_host);
  if(serverAddr.sin_addr.s_addr==INADDR_NONE)
  {
    struct hostent *host=gethostbyname(s_host);
    if (host==NULL)
    {
      char errmsg[256];
      sprintf(errmsg,"Unable to resolve server: %s\n", s_host);
      perror(errmsg);
      return -1;
    }
    memcpy(&(serverAddr.sin_addr),host->h_addr,host->h_length);
  }
  if(connect(c_socket,(struct sockaddr*)&serverAddr,sizeof(serverAddr))<0)
  {
    char errmsg[256];
    sprintf(errmsg,"Unable to connect to server %s:%d\n",s_host,s_port);
    perror(errmsg);
    return -1;
  }

  if(setsockopt(c_socket,6,TCP_NODELAY,(void*)&v,optlen) < 0)
    fprintf(stderr,"Cannot set TCP_NODELAY option on listen socket (%s)\n", strerror(errno));


  return c_socket;
}

int protocol_read(int socket,char* msg,int len)
{
  int v=0;
  int s=0;

  int counter=0;
  while(v<len)
  {
    int size=len-v;

    s=recv(socket,msg+v,size,0);
    if(s<0) { perror("Error when read from socket!"); return -1;}
    v+=s;
  }

  return v;
}

//return read length
int protocol_readint(int socket,int* i)
{
  int v=recv(socket,i,sizeof(int),0);
  if(v<0) { perror("Error when read int from socket!"); }
  *i=ntohl(*i);
  return v;
}

int protocol_write(int socket,char* msg,int len)
{
  char* ptr=msg;
  int v=0;

  v=send(socket,ptr,len,0);

  if(v<0)
  {
    switch(errno)
    {
      case EBADF:
        fprintf(stderr,"1\n");
        fflush(stderr);
        break;
      case EACCES:
        fprintf(stderr,"2\n");
        fflush(stderr);
        break;
      case ENOTSOCK:
        fprintf(stderr,"3\n");
        fflush(stderr);
        break;
      case EFAULT:
        fprintf(stderr,"4\n");
        fflush(stderr);
        break;
      case EMSGSIZE:
        fprintf(stderr,"5\n");
        fflush(stderr);
        break;
      case EAGAIN:
        fprintf(stderr,"6\n");
        fflush(stderr);
        break;
      case ENOBUFS:
        fprintf(stderr,"7\n");
        fflush(stderr);
        break;
      case EHOSTUNREACH:
        fprintf(stderr,"8\n");
        fflush(stderr);
        break;
      case EISCONN:
        fprintf(stderr,"9\n");
        fflush(stderr);
        break;
      case ECONNREFUSED:
        fprintf(stderr,"10\n");
        fflush(stderr);
        break;
      case EHOSTDOWN:
        fprintf(stderr,"11\n");
        fflush(stderr);
        break;
      case ENETDOWN:
        fprintf(stderr,"12\n");
        fflush(stderr);
        break;
      case EADDRNOTAVAIL:
        fprintf(stderr,"13\n");
        fflush(stderr);
        break;
      case EPIPE:
        fprintf(stderr,"14\n");
        fflush(stderr);
        break;
      default:
        fprintf(stderr,"OTHER ERROR!\n");
        fflush(stderr);

    }
  }
  return v;
}

//return write length
int protocol_writeint(int socket,int i)
{
  int tmp=htonl(i);
  int v=send(socket,&tmp,sizeof(int),0);
  if(v<0) { perror("Error when write int to socket!"); }
  return v;
}


int get_value(ulong128* v, int socket)
{
  int read_len=sizeof(ulong128);
  int len=0;

  if((len=protocol_read(socket,(char*)v,read_len))<0)
    fprintf(stderr,"Error in get_value!");

  return len;
}

//send a mpz value to peer
int send_value(ulong128* v,int socket)
{
  size_t size=0;
  int len=0;

  len=protocol_write(socket,(char*)v,sizeof(ulong128));

  return len;
}

/*
int readVec(char* vec,size_t size)
{
  int len=protocol_read(SP_p_socket,vec,size);
  return len;
}

int writeVec(char* vec,size_t size)
{
  int len=protocol_write(SP_p_socket,vec,size);
  return len;
}
*/
