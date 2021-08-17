/**
 Define functions to establish and end a connection. In addition, functions to send and receive data are also defined here.
 */
#include "smc.h"

#ifndef COMMUNICATE

#define COMMUNICATE

/**
 The server waits for the client's connection on socket
 @param socket The socket to wait
 @result The socket number for the connection or -1 if error occurs.
 */
int wait_connection(int socket);
/**
 Prepare network socket as a server which listens on a specific port.
 @param portno The port this socket will be used to wait client's connection
 @return The socket created or -1 if error occurs
 */
int setServer(int portno);
/**
 Set client socket and connect to the server.
 @param s_host Host IPv4 string
 @param s_port Host port to connect to
 @return The socket or -1 if failed
 */
int setClient(char* s_host,int s_port);

/**
 Read data from a given socket.
 @param socket The socket be used to read data
 @param msg The buffer to store the read data
 @param len The number of bytes to read.
 @return The numnber of bytes read from the socket or -1 if error occurs. 
 */
int protocol_read(int socket,char* msg,int len);
/**
 Write data to a given socket.
 @param socket The socket to be used to write data.
 @param msg The data to be sent.
 @param len The number of bytes to write.
 @return The number of bytes sent, or -1 if error occurs. 
 */
int protocol_write(int socket,char* msg,int len);

/**
 Write an integer using the socket.
 @param socket The socket to be used to write data.
 @param i The value to be written.
 @return The number of bytes (sizeof(int)) on successful write, or -1 if error occurs.
 */
int protocol_writeint(int socket,int i);
/**
 Read an integer using the socket.
 @param socket The socket to be used to read data.
 @param i The variable to store the integer value.
 @return The number of bytes (sizeof(int)) on successful read, or -1 if error occurs.
 */
int protocol_readint(int socket,int* i);

/**
 Get an 128-bits unsigned int using socket.
 @param v Where the value be stored in.
 @param socket The socket used to read the data
 @return The number of bytes read
 */
int get_value(ulong128* v, int socket);

/**
 Write an 128-bits unsigned int using socket.
 @param v The value to be written out.
 @param socket The socket used to write the data
 @return The number of bytes written
 */
int send_value(ulong128* v,int socket);

#endif
