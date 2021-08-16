#include "smc.h"
#include "two_parties/protocol/s_float.h"


/*
 K(X1,X2)=X1*X2, where * is the scalar product.
 Note that X1 is held by alice, and X2 is held by bob.
 X1 and X2 are length-len double arrays.
 Initially, these array are not pair-wised shared.
*/
SFloat Secure_Linear_Kernel(double* X1, double* X2,double param, int len,int client_type);

/*
 K(X1,X2)=e^{(-||X1-X2||^2)*delta}
 Note that X1 is held by alice, and X2 is held by bob.
 X1 and X2 are length-len double arrays that Alice and Bob hold their own arrays. Initially, these array are not pair-wised shared.
 delta is the parameter of RBF kernel.
*/
SFloat Secure_RBF(double* X1, double* X2,double delta, int len, int client_type);

//param is not used in this kernel
double Linear_Kernel(double* X1, double* X2, double param, int len);

//K(X1,X2)=e^{(-||X1-X2||^2)/(2*delta)}
double RBF(double* X1, double* X2,double param,int len);

SFloat* compute_K(double* X1,int len1,double* X2,int len2, int nCol, SFloat(*Secure_kernel)(double*,double*,double,int,int),double(*kernel)(double*,double*,double,int),double param,int client_type);
