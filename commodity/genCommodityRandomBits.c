#include<inttypes.h>

#include "commodity.h"

void usage(const char *fname)
{
    printf("Generate random bits and write to the file: commodity.dat\n");
    printf("commodity.dat is used for the commodity server. Random bits for both Alice and Bob are written in this file.\n");
    printf("usage: %s <dimension> <domain>\n", fname);
    printf("\t<dimension> is the max length of Ra/Rb.\n");
    printf("\t<domain> is the domain the integer operations perform on.\n\n");
    printf("The format of commodity.dat:\n");
    printf("There are four rows in commodity.dat. Row 1 and Row 2 are values in <domain>, separated by comma. The width of each value is 20 chars, totally <dimension> values.\n");
    printf("Values in Row 3 and Row 4 represents q and r that is the result of inner product of Row 1 and Row 2.\n");
    printf("Specifically, let q[i] and r[i] be respectively the i-th values of Row 3 and Row 4. A[i] and B[i] are the i-th values in Row 1 and Row 2. Then inner prodduct of (A[0],A[1],...,A[i]) and (B[0],B[1],...,B[i])=q[i]*Q+r[i]. Q is 2^16.\n");

}

int64_t inner_product(uint64_t* Ra,uint64_t* Rb,int dim,int64_t domain)
{
    uint64_t v=0;
    for(int i=0;i<dim;i++) v+=((Ra[i]*Rb[i])%domain);
    return v;
}

void write_seq(uint64_t* data,int dim,FILE* wp)
{
    for(int i=0;i<dim;i++)
    {
        fprintf(wp,"%20ld",data[i]);
        if(i+1<dim) fprintf(wp,",");
    }
    fprintf(wp,"\n");
}

int main(int argc, char** argv)
{
    if(argc!=3) { usage(argv[0]); return 1;}
    srandom(time(NULL));
    uint64_t domain=(uint64_t)strtoll(argv[2],NULL,10);
    int dim=(int)strtol(argv[1],NULL,10);
    uint64_t* Ra=(uint64_t*)malloc(sizeof(int64_t)*dim);
    uint64_t* Rb=(uint64_t*)malloc(sizeof(int64_t)*dim);
    uint64_t* q=(uint64_t*)malloc(sizeof(int64_t)*dim);
    uint64_t* r=(uint64_t*)malloc(sizeof(int64_t)*dim);

    FILE* cwp=fopen("commodity.dat","w");

    for(int i=0;i<dim;i++)
    {
        Ra[i]=random()%domain; Rb[i]=random()%domain;
    }

    uint64_t q0=0;
    uint64_t r0=0;
    for(int i=0;i<dim;i++)
    {
        uint64_t q1=Ra[i]*Rb[i]/Q;
        uint64_t r1=Ra[i]*Rb[i]%Q;
        q[i]=q1+q0;
        r[i]=r1+r0;
        while(r[i]>=Q)
        {
            q[i]++;
            r[i]-=Q;
        }
        q0=q[i]; r0=r[i];
    }
    
    write_seq(Ra,dim,cwp);
    write_seq(Rb,dim,cwp);
    write_seq(q,dim,cwp);
    write_seq(r,dim,cwp);
    
    fclose(cwp);
    free(Ra); free(Rb); free(q); free(r);

    return 0;
}
