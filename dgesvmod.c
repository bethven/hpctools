#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mkl_lapacke.h"

void *generate_matrix(int size);
void *spd_generate_matrix(int n, int *info);
void print_matrix(const char *name, double *matrix, int size);
int  check_result(double *bref, double *b, int size);
void my_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb);
void gauss(double *a, int size);

int main(int argc, char *argv[])
    {
        int size = atoi(argv[1]);
	int info;
        double *a, *aref;
        double *b, *bref;

        a = spd_generate_matrix(size, &info);
        aref = generate_matrix(size);        
        b = generate_matrix(size);
        bref = generate_matrix(size);

        memcpy( aref, a, sizeof(double) * size * size);
        memcpy( bref, b, sizeof(double) * size * size);


        
//      print_matrix("A", a, size);
//      print_matrix("B", b, size);

        // Using MKL to solve the system
        MKL_INT n = size, nrhs = size, lda = size, ldb = size;
        MKL_INT *ipiv = (MKL_INT *)malloc(sizeof(MKL_INT)*size);

        clock_t tStart = clock();
        info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
        printf("Time taken by MKL: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

        tStart = clock();    
        MKL_INT *ipiv2 = (MKL_INT *)malloc(sizeof(MKL_INT)*size);        
        my_dgesv(n, nrhs, a, lda, ipiv2, b, ldb);
        printf("Time taken by my implementation %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);


        if (check_result(bref,b,size)==1)
            printf("Result is ok!\n");
        else    
            printf("Result is wrong!\n");
        
//      print_matrix("X", b, size);
//      print_matrix("Xref", bref, size);
return 0;
    }

void *generate_matrix(int size)
{
    int i;
    double *matrix = (double *)malloc(sizeof(double) * size * size);
    srand(1);

    for (i = 0; i < size * size; i++)
    {
        matrix[i] = rand() % 100;
    }

    return matrix;
}

void *spd_generate_matrix(int n, int *info) {
    int i;
    double *matrix = generate_matrix(n);
    double *matrix_spd = (double *)malloc(sizeof(double) * n * n);
    for(int i=0; i<n; i++) {
       for(int j=0; j<n; j++) {
          double sum = 0;
	  for(int k=0; k<n; k++) {
	     sum += matrix[i*n+k] * matrix[j*n+k];
	  }
	  matrix_spd[i*n+j] = sum;
             
	}
      }
     free(matrix);
     return matrix_spd;
}




void print_matrix(const char *name, double *matrix, int size)
{
    int i, j;
   printf("matrix: %s \n", matrix);
 
    for (i = 0; i < size; i++)
    {
            for (j = 0; j < size; j++)
            {
                printf("%f ", matrix[i * size + j]);
            }
            printf("\n");
    }
}

int check_result(double *bref, double *b, int size) {
    int i;
    for(i=0;i<size*size;i++) {
       if ((bref[i]!=b[i])) return 1;

    }
    return 0;
}

void  my_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb) {
        gauss(a, n);
        LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, a, lda, ipiv, b, ldb);

}

void gauss(double *a, int size)
{
    int Col,C1,C2,A,cu;
    float p,v;

    for(Col=1;Col<=size;Col++){
        cu=0;
	A=Col;
        while(cu==0){
           if((a[A*size]>0.0000001)||((a[A*size]<-0.0000001))){
                cu=1;}
            else A++;}
        p=a[A*size];
        for(C1=1;C1<=(size+1);C1++){
            v=a[C1+Col*size];
            a[C1+Col*size]=a[C1+Col*size];
            a[C1+Col*size]=v/p;}
        for(C2=Col+1;C2<=size;C2++){
            v=a[C1+Col*size];
            for(C1=Col;C1<=(size+1);C1++){
                a[C1+Col*size]=a[C1+Col*size]-v*a[C1+Col*size];}
    }
}
    for(Col=size;Col>=1;Col--) for(C1=(Col-1);C1>=1;C1--){
        a[C1+Col*size+1]=a[C1+Col*size+1]-a[C1+Col*size]*a[C1+Col*size+1];
        a[C1+Col*size]=0;
    }
}


