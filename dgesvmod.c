#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "mkl_lapacke.h"

double *generate_matrix(int size)
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


void gauss(double *a, int size)
{
    int x,y,z,f,cu;
    float p,v;

    for(x=1;x<=size;x++){
        cu=0;
	f=x;
        while(cu==0){
           if((a[f*size]>0.0000001)||((a[f*size]<-0.0000001))){
                cu=1;}
            else f++;}
        p=a[f*size];
        for(y=1;y<=(size+1);y++){
            v=a[y+x*size];
            a[y+x*size]=a[y+x*size];
            a[y+x*size]=v/p;}
        for(z=x+1;z<=size;z++){
            v=a[y+x*size];
            for(y=x;y<=(size+1);y++){
                a[y+x*size]=a[y+x*size]-v*a[y+x*size];}
    }
}
        for(x=size;x>=1;x--) for(y=(x-1);y>=1;y--){
        a[y+x*size+1]=a[y+x*size+1]-a[y+x*size]*a[y+x*size+1];
        a[y+x*size]=0;
    }
}


int my_dgesv(int n, int nrhs, double *a, int lda, int *ipiv, double *b, int ldb) {
    clock_t t_ini, t_fin;
    double secs;
    t_ini = clock();
    gauss(a, n);
    t_fin = clock();
    secs = (double)(t_fin - t_ini) / CLOCKS_PER_SEC;
    printf("Time taken by my implementation: %.16g milisegundos\n", secs * 1000.0);
       
}

int main(int argc, char *argv[])
    {
        int size = atoi(argv[1]);
	int info;
        double *a, *aref;
        double *b, *bref;

        a = generate_matrix(size);
        aref = generate_matrix(size);        
        b = generate_matrix(size);
        bref = generate_matrix(size);
        
//      print_matrix("A", a, size);
//      print_matrix("B", b, size);

        // Using MKL to solve the system
        MKL_INT n = size, nrhs = size, lda = size, ldb = size;
        MKL_INT *ipiv = (MKL_INT *)malloc(sizeof(MKL_INT)*size);

        clock_t tStart = clock();
        info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
//      printf("Time taken by MKL: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

        tStart = clock();    
        MKL_INT *ipiv2 = (MKL_INT *)malloc(sizeof(MKL_INT)*size);        
        my_dgesv(n, nrhs, a, lda, ipiv2, b, ldb);
	
        if (check_result(bref,b,size)==1)
            printf("Result is ok!\n");
        else    
            printf("Result is wrong!\n");
        
//      print_matrix("X", b, size);
//      print_matrix("Xref", bref, size);
return 0;
    }
