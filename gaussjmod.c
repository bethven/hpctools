#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "mkl_lapacke.h"
#include <math.h>

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

void main(int argc, char *argv[])
{
  clock_t tStart = clock();
  int nc = atoi(argv[1]);
  int i, j, k, x, y;
  printf("Sistemas de Ecuaciones Gauss-Jordan\n\n");
  float matriz[nc][nc+1], ec[1][nc+1], ng[nc], bx0, bx1, bx2, ud;
  int size=nc;
  double *matrix, *aref;
  double *b, *bref;

  matrix=generate_matrix(nc);
// print_matrix("A", matrix, nc);

for (i = 0; i < size; i++)
{
     for (j = 0; j < size; j++)
     {
       matriz[i][j]=matrix[i * size + j];
     }

}

for(i=0;i<nc-1;i++)
{
  bx1 = matriz[i][i];
  for(j=i+1;j<nc;j++)
  {
       bx2 = matriz[j][i];
       for(k=i;k<nc+1;k++)
      {
         bx0 = matriz[i][k];
         ec[0][k] = bx0;
         bx0 = bx0*bx2*(-1);
         matriz[i][k] = bx0;
         bx0 = matriz[j][k];
         bx0 = bx0*bx1*1;
         matriz[j][k] = bx0;
      }
      for(y=0;y<nc+1;y++)
      {
         bx0 = matriz[i][y] + matriz[j][y];
         matriz[j][y] = bx0;
       }
      for(x=i;x<nc+1;x++)
      {
         bx0 = ec[0][x];
         matriz[i][x] = bx0;
      }
  }
}

k=1;
x=0;
for(i=nc-1;i>=0;i--)
{
    bx0 = matriz[i][i];
    bx1 = matriz[i][i+k];
    if(k==1)
        {
            bx2 = bx1/bx0;
            ng[x] = bx2;
            for(j=nc-1;j>=0;j--)
            {
                ud = matriz[j][nc-k];
                ud = ud*bx2;
                matriz[j][nc-k] = ud;
            }
            k++;
            x++;
        }
        else
        {
            for(y=i+1;y<nc;y++)
            {
                ud = matriz[i][y];
                bx1 = bx1-ud;
            }
            bx2 = bx1/bx0;
            ng[x] = bx2;
            for(j=i;j>=0;j--)
            {
                ud = matriz[j][nc-k];
                ud = ud*bx2;
                matriz[j][nc-k] = ud;
            }
            k++;
            x++;
        }
}
printf("\n FIN EJECUCION:  ");
printf("Time taken by my implementation: %.6fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

}
