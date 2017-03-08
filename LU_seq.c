#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static const enum pm { 
   NORM, L, U 
} PRINT_MODE;

void print_mx (size_t dim, float **M, enum pm mode);
float **gen_mx (size_t dim);

int main (int argc, char *argv[])
{
   srand(time(NULL));
   int mx_size = 0;
   if (argc < 2) {
      printf("Matrix size missing in the arguments\n");
      return EXIT_FAILURE;
   }
   mx_size = atol(argv[1]);
   float **A = gen_mx(mx_size);

   #ifdef ALU
   printf("[A]\n");
   print_mx(mx_size, A, NORM);
   #endif
   clock_t begin = clock();
   size_t i, j, k;
   for (i = 0; i < mx_size; i++) {
      for (j = i + 1; j < mx_size; j++) {
         float m = A[j][i] / A[i][i];
         for (k = i + 1; k < mx_size; k++) {
            A[j][k] -= m * A[i][k];
         }
         A[j][i] = m;
      }
   }
   clock_t end = clock();
   double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
   #ifdef ALU
   printf("\n\n[L]\n");
   print_mx(mx_size, A, L);
   printf("\n\n[U]\n");
   print_mx(mx_size, A, U);
   #endif
   printf("seq: %f s\n", time_spent);
   return 0;
}

float **gen_mx (size_t dim)
{
   int i, j;
   float **M = malloc(sizeof(float *) * dim);
   for (i = 0; i < dim; i++) {
      M[i] = malloc(sizeof(float) * dim);
      for (j = 0; j < dim; j++) {
         M[i][j] = rand() % 101 - 50;
      }
   }

   return M;
}

void print_mx (size_t dim, float **M, enum pm mode)
{
   size_t i, j;
   float z = 0, u = 1;
   for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
         if (mode == L && j >= i) {
            if (j == i) {
               printf("% *.*f\t", 4, 2, u);
            } else {
               printf("% *.*f\t", 4, 2, z);
            }
         } else if (mode == U && j < i) {
            printf("% *.*f\t", 4, 2, z);
         } else {
            printf("% *.*f\t", 4, 2, M[i][j]);
         }
      }
      putchar('\n');
   }
}
