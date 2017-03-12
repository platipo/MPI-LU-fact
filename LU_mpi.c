#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define ln() putchar('\n')
#define GENERIC_TAG (0)

float *gen_mx (size_t dim);
float *gen_row(size_t dim);
float *gen_row_ref (size_t dim, size_t ref);
void print_mx (float *M, size_t dim, size_t sep);
void forw_elim(float **origin, float *master_row, size_t dim);
void U_print (float *M, int dim);
void L_print (float *M, int dim);

int main(int argc, char *argv[])
{
   srand(time(NULL));
   const int root_p = 0;
   int mx_size = 0, p, id;
   if (argc < 2) {
      printf("Matrix size missing in the arguments\n");
      return EXIT_FAILURE;
   }
   mx_size = atol(argv[1]);
   float *A = gen_mx(mx_size);

   MPI_Init(NULL, NULL);
   MPI_Comm_size(MPI_COMM_WORLD, &p);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   if (id == root_p) {
      #ifdef ALU
      printf("[A]\n");
      print_mx(A, mx_size * mx_size, mx_size);
      ln();
      #endif
   }

   int i, j, tmp_size = mx_size - 1, diag_ref = 0;
   double start = MPI_Wtime();

   for (i = 0; i < tmp_size; i++, diag_ref++) {
      float *diag_row = &A[diag_ref * mx_size + diag_ref];
      for (j = diag_ref + 1; j < mx_size; j++) {
         if (j % p == id) {
            float *save = &A[j * mx_size + diag_ref];
            //printf("[%d] ", id);
            //print_mx(save, mx_size - diag_ref, mx_size - diag_ref);
            forw_elim(&save, diag_row, mx_size - diag_ref);
         }
      }
      //MPI_Barrier(MPI_COMM_WORLD);

      for (j = diag_ref + 1; j < mx_size; j++) {
         float *save = &A[j * mx_size + diag_ref];
         MPI_Bcast(save, mx_size - diag_ref, MPI_FLOAT, j % p, MPI_COMM_WORLD);
      }
      //MPI_Barrier(MPI_COMM_WORLD);

      /*
      if (id == root_p) {
         #ifdef ALU
         printf("(%d)\n", i);
         print_mx(A, mx_size * mx_size, mx_size);
         ln();
         #endif
      }
      */
   }

   double end = MPI_Wtime();

   if (id == root_p) {
      /*
      printf("[LU]\n");
      print_mx(A, mx_size * mx_size, mx_size);
      */
      #ifdef ALU
      printf("\n[L]\n");
      L_print(A, mx_size);
      printf("\n[U]\n");
      U_print(A, mx_size);
      #endif
      printf("mpi: %f s\n", end - start);
   }
   free(A);


   MPI_Finalize();
   return EXIT_SUCCESS;
}

/*
 * gen_mx - generate contiguous matrix
 *
 * @dim dim x dim matrix
 * @return matrix
 */
float *gen_mx (size_t dim)
{
   int i, j, tot = dim * dim;
   float *M = malloc(sizeof(float) * tot);
   for (i = 0; i < tot; i++) {
      M[i] = rand() % 101 - 50;
   }

   return M;
}

/*
 * mx_print - dumb matrix print function
 *
 * @M matrix/row
 * @dim matrix/row dimension
 * @sep where put separator
 */
void print_mx (float *M, size_t dim, size_t sep)
{
   int i, j;
   for (i = 0; i < dim; i++) {
      printf("% *.*f\t", 4, 2, M[i]);
      if ((i + 1) % sep == 0) {
         ln();
      }
   }
}

/*
 * forw_elim - forward Gauss elimination
 *
 * @origin row pointer by reference
 * @master_row row in which lays diagonal
 */
void forw_elim(float **origin, float *master_row, size_t dim)
{
   if (**origin == 0)
      return;

   float k = **origin / master_row[0];

   int i;
   for (i = 1; i < dim; i++) {
      (*origin)[i] = (*origin)[i] - k * master_row[i];
   }
   **origin = k;
}

/*
 * U_print - dumb U matrix print function
 */
void U_print (float *M, int dim)
{
   int i, j;
   float z = 0;
   for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
         if (j >= i) {
            printf("% *.*f\t", 4, 2, M[i * dim + j]);
         } else {
            printf("% *.*f\t", 4, 2, z);
         }
      }
      ln();
   }
}

/*
 * L_print - dumb L matrix print function
 */
void L_print (float *M, int dim)
{
   int i, j;
   float z = 0, u = 1;
   for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
         if (j > i) {
            printf("% *.*f\t", 4, 2, z);
         } else if (i == j) {
            printf("% *.*f\t", 4, 2, u);
         } else {
            printf("% *.*f\t", 4, 2, M[i * dim + j]);
         }
      }
      ln();
   }
}
