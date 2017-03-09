#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define ln() putchar('\n')
#define GENERIC_TAG (0)

typedef struct {
   size_t n_rows;
   size_t *cell;
} column_step;

float *gen_mx (size_t dim);
float *gen_row(size_t dim);
float *gen_row_ref (size_t dim, size_t ref);
void print_mx (float *M, size_t dim, size_t sep);
float *forw_elim(float *master_row, float *slave_row, size_t dim);
void U_print (float *M, int dim);
void L_print (float *M, int dim);

int main(int argc, char *argv[])
{
   const int root_p = 0;
   int mx_size = 0, p, id;
   if (argc < 2) {
      printf("Matrix size missing in the arguments\n");
      return EXIT_FAILURE;
   }
   mx_size = atol(argv[1]);
   MPI_Init(NULL, NULL);
   MPI_Comm_size(MPI_COMM_WORLD, &p);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   // we got some trouble with generating random numbers because processes
   // accessed at same time to rand function and it generated same number for all processes
   //srand(time(NULL) + id * p + rand());

   if (p < 2) {
      perror("Too few workers, minimum 2\n");
      MPI_Abort(MPI_COMM_WORLD, 0);
      return EXIT_FAILURE;
   }

   /*
    * map - link every row to a slave
    * save_point - memorize save point for each row
    */
   //size_t *map, *save_point;
   column_step **smap;
   float *A;
   int i;

   /*
    * Square matrix generator
    */
   if (id == root_p) {
      srand(time(NULL));
      A = gen_mx(mx_size);
      #ifdef ALU
      printf("[A]\n");
      print_mx(A, mx_size * mx_size, mx_size);
      ln();
      #endif
   }

   /*
    * LU factorization
    */
   if (id == root_p) {
      /*
       * To distriubure more efficiently the load, instead of computing segmented rows,
       * it's better to send continuos row block like
       * https://github.com/puneetar/Parallel-LU-Factorization-with-OpenMP-MPI/blob/master/MPI/MPI.c#L52
       */
      int tmp_size = mx_size - 1, p_no_root = p - 1, g = 0; // counter
      smap = malloc(sizeof(column_step *) * tmp_size);
      // all LU decomposiztion last mx_size * (mx_size - 1) / 2
      for (i = 0; i < tmp_size; i++) {
         int jp, kpoint, split_work = tmp_size - i, n_rows = 0;
         smap[i] = malloc(sizeof(column_step) * p_no_root);
         for (jp = 0; jp < p_no_root; jp++) {
            smap[i][jp].n_rows = split_work / (p_no_root - jp);
            if (smap[i][jp].n_rows > 0) {
               smap[i][jp].cell = malloc(sizeof(size_t) * smap[i][jp].n_rows);
               for (kpoint = 0; kpoint < smap[i][jp].n_rows; kpoint++, n_rows++) {
                  smap[i][jp].cell[kpoint] = (size_t) &A[mx_size * (i + n_rows + 1) + i];
                  //printf("p:%d %d/%d %f\n", jp+1, kpoint + 1, smap[i][jp].n_rows, *((float *) smap[i][jp].cell[kpoint]));
               }
            } else {
               smap[i][jp].cell = NULL;
            }
            split_work -= smap[i][jp].n_rows;
         }
      }
   }

   MPI_Barrier(MPI_COMM_WORLD);
   double start = MPI_Wtime();

   int g = 0; // counter
   for (i = mx_size; i > 1; i--) {
      int ld = i * sizeof(float);
      float *root_row;
      if (id == root_p) {
         // reference of diagonal
         root_row = &A[(mx_size - i) * mx_size + mx_size - i];
      } else {
         root_row = malloc(sizeof(float) * i);
      }

      MPI_Bcast(root_row, i, MPI_FLOAT, root_p, MPI_COMM_WORLD);

      /*
       * Worker
       * Every slave waits assigned rows, execute forward elimination between master and recived rows,
       * save the rows and wait end signal to send back work.
       */
      int slave_msg_size = 0;
      float *ready_messagge = NULL;
      if (id != root_p) {
         MPI_Status status;
         int workload;
         MPI_Recv(&workload, 1, MPI_INT, root_p, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
         if (workload > 0) {
            int worker_elem = i * workload;
            //printf("[%d] work %d\n", id, workload);
            float mailbox[worker_elem];
            MPI_Recv(&mailbox, worker_elem, MPI_FLOAT, root_p, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            //print_mx(&mailbox[0], workload * i, i);
            int h;
            for (h = 0; h < workload; h++) {
               float *tmp_forw_elim = forw_elim(root_row, &mailbox[h * i], i);
               memmove(&mailbox[h * i], tmp_forw_elim, i * sizeof(float));
               free(tmp_forw_elim);
            }
            MPI_Send(&mailbox, workload * i, MPI_FLOAT, root_p, GENERIC_TAG, MPI_COMM_WORLD);
         }
      }

      /*
       * Coordinator
       * At beginning master process send at each slave the mapped row. When it reaches last matrix
       * row it sends to all slaves and "empty message" with end tag waiting till they send back
       * two messagges: the former contains the number of the rows reduced and the latter
       * concatenated rows.
       */
      else {
         int ip;
         for (ip = 0; ip < p - 1; ip++, g++) {
            //printf("=>[%d] work %d\n", ip + 1, smap[mx_size - i][ip].n_rows);
            MPI_Send(&smap[mx_size - i][ip].n_rows, 1, MPI_INT, ip + 1, GENERIC_TAG, MPI_COMM_WORLD);
            if (smap[mx_size - i][ip].n_rows > 0) {
               int r, coord_elem = i * smap[mx_size - i][ip].n_rows;
               float wraped_rows[coord_elem];
               for (r = 0; r < smap[mx_size - i][ip].n_rows; r++) {
                  memmove(&wraped_rows[r * i], (float *) smap[mx_size - i][ip].cell[r], ld);
               }
               //print_mx(wraped_rows, smap[mx_size - i][ip].n_rows * i, i);
               MPI_Send(&wraped_rows, coord_elem, MPI_FLOAT, ip + 1, GENERIC_TAG, MPI_COMM_WORLD);
               MPI_Recv(&wraped_rows, coord_elem, MPI_FLOAT, ip + 1, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               //print_mx(wraped_rows, smap[mx_size - i][ip].n_rows * i, i);
               for (r = 0; r < smap[mx_size - i][ip].n_rows; r++) {
                  memmove((float *) smap[mx_size - i][ip].cell[r], &wraped_rows[r * i], ld);
               }
            }
         }
      }
      //MPI_Barrier(MPI_COMM_WORLD);
   }

   double end = MPI_Wtime();

   if (id == root_p) {
      int tmp_size = mx_size - 1, p_no_root = p - 1;
      for (i = 0; i < tmp_size; i++) {
         int jp;
         for (jp = 0; jp < p_no_root; jp++) {
            if (smap[i][jp].n_rows > 0) {
               free(smap[i][jp].cell);
            }
         }
         free(smap[i]);
      }
      free(smap);

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
      free(A);
      printf("mpi: %f s\n", end - start);
   }

   MPI_Finalize();
   return EXIT_SUCCESS;
}

/*
 * gen_row - method used before gen_row_ref
 * generate random array of dim elements
 *
 * @dim dimension of array
 * @return random array (matrix row)
 */
float *gen_row (size_t dim)
{
   int i;
   float *row = malloc(sizeof(float) * dim);
   for (i = 0; i < dim; i++) {
      row[i] = rand() % 101 - 50;
   }

   return row;
}

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
 * gen_row_ref - similar to gen_row
 * only difference is that generated array
 * is dim + 1 big because array[0] is the
 * row reference
 *
 * @dim number of random array elements
 * @ref matrix row number
 * @return  random array (matrix row) with row reference
 */
float *gen_row_ref (size_t dim, size_t ref)
{
   int i;
   float *row = malloc(sizeof(float) * dim);
   row[0] = ref;
   for (i = 1; i < dim + 1; i++) {
      row[i] = rand() % 20 - 10;
   }

   return row;
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
 * forw_elim - forward Gauss elimination between mster and slave rows of
 * dim size
 *
 * @master_row row sent from master
 * @slave_row row sent from slave
 * @return reduced row
 */
float *forw_elim(float *master_row, float *slave_row, size_t dim)
{
   int i;
   // alloc +1 to store later the index
   float *reduc_row = malloc(sizeof(float) * (1 + dim));
   float l_coeff = reduc_row[0] = slave_row[0] / master_row[0];
   /*
      printf("master: ");
      print_mx(master_row, dim, dim);
      printf("slave: ");
      print_mx(slave_row, dim, dim);
      */
   for (i = 1; i < dim; i++) {
      reduc_row[i] = slave_row[i] - master_row[i] * l_coeff;
   }

   return reduc_row;
}

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
