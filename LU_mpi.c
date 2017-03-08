#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

#define ln() putchar('\n')

const static enum {
   GENERIC_TAG, END_TAG
} CUSTOM_MPI_TAG;

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
   size_t *map, *save_point;
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
      #endif
   }

   /*
    * LU factorization
    */
   if (id == root_p) {
      // all LU decomposiztion last mx_size * (mx_size - 1) / 2
      int steps = mx_size * (mx_size - 1) / 2;
      //printf("st %d\n", steps);
      map = malloc(sizeof(size_t) * steps);
      save_point = malloc(sizeof(size_t) * steps);

      // compute save_points and map
      int g = 0; // counter
      for (i = 0; i < mx_size; i++) {
         int j;
         for (j = i + 1; j < mx_size; j++, g++) {
            save_point[g] = (size_t) &A[j * mx_size + i];
            map[g] = g % (p - 1) + 1;
            //printf("%d  %d  %d  %f\n", j, i, g, A[j * mx_size + i]);
            //printf("%d  %d  %d  %f\n", j, i, g, *((float *) save_point[g]));
         }
      }
   }

   MPI_Barrier(MPI_COMM_WORLD);
   double start = MPI_Wtime();

   int j = 0;
   for (i = mx_size; i > 1; i--) {
      int row_len = i + 1;
      int ld = row_len * sizeof(float);
      float *root_row;
      if (id == root_p) {
         // reference of diagonal
         root_row = &A[(mx_size - i) * mx_size + mx_size - i];
      } else {
         root_row = malloc(sizeof(float) * i);
      }

      MPI_Bcast(root_row, i, MPI_FLOAT, root_p, MPI_COMM_WORLD);

      /*
       * slave
       * Every slave waits assigned rows, execute forward elimination between master and recived rows,
       * save the rows and wait end signal to send back work.
       */
      int slave_msg_size = 0;
      float *ready_messagge = NULL;
      while (id != root_p) {
         float slave_row[row_len];
         MPI_Status status;
         MPI_Recv(&slave_row, row_len, MPI_FLOAT, root_p, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
         if (status.MPI_TAG == END_TAG) {
            //printf("[%d] ", id);
            //print_mx(ready_messagge, row_len * slave_msg_size, row_len * slave_msg_size);
            MPI_Send(&slave_msg_size, 1, MPI_INT, root_p, END_TAG, MPI_COMM_WORLD);
            if (slave_msg_size > 0) {
               MPI_Send(ready_messagge, row_len * slave_msg_size, MPI_FLOAT, root_p, END_TAG, MPI_COMM_WORLD);
            }
            free(ready_messagge);
            free(root_row);
            //printf("Process %d exiting work loop.\n", id);
            break;
         } else {
            float *red_row = forw_elim(root_row, slave_row, i);
            red_row[row_len - 1] = slave_row[row_len - 1];
            //printf("[%d] ", id);
            //print_mx(red_row, row_len, row_len);
            if (!slave_msg_size) {
               ready_messagge = malloc(ld);
            } else {
               ready_messagge = realloc(ready_messagge, (slave_msg_size + 1) * ld);
            }
            memmove(&ready_messagge[slave_msg_size * row_len], red_row, ld);
            slave_msg_size++;
            free(red_row);
         }
      }

      /*
       * MASTER
       * At beginning master process send at each slave the mapped row. When it reaches last matrix
       * row it sends to all slaves and "empty message" with end tag waiting till they send back
       * two messagges: the former contains the number of the rows reduced and the latter
       * concatenated rows.
       */
      if (id == root_p) {
         int h; 
         for (h = 0; h < i - 1; h++, j++) {
            float work_row[row_len];
            work_row[row_len - 1] = j;
            memcpy(work_row, (float *) save_point[j], i * sizeof(float));
            //printf("%d diag: %f\n", j, *root_row);
            //print_mx(work_row, row_len, row_len);
            MPI_Send(work_row, row_len, MPI_FLOAT, map[j], GENERIC_TAG, MPI_COMM_WORLD);
         }

         int msg_len, end, cpy;
         for (int end = 1; end < p; end++) {
            MPI_Send(A, row_len, MPI_FLOAT, end, END_TAG, MPI_COMM_WORLD);
            MPI_Recv(&msg_len, 1, MPI_INT, end, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            if (msg_len > 0) {
               int buff_size = row_len * msg_len;
               float *msg_buffer = malloc(sizeof(float) * buff_size);
               MPI_Recv(msg_buffer, buff_size, MPI_FLOAT, end, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
               //printf("[%d] ", i);
               //print_mx(msg_buffer, buff_size, buff_size);
               for (cpy = 0; cpy < msg_len; cpy++) {
                  int save_index = msg_buffer[(cpy * row_len) + (row_len - 1)];
                  memcpy((float *) save_point[save_index], &msg_buffer[row_len * cpy], i * sizeof(float));
               }
               free(msg_buffer);
            }
         }
      }
      //MPI_Barrier(MPI_COMM_WORLD);
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
   float *row = malloc(sizeof(float) * (dim + 1));
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
