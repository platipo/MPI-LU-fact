# MPI LU decomposition of square matrix

> A first approach to Open MPI with LU decomposition

This project was created to be as a case study for Open MPI and matrices.

I found [this slides](http://www.cse.buffalo.edu/faculty/miller/Courses/CSE633/Tummala-Spring-2014-CSE633.pdf) on LU decomposition fairly straightforward.

### Usage

``make`` arguments you can specify are matrix size with ``N`` and process number with ``P``. For
small matices, size less than 16, both sequential and parallel program will print the matrix.
The program will always output running time.

Example with mpi:

```
$ make mpi P=2 N=3
mpicc LU_mpi.c -DALU -o par.out
mpirun -np 2 ./par.out 3
[A]
 8.00	 18.00	 5.00	
 24.00	 10.00	 2.00	
-11.00	-45.00	-4.00	

[L]
 1.00	 0.00	 0.00	
 3.00	 1.00	 0.00	
-1.38	 0.46	 1.00	

[U]
 8.00	 18.00	 5.00	
 0.00	-44.00	-13.00	
 0.00	 0.00	 8.86	

mpi: 0.000252 s
```

### Part one: sequential

Sequential algorithm is implemented in `LU_seq.c` and below there is pseudocode:

```
FOR ii = 0,..,n-1
   FOR jj = ii+1,...,n-1
      m = a(jj,ii) / a(ii,ii)
      FOR kk = ii+1, ...,n-1
         a(jj,kk) -= m * a(ii,kk)
      END FOR
      a(jj,ii) = m
   END FOR
END FOR
```

### Part two: parallel

The parallel implementation was tougher than expected, in fact it took me two days only to understand the MPI inner working. I left [below](#Useful-resources) some links I found critical to deepen my MPI knowledge.

Some notes on matrix A:

* it's stored as a contiguous pointer and accessed as ``M[ii * mx_dim + jj]`` where ``ii`` is the row and ``jj`` the column.
* The matices L and U are both saved in A and in the end A is made up like
```
      u u u u u u u
      l u u u u u u
      l l u u u u u
  A = l l l u u u u
      l l l l u u u
      l l l l l u u
      l l l l l l u 
```

#### [v1](https://github.com/platipo/MPI-LU-fact/releases/tag/v1.0)

The idea that stands behind the algorithm is a master-slaves mechanism.

Indeed there is a coordinator, process 0, that splits work distribution among other processes and saves it.
Master divistribute the row j to the worker ``j % (proc num - 1) + 1``, that generates numbers from 1 to num proc because 0 coordinates only. 
After sending last row coordinator forward to each worker an "empty message" with end tag and copy recived rows in the original position (each row length is increased by one to store row number in the end).

On the other hand the worker receives load, decompose each row, "wraps" rows in a single message the work done and send it back to the master process when recives end tag message.

Of course this isn't the best implementation because master sends too many messages, one for each "raw" row.

#### [v2](https://github.com/platipo/MPI-LU-fact/releases/tag/v1.1)

Save idea behind [v1](https://github.com/platipo/MPI-LU-fact/releases/tag/v1.0) but in this implementation wrapping more rows together to lower messagges is a key strength.
As matter of fact master process evaluate for each column step which rows will be assigned to which worker and store their pointers. Therefore coordinator takes advantage of the previous abstraction to creates bigger messagges with multiple contiguous rows and send them to workers.

Slaves remain fundamentally unchanged.

Not sure, but non blocking calls seems to be a problem in really big matrices (need more testing).

#### [v3](https://github.com/platipo/MPI-LU-fact/releases/tag/v1.2)

Third implementation diverges from first two, expecially because there is no master, all processes are peer workers.
Every process starts with a clone of the original matrix and rows are decomposed by ``j % proc`` process. Finished the column every process updates the others by broadcasting their work row by row.

### Useful resources

* http://mpitutorial.com/
* https://github.com/wesleykendall/mpitutorial/tree/gh-pages/tutorials
* http://stackoverflow.com/questions/10017301/mpi-blocking-vs-non-blocking#10017425
* http://condor.cc.ku.edu/~grobe/docs/intro-MPI-C.shtml
* http://stackoverflow.com/questions/10490983/mpi-slave-processes-hang-when-there-is-no-more-work
* http://www.cs.indiana.edu/classes/b673/notes/mpi4.html
* http://stackoverflow.com/questions/21512975/what-is-the-difference-between-isend-and-issend-in-openmpi

