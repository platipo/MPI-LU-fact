# MPI LU decomposition of square matrix

> A first approach to Open MPI with LU decomposition

This project was created to be as a case study for Open MPI and matrices.

I found [this slides](http://www.cse.buffalo.edu/faculty/miller/Courses/CSE633/Tummala-Spring-2014-CSE633.pdf) on LU decomposition fairly straightforward.

### Part one: sequential

Sequential algorithm is implemented in `LU_seq.c`
 and below there is pseudocode:

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

The parallel implementation was tougher than expected, in fact it took me two days only to understand the MPI inner working.

The idea that stands behind the algorithm is a master-slaves concept. 

Indeed there is a coordinator, process 0, that divides and distributes the rows to the workers, waits the forward elimination processing, gets back rows and copy them in the original position.

On the other hand the worker receives rows untill end signal arrives, decompose each row, "wraps" rows in a single message and send back to the master process.

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

### Useful resources

* http://mpitutorial.com/
* https://github.com/wesleykendall/mpitutorial/tree/gh-pages/tutorials
* http://stackoverflow.com/questions/10017301/mpi-blocking-vs-non-blocking#10017425
* http://condor.cc.ku.edu/~grobe/docs/intro-MPI-C.shtml
* http://stackoverflow.com/questions/10490983/mpi-slave-processes-hang-when-there-is-no-more-work

