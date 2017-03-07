CC=gcc
MPI=mpicc
PAR=LU_mpi.c
EX_PAR=par.out
SEQ=LU_seq.c
EX_SEQ=seq.out
SAMPLE=100

default: sequential mpi

sequential: $(SEQ)
	$(CC) $(SEQ) -o $(EX_SEQ)
	./$(EX_SEQ) $(SAMPLE)

mpi: $(PAR)
	$(MPI) $(PAR) -o $(EX_PAR)
	mpirun -np 5 ./$(EX_PAR) $(SAMPLE)

clean:
	$(RM) -f *.out
